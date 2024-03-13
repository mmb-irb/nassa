import pandas as pd
import numpy as np
import numpy.ma as ma 
from .base import Base
from .utils.heatmaps import arlequin_plot
from ..loaders.sequence import load_sequence
from ..loaders.trajectory import load_serfile
import numpy.ma as ma
from collections import Counter
class StiffnessDistributions(Base):

    def __init__(
            self,
            *args,
            **kwargs):
        super().__init__(
            *args,
            **kwargs)

    def extract(self):
        extracted = {}
        sequences = []
        for seq_file in self.sequence_files:
            seq = load_sequence(
                seq_file,
                unit_name=self.unit_name,
                unit_len=self.unit_len)
            sequences.append(seq)
        self.Ac = seq.Ac
        self.logger.debug(f"Adenine complement set to: <{self.Ac}>")
        extracted["sequences"] = sequences
        self.logger.info(f"loaded {len(sequences)} sequences")

        for coord in self.coordinate_info.keys():
            crd_data = []
            for ser_file in self.coordinate_info[coord]:
                ser = load_serfile(
                    ser_file,
                    self.tail,
                    self.n_lines)
                crd_data.append(ser)
            extracted[coord.lower()] = crd_data
            self.logger.info(
                f"loaded {len(crd_data)} files for coordinate <{coord}>")
        return extracted

    def transform(self, data):
        sequences = data.pop("sequences")
        results = {"stiffness": [], "covariances": {}, "constants": {}}
        for traj, seq in enumerate(sequences):
            traj_series = {coord.lower(): data[coord][traj]
                           for coord in data.keys()}
            traj_results = self.get_stiffness(
                seq,
                traj_series)
            results["stiffness"].append(traj_results["stiffness"])
            results["covariances"].update(traj_results["covariances"])
            results["constants"].update(traj_results["constants"])
        stiffness_df = pd.concat(results["stiffness"], axis=0)
        stiffness_df = stiffness_df.drop_duplicates(subset=[self.unit_name])
        stiffness_df = stiffness_df.set_index(self.unit_name)
        results["stiffness"] = stiffness_df
        return results

    def get_stiffness(
            self,
            sequence,
            series_dict):
        # get stiffness table for a given trajectory
        coordinates = list(series_dict.keys())
        results = {"stiffness": {}, "covariances": {}, "constants": {}}
        diagonals = {}
        start = 2 + sequence.flanksize
        end = sequence.size - (2 + sequence.baselen + sequence.flanksize - 1)
        for i in range(start, end):
            tetramer = sequence.get_subunit(i)
            ic_tetramer = sequence.inverse_complement(tetramer)
            cols_dict = {coord: series_dict[coord][i+1]
                         for coord in series_dict.keys()}
            stiffness_diag, cte, cov_df = self.get_subunit_stiffness(
                cols_dict,
                coordinates)
            diagonals[tetramer] = np.append(
                stiffness_diag,
                [np.product(stiffness_diag), np.sum(stiffness_diag)])
            diagonals[ic_tetramer] = np.append(
                stiffness_diag,
                [np.product(stiffness_diag), np.sum(stiffness_diag)])
            #results["covariances"][tetramer] = cov_df
            results["covariances"][ic_tetramer] = cov_df
            #results["constants"][tetramer] = cte
            results["constants"][ic_tetramer] = cte
        # build stiffness table
        columns = [sequence.unit_name] + coordinates + ["product", "sum"]
        results["stiffness"] = pd.DataFrame.from_dict(
            diagonals,
            orient="index").reset_index()
        results["stiffness"].columns = columns
        return results

    def get_subunit_stiffness(
            self,
            cols_dict,
            coordinates,
            scaling=[1, 1, 1, 10.6, 10.6, 10.6],
            KT=0.592186827):
        if (self.unit_len % 2) == 0:
            SH_av = cols_dict["shift"].mean()
            SL_av = cols_dict["slide"].mean()
            RS_av = cols_dict["rise"].mean()
            TL_av = self.circ_avg(cols_dict["tilt"])
            RL_av = self.circ_avg(cols_dict["roll"])
            TW_av = self.circ_avg(cols_dict["twist"])
        elif (self.unit_len % 2) == 1:
            SH_av = cols_dict["shear"].mean()
            SL_av = cols_dict["stretch"].mean()
            RS_av = cols_dict["stagger"].mean()
            TL_av = self.circ_avg(cols_dict["buckle"])
            RL_av = self.circ_avg(cols_dict["propel"])
            TW_av = self.circ_avg(cols_dict["opening"])
        cols_arr = [cols_dict[coord] for coord in coordinates]
        cols_arr = np.array(cols_arr).T

        cv = ma.cov(ma.masked_invalid(cols_arr), rowvar=False)
        cv.filled(np.nan)

        cov_df = pd.DataFrame(cv, columns=coordinates, index=coordinates)
        stiff = np.linalg.inv(cv) * KT
        last_row = [SH_av, SL_av, RS_av, TL_av, RL_av, TW_av]
        stiff = np.append(stiff, last_row).reshape(7, 6)
        stiff = stiff.round(6)
        stiff_diag = np.diagonal(stiff) * np.array(scaling)

        cte = pd.DataFrame(stiff)
        cte.columns = coordinates
        cte.index = coordinates + ["avg"]
        return stiff_diag, cte, cov_df

    @ staticmethod
    def circ_avg(xarr, degrees=True):
        n = len(xarr)
        if degrees:
            # convert to radians
            xarr = xarr * np.pi / 180
        av = np.arctan2(
            (np.sum(np.sin(xarr)))/n,
            (np.sum(np.cos(xarr)))/n) * 180 / np.pi
        return av

    def unimod_labels(self, df):
        # get limiting values for col1 and col2
        global_mean = df.mean(axis=0)
        global_std = df.std(axis=0)
        l1 = global_mean + global_std
        l2 = global_mean - global_std

        # 1: above mean+std (blue)
        # 0: between mean-std and mean+std (white)
        # -1: below mean-std (red)
        labeled_df = (df < l2) * -1 + (df > l1)
        labeled_df.loc["g_mean"] = global_mean
        labeled_df.loc["g_std"] = global_std

        return labeled_df

    def make_tables(self, dataset):
        # stiffness
        stiffness_data = dataset["stiffness"]
        stiffness_data.to_csv(self.save_path / "stiffness.csv")
        # covariances
        covariances_path = self.save_path / "covariances"
        covariances_path.mkdir(exist_ok=True)
        for key, val in dataset["covariances"].items():
            val.to_csv(covariances_path / f"{key}.csv")
        # constants
        constants_path = self.save_path / "constants"
        constants_path.mkdir(exist_ok=True)
        for key, val in dataset["constants"].items():
            val.to_csv(constants_path / f"{key}.csv")

    def make_plots(self, dataset):
        stiffness_data = dataset["stiffness"]
        labeled_df = self.unimod_labels(stiffness_data)
        print(stiffness_data)
        for col in labeled_df.columns:
            df = labeled_df[col]
            g_mean = df.loc["g_mean"]
            g_std = df.loc["g_std"]
            df = df.iloc[:-2]
            df = df.rename("col1")
            df = df.reset_index()
            df["col2"] = df["col1"].copy()
            arlequin_plot(
                df,
                g_mean,
                g_std,
                col,
                self.save_path,
                base=self.Ac,
                unit_name=self.unit_name,
                unit_len=self.unit_len)
