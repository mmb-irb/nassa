from itertools import combinations_with_replacement

import pandas as pd
import numpy as np

from .base import Base
from ..loaders.sequence import load_sequence
from ..loaders.trajectory import load_serfile
from .utils.heatmaps import correlation_plot


class CoordinateCorrelation(Base):
    """
    Execution plan and methods for correlation analyses
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

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
        # iterate over trajectories
        corr_results = {}
        for traj, seq in enumerate(sequences):
            trajectory_series = {coord.lower(): data[coord][traj]
                                 for coord in data.keys()}
            correlations = self.iterate_trajectory(
                seq, trajectory_series)
            corr_results[traj] = correlations
        return corr_results

    def iterate_trajectory(self, sequence, coordinates):
        corrtype = {
            "shift": "linear",
            "slide": "linear",
            "rise": "linear",
            "tilt": "circular",
            "roll": "circular",
            "twist": "circular"}
        start = 3 + sequence.flanksize
        end = sequence.size - sequence.flanksize - sequence.baselen
        all_subunits = sequence.all_subunits[2:-2]
        results_df = {}
        for crd1, crd2 in combinations_with_replacement(
                coordinates.keys(),
                r=2):
            crd1_series = pd.DataFrame(data=coordinates[crd1][start:end]).T
            crd2_series = pd.DataFrame(data=coordinates[crd2][start:end]).T
            crd1_series.columns = all_subunits
            crd2_series.columns = all_subunits
            method = self.get_corr_method(corrtype[crd1], corrtype[crd2])
            corr_matrix = pd.DataFrame({
                col: crd1_series.corrwith(
                    crd2_series[col],
                    method=method
                ) for col in crd2_series.columns})
            crd1_inner_dict = results_df.get(crd1, {})
            crd1_inner_dict[crd2] = corr_matrix
            results_df[crd1] = crd1_inner_dict

        # build complete dataset
        coords = list(results_df.keys())
        for crd1, inner_dict in results_df.items():
            missing_coords = [
                crd for crd in coords if crd not in inner_dict.keys()]
            for crd2 in missing_coords:
                results_df[crd1][crd2] = results_df[crd2][crd1]
        dfs = {k: pd.concat(v) for k, v in results_df.items()}
        complete_df = pd.concat(dfs, axis=1)
        results_df["complete"] = complete_df
        return results_df

    def get_corr_method(self, corrtype1, corrtype2):
        if corrtype1 == "circular" and corrtype2 == "linear":
            method = self.circlineal
        if corrtype1 == "linear" and corrtype2 == "circular":
            method = self.circlineal
        elif corrtype1 == "circular" and corrtype2 == "circular":
            method = self.circular
        else:
            method = "pearson"
        return method

    @staticmethod
    def circular(x1, x2):
        x1 = x1 * np.pi / 180
        x2 = x2 * np.pi / 180
        diff_1 = np.sin(x1 - x1.mean())
        diff_2 = np.sin(x2 - x2.mean())
        num = (diff_1 * diff_2).sum()
        den = np.sqrt((diff_1 ** 2).sum() * (diff_2 ** 2).sum())
        return num / den

    @staticmethod
    def circlineal(x1, x2):
        x2 = x2 * np.pi / 180
        rc = np.corrcoef(x1, np.cos(x2))[1, 0]
        rs = np.corrcoef(x1, np.sin(x2))[1, 0]
        rcs = np.corrcoef(np.sin(x2), np.cos(x2))[1, 0]
        num = (rc ** 2) + (rs ** 2) - 2 * rc * rs * rcs
        den = 1 - (rcs ** 2)
        correlation = np.sqrt(num / den)
        if np.corrcoef(x1, x2)[1, 0] < 0:
            correlation *= -1
        return correlation

    def make_tables(self, dataset):
        for traj, data in dataset.items():
            # make a new directory to save data for each trajectory
            save_subpath = self.save_path / f"traj_{traj}"
            save_subpath.mkdir(exist_ok=True)
            # data.to_csv(save_subpath / "coordinate_corr.csv")
            # coords = data.columns.get_level_values(0).unique()
            for crd1, inner_dict in data.items():
                if crd1 == "complete":
                    inner_dict.to_csv(save_subpath / "coordinate_corr.csv")
                else:
                    for crd2, df in inner_dict.items():
                        df.to_csv(save_subpath / f"{crd1}_{crd2}.csv")

    def make_plots(self, dataset):
        complete_dataset = pd.concat(
            [value['complete'] for value in dataset.values()],
            axis=1)
        save_subpath = self.save_path / "complete"
        save_subpath.mkdir(exist_ok=True)
        correlation_plot(
            complete_dataset,
            "coordcorr",
            save_subpath)
