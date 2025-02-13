import pandas as pd
import numpy as np

from .base import Base
from .utils.bibitransformer import BiBiTransformer
from .utils.heatmaps import arlequin_plot
from ..loaders.sequence import load_sequence
from ..loaders.trajectory import load_serfile


class CoordinateDistributions(Base):
    """
    Execution plan and methods for coordinate distributions analysis pipeline
    """

    def __init__(
            self,
            max_iter=400,
            tol=1e-5,
            **kwargs):
        super().__init__(**kwargs)
        self.max_iter = max_iter
        self.tol = tol
        self.logger.debug(f"max_iter: {self.max_iter}")
        self.logger.debug(f"tol: {self.tol}")

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
        subunits_classification = {}
        # get dataframe for each coordinate
        for coordinate, coord_dataset in data.items():
            coordinate_df = self.coordinate_iteration(
                sequences,
                coordinate,
                coord_dataset)
            coordinate_df = coordinate_df.drop_duplicates(
                subset=[self.unit_name])
            coordinate_df, global_mean, global_std = self.add_modality_labels(
                coordinate_df)
            subunits_classification[coordinate] = [
                coordinate_df, global_mean, global_std]
        return subunits_classification

    def coordinate_iteration(
            self,
            sequences,
            coordinate,
            coord_dataset):
        coordinate_df = []
        for seq, dataseries in zip(sequences, coord_dataset):
            trajectory_df = self.trajectory_iteration(
                dataseries,
                seq,
                coordinate)
            coordinate_df.append(trajectory_df)
        # concatenate all dataframes
        coordinate_df = pd.concat(
            coordinate_df,
            axis=0)
        return coordinate_df

    def trajectory_iteration(
            self,
            dataseries,
            sequence,
            coordinate):
        trajectory_info = []
        # iterate over subunits
        start = 2 + sequence.flanksize
        end = sequence.size - (2 + sequence.baselen + sequence.flanksize - 1)
        for idx in range(start, end):
            # get unit and inverse-complement unit
            subunit = sequence.get_subunit(idx)
            ic_subunit = sequence.inverse_complement(subunit)
            self.logger.info(
                f"analyzing {self.unit_name} {subunit}/{ic_subunit}...")
            # add 1 to idx since .ser table includes an index
            ser = dataseries[idx + 1]
            ser = ser[~np.isnan(ser)].to_numpy()
            if ser.shape[0] < 2:
                self.logger.info("skipping because of insufficient data!")
                subunit_information = dict(
                    coordinate=coordinate,
                    binormal=False,
                    uninormal=True,
                    insuf_ev=False,
                    unimodal=True,
                    bics=[np.nan, np.nan],
                    mean1=np.nan,
                    mean2=np.nan,
                    var1=np.nan,
                    var2=np.nan,
                    w1=np.nan,
                    w2=np.nan)
                subunit_information[self.unit_name] = subunit
            else:
                # reshape dataset
                ser = ser.reshape(-1, 1)
                subunit_information = self.subunit_iteration(
                    ser,
                    subunit,
                    coordinate)
            if not self.bimod:
                subunit_information["unimodal"] = True
            trajectory_info.append(subunit_information)
            # add inverse-complement
            ic_subunit_information = subunit_information.copy()
            if ic_subunit_information["coordinate"] == 'shift' or ic_subunit_information["coordinate"] == 'tilt':
                ic_subunit_information["mean1"] = ic_subunit_information["mean1"]*-1
                ic_subunit_information["mean2"] = ic_subunit_information["mean2"]*-1
            ic_subunit_information[self.unit_name] = ic_subunit
            trajectory_info.append(ic_subunit_information)
        # create dataframe from list of dictionaries
        trajectory_df = pd.DataFrame.from_dict(trajectory_info)
        return trajectory_df

    def subunit_iteration(
            self,
            ser,
            subunit,
            coordinate):
        # model subunit's series
        bibi = BiBiTransformer(max_iter=self.max_iter, tol=self.tol)
        subunit_info = bibi.fit_transform(ser)
        # add subunit name and coordinate to info dictionary
        subunit_info[self.unit_name] = subunit
        subunit_info["coordinate"] = coordinate
        return subunit_info

    def add_modality_labels(self, data):
        df = data.set_index(self.unit_name)

        # Series with global mean value for each tetramer
        weighted_mean = df.apply(
            lambda t: t["mean1"] if (
                t["uninormal"] or t["insuf_ev"] or np.isnan(t["mean2"])
            ) else (t["mean1"]*t["w1"]+t["mean2"]*t["w2"]),
            axis=1)

        # Series with comparison values
        comparison_1 = df.apply(
            lambda t: t["mean1"] if (
                t["uninormal"] or t["insuf_ev"] or not t["unimodal"]
            ) else t["mean1"]*t["w1"]+t["mean2"]*t["w2"],
            axis=1)

        # mean1: uninormal+unimodal
        # u1*w1+u2*w2: binormal+unimodal
        # mean2: binormal+bimodal
        comparison_2 = df.apply(
            lambda t: np.nan if t["unimodal"] else t["mean2"],
            axis=1)
        comparison_2 = comparison_2.fillna(comparison_1)

        # get limiting values for col1 and col2
        global_mean = weighted_mean.mean()
        global_std = weighted_mean.std()
        l1 = (global_mean + global_std)
        l2 = (global_mean - global_std)

        # 1: above mean+std (blue)
        # 0: between mean-std and mean+std (white)
        # -1: below mean-std (red)
        col1 = comparison_1.apply(
            lambda t: -1 if (t < l2) else (1 if (t > l1) else 0))
        col1 = col1.rename("col1")
        col2 = comparison_2.apply(
            lambda t: -1 if (t < l2) else (1 if (t > l1) else 0))
        col2 = col2.rename("col2")

        df = pd.concat(
            [df, col1, col2],
            axis=1).reset_index()

        return df, global_mean, global_std

    def make_tables(self, dataset, index=True):
        for coordinate, dataset in dataset.items():
            dataset[0].to_csv(
                self.save_path / f"{coordinate}.csv", index=index)

    def make_plots(self, dataset):
        for coordinate, data_dict in dataset.items():
            dataframes = data_dict[0]
            global_mean = data_dict[1]
            global_std = data_dict[2]
            arlequin_plot(
                dataframes,
                global_mean,
                global_std,
                coordinate,
                self.save_path,
                base=self.Ac,
                unit_name=self.unit_name,
                unit_len=self.unit_len)
