from itertools import product

import pandas as pd
import numpy as np
from ..entities.nucleicacid import NucleicAcid
from .base import Base
from ..loaders.sequence import load_sequence
from ..loaders.trajectory import load_serfile
from .utils.heatmaps import basepair_plot
import os
import glob
import shutil
from Bio.Seq import Seq
import numpy.ma as ma


class BasePairCorrelation(Base):
    """
    Execution plan and methods for basepair correlation analysis.
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

        #Iterate over the different parameters: shift, twist, roll, slide, tilt, rise
        for coord in self.coordinate_info.keys():
            crd_data = []
            for ser_file in self.coordinate_info[coord]:
                ser = load_serfile(
                    ser_file,
                    self.tail,
                    self.n_lines)
                crd_data.append(ser)
            extracted[coord.lower()] = crd_data
            self.logger.info(f"loaded {len(crd_data)} files for coordinate <{coord}>")
            
        return extracted

    def extraer_tetramero_central(hexamero):
                return hexamero[1:5]
    
    def transform(self, data):
        sequences = data.pop("sequences")
        # iterate over trajectories
        corr_results = {}
        for traj, seq in enumerate(sequences):
            trajectory_series = {coord.lower(): data[coord][traj]
                                 for coord in data.keys()}
            coordinate_corr = self.iterate_trajectory(
                seq, trajectory_series)
            corr_results[seq.sequence] = coordinate_corr
        joined_df = []
        for seq, val in corr_results.items():
            df = pd.DataFrame.from_dict(val).T
            joined_df.append(df)
        joined_df = pd.concat(joined_df)
        return joined_df

    def iterate_trajectory(self, sequence, coordinates):
        coordinate_correlations = {}
        # iterate over subunits
        start = 2 + sequence.flanksize
        end = sequence.size - (2 + sequence.baselen + sequence.flanksize - 1)
        # subtract 1 from end in order to stop at last sub1unit
        for idx in range(start, end-1):
            subunit = sequence.get_subunit(idx)
            next_subunit = sequence.get_subunit(idx+1)
            if subunit in coordinate_correlations:
                self.logger.info(
                    f"skipping repeated {self.unit_name} {subunit}...")
                continue
            self.logger.info(
                f"analyzing {self.unit_name} {subunit}...")
            # add 1 to idx since .ser table includes an index
            unit_df = pd.DataFrame(
                {coord: coordinates[coord][idx+1]
                 for coord in coordinates.keys()})
            next_unit_df = pd.DataFrame(
                {coord: coordinates[coord][idx+2]
                 for coord in coordinates.keys()})
            crd_corr = self.get_correlation(next_unit_df, unit_df) # inverse order ? 
            coordinate_correlations[f"{subunit}/{next_subunit}"] = crd_corr
        return coordinate_correlations

    def get_correlation(self, unit, next_unit):
        method = {
            "shift": "linear",
            "slide": "linear",
            "rise": "linear",
            "tilt": "circular",
            "roll": "circular",
            "twist": "circular"}
        coordinates = method.keys()
        combos = product(coordinates, repeat=2)
        result = {}
        for crd1, crd2 in combos:
            method1 = method[crd1]
            method2 = method[crd2]
            arr1 = unit[crd1]
            arr2 = next_unit[crd2]
            value = self.get_corr_by_method(
                method1,
                method2,
                arr1,
                arr2)
            result[f"{crd1}/{crd2}"] = value
        return result

    def get_corr_by_method(self, method1, method2, arr1, arr2):
        if method1 == "circular" and method2 == "linear":
            value = self.circlineal(arr2, arr1)
        if method1 == "linear" and method2 == "circular":
            value = self.circlineal(arr1, arr2)
        elif method1 == "circular" and method2 == "circular":
            value = self.circular(arr1, arr2)
        else:
            value = ma.corrcoef(ma.masked_invalid(arr1), ma.masked_invalid(arr2))[1, 0]
        return value

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
        x1 = x1.to_numpy()
        x2 = x2.to_numpy()
        rc = ma.corrcoef(ma.masked_invalid(x1), ma.masked_invalid(np.cos(x2)))[1, 0]
        rs = ma.corrcoef(ma.masked_invalid(x1), ma.masked_invalid(np.sin(x2)))[1, 0]
        rcs = ma.corrcoef(ma.masked_invalid(np.sin(x2)), ma.masked_invalid(np.cos(x2)))[1, 0]
        num = (rc ** 2) + (rs ** 2) - 2 * rc * rs * rcs
        den = 1 - (rcs ** 2)
        correlation = np.sqrt(num / den)
        if ma.corrcoef(ma.masked_invalid(x1), ma.masked_invalid(x2))[1, 0] < 0:
            correlation *= -1
        return correlation

    def make_tables(self, dataset):
        dataset.to_csv(self.save_path / "all_basepairs.csv")

    def make_plots(self, dataset):
        basepair_plot(dataset, "all_basepairs", self.save_path)
    
