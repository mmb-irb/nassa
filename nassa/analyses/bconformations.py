import pandas as pd

from .base import Base
from .utils.heatmaps import bconf_heatmap
from .utils.angle_utils import fix_angle_range
from ..loaders.sequence import load_sequence
from ..loaders.trajectory import load_serfile


class BConformations(Base):
    """Execution plan and methods for BI/BII conformations analysis pipeline"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        valid_coordinates = ["epsilC", "zetaC", "epsilW", "zetaW"]
        for coord in self.coordinate_info.keys():
            try:
                assert coord in valid_coordinates
            except AssertionError as e:
                raise ValueError(
                    f"{coord} is not a valid coordinate! "
                    "Please rename coordinates in your configuration file "
                    f"to match any of {valid_coordinates}") from e

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
        angles_df = []
        # get dataframe for each coordinate
        for traj, seq in enumerate(sequences):
            # start reading from the 4th column:
            # skip index, first two bases and first flanks
            start = 3 + seq.flanksize
            # skip last two bases
            end = seq.size - (seq.baselen + seq.flanksize)
            # select relevant subset of columns
            epsilC = data["epsilc"][traj][start:end]
            zetaC = data["zetac"][traj][start:end]
            epsilW = data["epsilw"][traj][start:end]
            zetaW = data["zetaw"][traj][start:end]
            traj_df = self.get_angles_difference(
                seq,
                epsilC,
                zetaC,
                epsilW,
                zetaW)
            angles_df.append(traj_df)
        angles_df = pd.concat(angles_df, axis=1)
        # percentages BI
        B_I = (angles_df < 0).sum(axis=0) * 100 / self.n_lines
        # clean dataset
        B_I = B_I[~B_I.index.duplicated(keep='first')]
        B_I = B_I.reset_index()
        B_I = B_I.rename(columns={"index": self.unit_name, 0: "pct"})
        return B_I

    def get_angles_difference(self, seq, epsilC, zetaC, epsilW, zetaW):
        # get list of tetramers, except first and last two bases
        all_subunits = seq.all_subunits[2:-2]
        all_ic_subunits = seq.all_ic_subunits[2:-2]

        # concatenate zeta and epsil arrays
        zeta = pd.concat([
            pd.DataFrame(zetaW.T),
            pd.DataFrame(zetaC[::-1].T)],
            axis=1)
        zeta.columns = all_subunits * 2
        epsil = pd.concat([
            pd.DataFrame(epsilW.T),
            pd.DataFrame(epsilC[::-1].T)],
            axis=1)
        epsil.columns = all_subunits * 2
        # difference between epsilon and zeta coordinates
        diff = epsil - zeta
        diff = diff.applymap(lambda x: fix_angle_range(x, domain=[-180, 180]))

        # repeat with inverse-complementary sequence
        zeta_ic = pd.concat([
            pd.DataFrame(zetaW.T),
            pd.DataFrame(zetaC[::-1].T)],
            axis=1)
        zeta_ic.columns = all_ic_subunits * 2
        epsil_ic = pd.concat([
            pd.DataFrame(epsilW.T),
            pd.DataFrame(epsilC[::-1].T)],
            axis=1)
        epsil_ic.columns = all_ic_subunits * 2
        diff_ic = epsil_ic - zeta_ic
        diff_ic = diff_ic.applymap(
            lambda x: fix_angle_range(x, domain=[-180, 180]))

        diff_df = pd.concat([diff, diff_ic], axis=1)

        return diff_df

    def make_tables(self, B_I):
        B_I.to_csv(self.save_path / "BI.csv", index=False)
        # create BII
        B_II = B_I.copy()
        B_II["pct"] = 100 - B_I["pct"]
        B_II.to_csv(self.save_path / "BII.csv", index=False)

    def make_plots(self, B_I):
        bconf_heatmap(B_I, "BI", self.save_path, self.unit_len, self.Ac)
        B_II = B_I.copy()
        B_II["pct"] = 100 - B_I["pct"]
        bconf_heatmap(B_II, "BII", self.save_path, self.unit_len, self.Ac)
