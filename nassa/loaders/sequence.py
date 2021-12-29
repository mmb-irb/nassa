import pathlib

from ..entities.nucleicacid import NucleicAcid


def load_sequence(seqfile, unit_len, unit_name=None):
    """
    Load single text file containing forward and inverse-complementary sequences.

    :param str seq_file: a string with the path to the sequences file.
    :param str unit_name: name of subunit to be analyzed.
    :param int unit_len: length of subunit to be analyzed.
    :returns: NASequence object for loaded sequence.
    """
    assert isinstance(seqfile, str) or isinstance(seqfile, pathlib.Path)
    sequences = pathlib.Path(seqfile).read_text().split()
    if len(sequences) < 2:
        sequences.append(None)
        try:
            assert len(sequences) == 2
        except AssertionError:
            raise AssertionError("Error in sequence file! Check its not empty")
    nucleic_acid = NucleicAcid(
        sequence=sequences[0],
        ic_sequence=sequences[1],
        unit_name=unit_name,
        unit_len=unit_len)
    return nucleic_acid


def write_sequence(nucleic_acid, filename):
    assert isinstance(nucleic_acid, NucleicAcid)
    output = f"{nucleic_acid.sequence}\n{nucleic_acid.ic_sequence}"
    pathlib.Path(filename).write_text(output)
