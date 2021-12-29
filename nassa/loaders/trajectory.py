from collections import deque

import pandas as pd


def load_serfile(ser_file, tail=True, n_lines=None):
    """
    Load single file containing a coordinate's series.

    :param str ser_file: path to .ser file.
    :param bool tail: (Default True) read the last ``n_lines`` of the file. Otherwise, read the first ``n_lines``.
    :param int n_lines: number of rows to read.
    :returns pandas.DataFrame: .ser file converted into a pandas.DataFrame table
    """
    if tail:
        with open(ser_file, "r") as f:
            # read last line of file, get index number for that line
            total_lines = int(deque(f, 1)[0].split()[0])
        extra_kwargs = dict(skiprows=total_lines - n_lines)
    else:
        extra_kwargs = dict(nrows=n_lines)
    ser_data = pd.read_csv(
        ser_file,
        header=None,
        sep='\s+',
        index_col=0,
        **extra_kwargs)
    return ser_data


def write_serfile(data, filename, indent=8, decimals=2, transpose=True):
    """Write data to same format as .ser file.
    By default, data is asumed to be in shape (n_cols, n_frames), and it's written in 8-spaced columns with values rounded to two decimals.

    :param numpy.ndarray data: output data
    :param str filename: dataset's filename
    :param indent: width of columns, defaults to 8
    :type indent: int, optional
    :param decimals: number of rounding decimals, defaults to 2
    :type decimals: int, optional
    :param transpose: transpose data array before writing. It should be used so array shape is (n_frames, n_cols). Defaults to True
    :type transpose: bool, optional
    """
    if transpose:
        data = data.T
    with open(filename, "w") as f:
        for row in data:
            s = f"{int(row[0]):>indent}"
            for elem in row[1:]:
                elem = round(elem, decimals)
                s += f"{elem:>indent}"
            s += "\n"
            f.write(s)
