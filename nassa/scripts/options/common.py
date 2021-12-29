import pathlib

import click


def common_options(function):
    function = click.option(
        '--nlines', "n_lines",
        default=2500,
        type=int,
        help="number of lines to read from files")(function)
    function = click.option(
        '--tail/--head',
        help="read starting from tail/head of file")(function)
    function = click.option(
        '--bimod/--unimod',
        help="indicate if distribution is "
        "to be considered bi/unimodal")(function)
    function = click.option(
        '-o', 'save_path',
        type=pathlib.Path,
        help="path where to save tables and plots")(function)
    function = click.option(
        '--tables', "save_tables",
        is_flag=True,
        help="save data tables as .csv files")(function)
    function = click.option(
        '--viz', "save_plots",
        is_flag=True,
        help="save data visualizations as .pdf")(function)
    function = click.option(
        '-v',
        '--verbose',
        is_flag=True,
        help="set logs to verbose")(function)
    return function
