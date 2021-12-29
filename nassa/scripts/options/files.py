import pathlib

import click


def seq_options(function):
    function = click.option(
        "-seq",
        type=str,
        multiple=True,
        help="sequence filenames. "
        "If --seqpath is specified, it is used as "
        "the path prefix for each sequence filename given.")(function)
    function = click.option(
        "-seqpath",
        type=pathlib.Path,
        help="Path to directory where "
        "sequence files can be found")(function)
    return function


def crd_options(function):
    function = click.option(
        "-crd",
        type=str,
        multiple=True,
        help="coordinate filenames. "
        "If --crdpath is specified, it is used as "
        "the path prefix for each coordinate filename given.")(function)
    function = click.option(
        "-crdpath",
        type=pathlib.Path,
        help="Path to directory where "
        "coordinate files can be found")(function)
    return function


def config_options(function):
    function = click.option(
        "--config",
        type=pathlib.Path,
        help="Path for YAML configuration file")(function)
    return function
