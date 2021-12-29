import click


def subunit_options(function):
    function = click.option(
        '-su', "unit_name",
        default="tetramer",
        type=str,
        help="unit name (trimer, tetramer, etc)")(function)
    function = click.option(
        '-lsu', "unit_len",
        default=4,
        type=int,
        help="unit length (3 for trimer, 4 for tetramer, etc)")(function)
    return function
