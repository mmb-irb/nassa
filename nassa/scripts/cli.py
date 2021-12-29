import pathlib

import yaml
import click
from numpy import vectorize

from ..analyses.coordinatedistributions import CoordinateDistributions
from ..analyses.stiffnessdistributions import StiffnessDistributions
from ..analyses.bconformations import BConformations
from ..analyses.coordinatecorrelation import CoordinateCorrelation
from ..analyses.basepairscorrelation import BasePairCorrelation
from ..analyses.utils.angle_utils import fix_angle_range
from ..loaders.trajectory import load_serfile, write_serfile
from .options import (
    common_options,
    subunit_options,
    seq_options,
    crd_options,
    config_options)


def parse_pathfiles(names, prefix):
    if prefix is not None:
        names = [prefix / c for c in names]
    else:
        names = [pathlib.Path(c) for c in names]
    return names


@click.command()
@click.pass_obj
def coordist(kwargs):
    """Execute Coordinate Distributions analysis."""
    click.echo("Executing coordinate distributions analysis")
    cdist = CoordinateDistributions(**kwargs)
    cdist.run()


@click.command()
@click.pass_obj
def stiff(kwargs):
    """Execute Stiffness Distributions analysis."""
    click.echo("Executing stiffness distributions analysis")
    stiffnessdist = StiffnessDistributions(**kwargs)
    stiffnessdist.run()


@click.command()
@click.pass_obj
def bconf(kwargs):
    """Execute BI/BII conformations analysis."""
    click.echo("Executing BI/BII conformations analysis")
    bconcoformations = BConformations(**kwargs)
    bconcoformations.run()


@click.command()
@click.pass_obj
def bpcorr(kwargs):
    """Execute basepair correlations analysis."""
    click.echo("Executing basepair correlations analysis")
    correlations = BasePairCorrelation(**kwargs)
    correlations.run()


@click.command()
@click.pass_obj
def crdcorr(kwargs):
    """Execute coordinate correlations analysis."""
    click.echo("Executing coordinate correlations analysis")
    correlations = CoordinateCorrelation(**kwargs)
    correlations.run()


@click.command()
@config_options
def makecfg(config):
    """Create a new template configuration file.
    If filename is not given, config.yaml is created in the current directory.

    :param str config: relative path to configuration filename"""
    if not config:
        config = pathlib.Path.cwd() / "config.yaml"
    template = [
        "unit_name: ",
        "unit_len: ",
        "dna: ",
        "n_lines: ",
        "tail: ",
        "bimod: ",
        "save_tables: ",
        "save_plots: ",
        "save_path: ",
        "sequence_files: ",
        "    - # path to sequence file",
        "coordinate_info: ",
        "    coord1: ",
        "        - # path to coord1 file"
    ]
    click.echo(
        f"writing configuration template file to {config}")
    with config.open("w") as configfile:
        for line in template:
            configfile.write(line)
            configfile.write("\n")


@click.group()
@common_options
@subunit_options
@crd_options
@seq_options
@config_options
@click.pass_context
def entry_point(ctx, **kwargs):
    # pop arguments not used in workflows
    configfile = kwargs.pop("config")
    seq = kwargs.pop("seq")
    crd = kwargs.pop("crd")
    seqpath = kwargs.pop("seqpath")
    crdpath = kwargs.pop("crdpath")
    if configfile:
        # parse config
        click.echo("Reading input from configuration file")
        try:
            with configfile.open("r") as ymlfile:
                base_config = yaml.load(
                    ymlfile, Loader=yaml.FullLoader)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"Configuration file {configfile} not found!")
        kwargs.update(base_config)
    else:
        click.echo("Reading input from arguments")
        # create coordinate info and sequence files
        seq = parse_pathfiles(seq, seqpath)
        crd = parse_pathfiles(crd, crdpath)
        crd = {"coord": crd}
        kwargs.update(dict(
            coordinate_info=crd,
            sequence_files=seq
        ))
    ctx.obj = kwargs


entry_point.add_command(coordist)
entry_point.add_command(bconf)
entry_point.add_command(bpcorr)
entry_point.add_command(crdcorr)
entry_point.add_command(stiff)
entry_point.add_command(makecfg)
