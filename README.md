Nucleic Acid Sequence Structural Analysis (NASSA) 
=================================================

Structural analyses of nucleic acids. Made by the MMB Lab at IRB Barcelona.


Quickstart
------------

In order to use the library, first clone the repository

`git clone https://github.com/svletana/NASSA.git`

Then, you can install the library from the directory where you cloned the repository by doing:

`python setup.py install`

Install required libraries with `pip` by doing:

`pip install -r requirements.txt`

After installation you can use it like any other Python module, or a command line application.

Print help message with:

`nassa --help`

## Usage 

There are 7 analyses/transformations that can be executed with `nassa [OPTIONS] <analysis name> [ANALYSIS ARGS]`.

```
Commands:
  bconf      Execute BI/BII conformations analysis.
  bpcorr     Execute basepair correlations analysis.
  crdcorr    Execute coordinate correlations analysis.
  coordist   Execute Coordinate Distributions analysis.
  stiff      Execute Stiffness Distributions analysis.
  fixangles  Fix angles degree values for the given coordinate files.
  makecfg    Create a new template configuration file.
```

```
The common options for all commands are:
  -v, --verbose       set logs output to verbose
  --viz               save data visualizations as .pdf
  --tables            save data tables as .csv files
  -o PATH             path where to save tables and plots
  --bimod / --unimod  indicate if distribution analysed is to be considered bimodal/unimodal
  --tail / --head     read .ser files starting from tail/head of file
  --nlines INTEGER    number of lines to read from .ser files
  -lsu INTEGER        unit length (3 for trimers, 4 for tetramers, etc)
  -su TEXT            unit name (trimer, tetramer, etc)
  -crdpath PATH       Path to directory where coordinate files can be found
  -crd TEXT           coordinate filenames. If --crdpath is specified, it is
                      used as the path prefix for each coordinate filename
                      given.

  -seqpath PATH       Path to directory where sequence files can be found
  -seq TEXT           sequence filenames. If --seqpath is specified, it is
                      used as the path prefix for each sequence filename
                      given.

  --config PATH       Path for YAML configuration file with the options listed
  --help              Show this message and exit
```

Input files for each analysis workflow should be coordinate files (.ser files for a given helical parameter coordinate, depends on the analysis) and sequence files that correspond to each file.

The analysis could have additional arguments the go after the analysis name (options are stated before analysis name). Use the `--help` option to see them.

The options can be input manually, or with a `yaml` configuration file, like this:

`nassa --config /path/to/config.yaml <analysis_name>`

This is an example configuration file:

```
unit_name: tetramer
unit_len: 4
n_lines: 25000
tail: True
bimod: True
save_tables: True
save_plots: True
save_path: /path/to/output/directory

sequence_files:
    - /path/to/sequence/file/1
    - /path/to/sequence/file/2
    ...

coordinate_info:
    helical_parameter_1: 
        - /path/to/helical_parameter_1/file/1
        - /path/to/helical_parameter_1/file/2
        ...
    helical_parameter_2: 
        - /path/to/helical_parameter_2/file/1
        - /path/to/helical_parameter_2/file/2
        ...
```

The first file listed for each helical parameter corresponds to the first sequence file listed. Meaning, `/path/to/helical_parameter_1/file/1` and `/path/to/helical_parameter_2/file/1` both correspond to the sequence in file `/path/to/sequence/file/1`.

You can create a template for the configuration file with the command:

`nassa makecfg --config <config_filename>.yaml`

Helical parameter files are `.ser` files obtained from Curves+ and Canal softwares. Sequence files are expected to be text files with the sequence in the first line, and optionally the inverse-complement sequence in the second line, like so:

```
CGCGGACGTTCGCG
CGCGAACGTCCGCG
```

## Basepair correlation analysis
Given parameter and sequence files, calculate the correlation between basepairs, with the method corresponding to the helical parameter. That is:

```
shift: "linear"
slide: "linear"
rise: "linear"
tilt: "circular"
roll: "circular"
twist: "circular"
```

The correlation between two linear helical parameters is the Pearson correlation, between two circulars, or a linear and a circular, performs a sine/cosine transformation on the non-linear coordinate. **Circular coordinates are asumed to be in degrees**.

The result will be a heatmap and a .csv table, with the pairs of helical parameters in the x-axis and the neighboring basepairs indicated in the y-axis. The colorbar on the right indicates the correlation magnitude. The output plot won't indicate the neighboring basepairs, for that it's best to look at the generated .csv file, or use [biobb_dna](https://github.com/bioexcel/biobb_dna) and perform a single-basepair analysis.

## Coordinate correlation analysis

Given all the helical parameter files obtained from one trajectory, calculate the correlation between h.p. coordinates for each subunit. What's different from the basepair correlation analysis, is that it only compares helical parameters to calculate correlation, but the subunit is the same for each value.

The input for this correlation analysis is the same as in the previous correlation analysis. The output is a .csv table with the correlation magnitudes, and a heatmap with subunits on both axes, each axis corresponding to different h.p. coordinates.

## Coordinate distribution analysis

Performs distribution calculations of each helical parameter, outputs a table with all the distribution infrmation and creates heatmaps, or *arlequin plots*. If `--bimod` option is set, the heatmap cells of the plot will indicate if bimodality is present in each subunit, for that helical parameter, by dividing it between its upper and lower states. Otherwise, it will treat all subunits as having a single state. Colors in the heatmap indicate if the subunit's cordinate value is within the average +/- 1 stdev, above avg + 1 stdev, or below avg - 1 stdev.


## Stiffness distribution analysis

Calculates distribution of stiffness value for each coordinate. The input and results is the same as with the coordinate distribution analysis.

## BI/BII conformations analysis

Calculates distribution of BI/BII conformation analysis for each subunit. The results are a table with the percentage of BI/BII states present in each subunit, as well as a heatmap.

