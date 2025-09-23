# ICON COMBLE MIP

Author: [Peter Kuma](peter@peterkuma.net) based on code inherited from [Anna
Possner](mailto:apossner@iau.uni-frankfurt.de).

## Overview

This repository contains files for running the [COMBLE MIP
experiments](https://arm-development.github.io/comble-mip/README.html) with ICON
LES.

The following experiments are supported:

- **FixN**
- **FixN_def_z0**
- **FixN_noice**
- **ProgNa** (alt) [see `ALT_PROGNA` in `bin/postproc_steps/dephy` for the
  description of the alt case]

Some information about the ICON runs for the COMBLE MIP is in the poster
[Exploring sensitivity to ice nucleating particles and secondary ice production
during COMBLE in idealised ICON large eddy
simulations](https://zenodo.org/records/15174379).

## Anna's files

Anna's original files for the COMBLE MIP experiment are under:

- `/home/b/b380903/code/ICON_DWD_NWP_DEPHY_FixN/icon-nwp`: ICON code.
- `/work/bb1358/possnera/COMBLE`: Run scripts, namelists, input files, and
simulation output.
`/home/b/b380903/postproc/COMBLE`: Postprocessing.

This repository inherits from the above, but the process is more consolidated,
various fixes have been applied to the ICON source code (Section "ICON source
code"), and 2D DEPHY output has been added. This repository largely inherits
from and supersedes Anna's files, but there may be some useful parts that are
under Anna's files but not here.

## Code repository

This directory is a Git repository, but not all of the files are committed to the
Git history, because the remaining files are either too large or not suitable to
be committed. Notably, scripts under `bin`, postprocessing code under `src`, and
experiment configuration under `exp` are committed, while the model output under
`runs` and `submission` is not, nor are ICON repositories under
`src/ICON_DWD_NWP_DEPHY` and the COMBLE MIP repository under `arm-comble-mip`.

## Required software

To run the scripts under `bin`, some Python packages are required:

- [ds-format](https://ds-format.peterkuma.net)
- [aquarius-time](https://github.com/peterkuma/aquarius-time/)
- [rsool](https://github.com/peterkuma/rstool)

These can be installed from PyPI with `pip3 install ds-format aquarius-time
rstool`.

## Experiment configuration

The experiment configuration is in files under `exp`, one file per experiment.
The common configuration (ICON namelists) for all experiments is in `bin/les`
and `input/case_definition_dephy`. Only a few experiments were used for the
submission to COMBLE MIP:

- `FixN_as0.68`: The **FixN** experiment. This sets the `alpha_spacefilling`
  parameter to 0.68 to reduce graupel formation.
- `FixN_noice`: The **FixN noice** experiment. The `alpha_spacefilling`
  parameter does not need to be set (no graupel).
- `ProgNa_as0.68`: The **ProgNa** experiment.

The experiment configuration files are shell scripts included by the `bin/les`
script. The experiment configuration supports basic inheritance through the
inclusion of another configuration, such as `. exp/FixN`. The experiments marked
as `small` are for the "small" grid `input/Torus_Triangles_10x10_100m.nc` (1 km
domain), while the other ones are for the "large" grid
`input/Torus_Triangles_182x182_141m.nc`. Note that the large grid is considered
a small "toy" domain in the context of the COMBLE MIP manuscript.

## ICON model input and output

The model experiment run directories are stored in `runs` and under
`runs_work`. The former is on the scratch partition and subject to old file
removal policies on Levante but is faster for reading and writing data.
`runs_work` is a copy of `runs` on the work partition (project bb1311) made with
rsync. This is permanent (as long as the project exists) but slower. After
changes in `runs`, rsync should be run in due time (within several weeks) to
ensure that the output is preserved: `rsync -av runs/ runs_work/`.

ICON is configured to output three NetCDF files: `mean.nc`, `hourly.nc`, and
`inst.nc`. `mean.nc` is 10-min means, `hourly.nc` is hourly instantaneous, and
`inst.nc` is 10-min instantaneous. The exact output configuration is in
`bin/les`. Two ICON namelists, `icon_master.namelist` and `icon.namelist`,
containing the run configuration, are created by `bin/les` in the experiment's
run directory. Other ICON input files are symlinked in the run directory.
Notably, the grids are located in `input/Torus_Triangles_`\*, and the LES
forcing file is `input/COMBLE_INTERCOMPARISON_FORCING_V2.5.nc`. These are
symlinked in the experiment's run directory as `grid.nc` and `init_SCM.nc`,
respectively. The job's standard and error output are stored in `output`, and
for the postprocessing job, in `postproc_output`.

## Postprocessing

The postprocessing involves calculation of optical depth by
`bin/postproc_steps/optdepth` and conversion of the model output to the DEPHY
format by `bin/postproc_steps/dephy`, as required by the COMBLE MIP. The optical
depth calculation creates files `opt_mean.nc` and `opt_inst.nc` in the
experiment's run directory with optical depth calculated from the 10-min mean
output and 10-min instantaneous output, respectively. The final DEPHY files are
`dephy_mean.nc` and `dephy_2d.nc` for the "LES: Domain-Mean Variables" and "LES:
Two-Dimensional Fields" COMBLE MIP output, respectively. The postprocessing
steps also can be run interactively from the command line in an SSH session
started with `salloc -A bb1311 -p interactive -t 8:00:00 --mem=8G` or similar.

The Fortran postprocessing source code is located under `src/postproc`. It can
be compiled with `make` (run in the same directory as the code) and requires
the gcc module to be loaded with `module load gcc` beforehand.

The postprocessing steps `add_optdepth`, `add_optdepth_vsedi`, and `vsedi`
are currently not applied by `bin/postproc`. The work done by `add_optdepth`
is now done by `dephy`. I don't know if vsedi will be needed in the future.

## Submitting experiments

The experiments can be submitted with `bin/submit_les` *exp*, where *exp* is the
experiment name, as defined in the `exp` directory. At the end of the run, a
postprocessing job is submitted automatically. It can also be submitted manually
with `bin/submit_postproc` *exp*. The model and postprocessing can also be run
in the current interactive session with `bin/les` *exp* and `bin/postproc`
*exp*. The individual postprocessing steps are implemented in scripts under
`bin/postproc_steps`, and these can also be run individually in an interactive
session. The programs under `bin` should always be run from the main repository
directory as the current work directory.

## Plotting

Some plotting scripts are under `bin`. They can be run for all run results in
`runs` with `./run plot_runs` and for all the MIP models under `mip/output_les`
with `./run plot_mip`. The following plotting scripts are available:

- `plot_2d`: Plot 2D surface fields from DEPHY 2D files or directly from the
  ICON model output files (`mean.nc`, `hourly.nc`, or `inst.nc`). The fields
  plotted depend on the model output available in the file.

- `plot_curtain`: Plot "curtain" (time vs. height) plots from the DEPHY mean
  files are `mean.nc`.

- `plot_timeseries`: Plot time series plots from the DEPHY mean files or
  `mean.nc`.

The plotting script and the `run` script should always be executed from the main
repository directory as the current working directory.

## Input files

Some ancillary input files are stored under `input`:

- `case_definition_dephy`: ICON configuration sourced by `bin/les`.
- `COMBLE CAO LES-SCM Intercomparison Output Variables.xlsx`: Original variables
  spreadsheet from COMBLE MIP.
- `COMBLE_INTERCOMPARISON_FORCING_V2.5.nc`: ICON forcing, as defined by COMBLE
  MIP.
- `comble_vars_2d.csv`: Definition of 2D DEPHY variables needed by
  `bin/postproc_steps/dephy`.
- `comble_vars_mean.csv`: Definition of mean DEPHY variables needed by
  `bin/postproc_steps/dephy`.
- `Torus_Triangles_*.nc`: ICON grid definitions.
- `track.nc`: Track of the COMBLE MIP cold air outbreak simulation.

## ICON source code

The ICON model source code is located under `src/ICON_DWD_NWP_DEPHY`. Each
of the main experiments has its own branch in its own subdirectory. The
directory `src/ICON_DWD_NWP_DEPHY/bare` is a bare Git repository that is
shared for all of them and can be pushed into with `git push origin` *branch*
from each of the experiment ICON source code subdirectories.

- `bare`: Git "bare" repository in which all of the branches for the experiments
  are stored by `git push` from the repositories listed below.
- `baseline`: Baseline for all of the experiments. This inherits from the
  `icon-nwp-scm` branch of the [main ICON
  repository](https://gitlab.dkrz.de/icon/icon-nwp) and adds a fix for surface
  flux calculation (commit b8e5168), use of `parturs` instead of `set_scm_bnd`
  for surface flux calculations (commit 576ab88), and changes required for the
  COMBLE forcing and vertical grid (commits dd6819 and c76eac).
- `baseline_fix_z0`: This inherits from `baseline` but also adds fixed surface
  roughness `z0h` (commit 446139).
- `FixN`: Repository for the FixN experiment. This is based on `baseline_fix_z0`
  and adds fixed hydrometeor number concentration (commit 86c788) and turns off
  ice multiplication (commit 108088). The latter makes little difference, but
  secondary ice production (SIP) is supposed to be off for FixN.
- `FixN_def_z0`: The same as above, but for the default z0 experiment, which
  means it does not include the commit for fixed surface roughness (commit
  446139).
- `icon-nwp-scm`: The original branch of the same name from [ICON
  GitLab](icon-nwp-scm).
- `ProgNa`: The ProgNa experiment. This is based on `baseline_fix_z0` and also
  sets initial aerosol parameters required for this experiment (commit a1d21a).
  TODO: This needs some more attention before it complies with the COMBLE MIP
  requirements for the experiment.
- `ProgNaNi`: The same as above, but for the ProgNaNi experiment. TODO: This
  needs some more attention before it complies with the COMBLE MIP requirements
  for the experiment.

## Submission

The `submission` directory contains the final ICON model output files as
submitted to COMBLE MIP. They are the same as the corresponding `dephy.nc` files
under `runs`, except that they also contain modified metadata set with
`bin/set_dephy_meta`. The current version of `bin/postproc_steps/dephy` now sets
the metadata correctly, and `bin/set_dephy_meta` does not need to be applied. The
files under `submission` have been committed to the [COMBLE MIP GitHub
repository](https://github.com/ARM-Development/comble-mip/).

## TODO

As of 23 Sep 2025, there is an ongoing discussion with [Gaurav
Dogra](mailto:gaurav.dogra@ipsl.fr) and [Tomi
Raatikainen](Tomi.Raatikainen@fmi.fi) about collaboration regarding SIP
simulations for COMBLE MIP. They should be using their respective models for the
simulation.

Anna has wanted to add the McCluskey scheme to the **ProgNa** and **ProgNaNi**
experiments. The code is supposedly somewhere in the same directory as FixN in
Anna's COMBLE MIP files.

The **ProgNa** and **ProgNaNi** may also need some additional work and checks
before submitting to COMBLE MIP regarding how diagnostic and prognostic ice
should be specified.

The SIP experiments should be conducted with the SIP implementation as in the
[ICON MICRO](https://github.com/apossner/ICON_MICRO) repository produced by Anna
and Kevin Pfannkuch. [Another version](https://github.com/peterkuma/ICON_MICRO)
of the repository exists, which also includes the whole Git history of ICON as on
the DKRZ GitLab. Both of these repositories are private and require an invite.
Contact [Peter Kuma](mailto:peter@peterkuma.net) or Anna.

## Other stuff

`src/grid-generator` is a grid generator that can produce torus grids. This is
not needed unless grids other than `input/Torus_Triangles_*.nc` are needed.

The `diff` directory contains some diff files for the ICON code. This is already
incorporated in `src/ICON_DWD_NWP_DEPHY`.

The `arm-comble-mip` directory is the [COMBLE MIP
repository](https://github.com/arm-development/comble-mip/) as on GitHub.

Some old code is located in `bin/_`, `bin/postproc_step/_`, and `src/_`. Use at
your own peril.

`runs_AP` is a symlink to the original runs performed by Anna, and `runs_AP_rw`
are some additional files/plots produced for the former, mostly for testing
purposes only.

The `plots` directory contains some old plots.

## Contact

For more information, you can contact either [Anna
Possner](mailto:apossner@iau.uni-frankfurt.de) or [Peter
Kuma](mailto:peter@peterkuma.net).
