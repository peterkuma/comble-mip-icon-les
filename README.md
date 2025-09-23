# ICON COMBLE MIP

## Overview

This repository contains files for running [COMBLE MIP experiments](https://arm-development.github.io/comble-mip/README.html) with ICON LES.

The following experiments are supported:

- FixN
- FixN_def_z0
- FixN_noice
- ProgNa (alt) [see `bin/postproc_steps/dephy` for the description of the alt case]

The experiment configuration is in files under `exp`, one file per experiment.
The common configuration (ICON namelists) for all experiments is in `bin/les`
and `input/case_definition_dephy`. The experiments can be submitted with
`bin/submit_les` *exp*, where *exp* is the experiment name, as defined in the
`exp` directory. At the end of the run, a postprocessing job is submitted
automatically. It can be also submitted manually with `bin/submit_postproc`
*exp*. The model and postprocessing can also be run in the current interactive
session with `bin/les` *exp* and `bin/postproc` *exp*. The individual
postprocessing steps are implemented in scripts under `bin/postproc_steps`. The
programs under `bin` should be always run from the main repository directory as
the current work directory. The grids are located in `input/Torus_Triangles_`\*,
and the LES forcing file is `input/COMBLE_INTERCOMPARISON_FORCING_V2.5.nc`.
These are also symlinked in the experiment's run directory as `grid.nc` and
`init_SCM.nc`, respectively.

ICON is configured to output three NetCDF files: `mean.nc`, `hourly.nc`, and
`inst.nc`. `mean.nc` is 10-min means, `hourly.nc` is hourly instantaneous, and
`inst.nc` is 10-min instantaneous. The exact output configuration is in
`bin/les`. Two ICON namelists `icon_master.namelist` and `icon.namelist`
containing the run configuration are created by `bin/les` in the experiment's
run directory. The job's standard and error output are stored in `output`, and
for the postprocessing job in `postproc_output`.

The postprocessing involve calculation of optical depth by
`bin/postproc_steps/optdepth` and conversion of the model output to the DEPHY
format by `bin/postproc_steps/dephy`, as required by COMBLE MIP. The optical
depth calculation creates a file `opt_mean.nc` and `opt_inst.nc` in the
experiment's run directory with optical depth calculated from the 10-min mean
output and 10-min instantaneous output, respectively. The final DEPHY files are
`dephy_mean.nc` and `dephy_2d.nc` for the "LES: Domain-Mean Variables" and "LES:
Two-Dimensional Fields" COMBLE MIP output, respectively. The postprocessing
steps also can be run interactively from the command line in an SSH session
started with `salloc -A bb1311 -p interactive -t 8:00:00 --mem=8G` or similar.

Additional plotting scripts are under `bin`. They can be run for all run results
in `runs` with `./run plot_runs`, and for all the MIP models under
`mip/output_les` with `./run plot_mip`. This plots curtain plots and timeseries.
The `bin/plot_curtain` and `bin/plot_timeseries` plotting script can also be
used individually.

The ICON model source code is located under `src/ICON_DWD_NWP_DEPHY`. Each
of the main experiments has its own branch in its own subdirectory. The
directory `src/ICON_DWD_NWP_DEPHY/bare` is a bare git repository which is
shared for all of them and can be pushed into with `git push origin` *branch*
from each of the experiment ICON source code subdirectories.

The Fortran postprocessing source code is located under `src/postproc`. It can
be complied with `make` (run in the same directory as the code), and requires
the gcc module to be loaded `module load gcc` beforehand.

The `submission` directory contains the final ICON model output files as submitted
to COMBLE MIP. They are the same as the corresponding `dephy.nc` files under
`runs`, except that they also contain modified metadata set with
`bin/set_dephy_meta`.

This directory is a Git repository, but not all of the files are commited to the
Git history, because the remaining files are either too large or not suitable
to be commited.
