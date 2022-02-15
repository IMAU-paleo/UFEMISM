Last updated: 2021-03-16 by Tijn Berends (c.j.berends@uu.nl)

This repository contains the source code of UFEMISM, some scripts for compiling and running the model, a collection of Matlab scripts for analysing and plotting output data, and the complete documentation.

A step-by-step instruction on how to download, compile, and run the model is provided in the documentation (documentatation/UFEMISM_documentation.pdf).


Setup on Snellius
-----------------


```bash
  module load 2021
  module load eb/4.5.2
  eblocalinstall PETSc-3.15.1-foss-2021a.eb
  module load foss/2021a
  module load netCDF-Fortran/4.5.3-gompi-2021a
  module load PETSc/3.15.1-foss-2021a
```

point to `Makefile_include_snellius.txt` in the Makefile

Profiling on Snellius
---------------------

```bash
  module load VTune/2021.6.0
  source ${EBROOTVTUNE}/setvars.sh
  vtune -collect hotspots -r vtune_output ./UFEMISM_program config_test
```
