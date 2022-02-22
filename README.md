Last updated: 2022-02-18 by Jorjo Bernales (j.a.bernalesconcha@uu.nl)

This repository contains the source code of UFEMISM, some scripts for
compiling and running the model, a collection of Matlab scripts for
analysing and plotting output data, and the complete documentation.

A step-by-step instruction on how to download, compile, and run the model
is provided in the documentation (documentation/UFEMISM_documentation.pdf).

Quick Start
-----------

This assumes a Linux or MacOS system with a working installation of
NetCDF, PETSc, an MPI implementation (e.g. OpenMPI), and LAPACK.

Example for bash:

```bash
  # 1. Copy compilation and run files from templates/ to UFEMISM/
  cp templates/compile_all.sh .
  cp templates/run_UFEMISM.sh .

  # 2. Modify src/Makefile_include_local.txt to your local settings and
  #    compilation preferences    

  # 3. Modify src/Makefile so it points to Makefile_include_local.txt

  # 4. Compile the model
  ./compile_all.sh

  # 5. Run a test simulation
  ./run_UFEMISM.sh
```

Setup on Snellius
-----------------

To compile and run the model on the Snellius supercomputer, the following
modules/tools need to be superloaded: 

```bash
  module load 2021
  module load eb/4.5.2
  eblocalinstall PETSc-3.15.1-foss-2021a.eb # ==> only necessary once <==
  module load foss/2021a
  module load netCDF-Fortran/4.5.3-gompi-2021a
  module load PETSc/3.15.1-foss-2021a
  module load imkl/2021.2.0-iompi-2021a
  module load OpenMPI/4.1.1-GCC-10.3.0
```

These modules are loaded within the compilation and running scripts as a safety
measure, but their loading can be outsourced to a separate file to avoid
repetitive reloading of the modules under frequent recompilations.

Then follow the steps of the Quick Start above, using the `*snellius*` files
from `templates/` instead, and modify `src/Makefile` so it points to
`Makefile_include_snellius.txt`. Adapt `run_UFEMISM_snellius.sh` to your needs,
but instead send your test simulation to the batch system by executing:

```bash
  # 5. Run a test simulation
  ./submit_UFEMISM_snellius.sh
```

Profiling on Snellius
---------------------

```bash
  module load VTune/2021.6.0
  source ${EBROOTVTUNE}/setvars.sh
  vtune -collect hotspots -r vtune_output ./UFEMISM_program config_test
```
