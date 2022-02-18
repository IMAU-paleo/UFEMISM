Last updated: 2021-03-16 by Tijn Berends (c.j.berends@uu.nl)

This repository contains the source code of UFEMISM, some scripts for compiling and running the model, a collection of Matlab scripts for analysing and plotting output data, and the complete documentation.

A step-by-step instruction on how to download, compile, and run the model is provided in the documentation (documentatation/UFEMISM_documentation.pdf).


Setup on Snellius
-----------------


```bash
  module load 2021
  module load eb/4.5.2
  eblocalinstall PETSc-3.15.1-foss-2021a.eb # only necessary once
  module load foss/2021a
  module load netCDF-Fortran/4.5.3-gompi-2021a
  module load PETSc/3.15.1-foss-2021a
  module load imkl/2021.2.0-iompi-2021a
  module load OpenMPI/4.1.1-GCC-10.3.0
```

point to `Makefile_include_snellius.txt` in the Makefile
