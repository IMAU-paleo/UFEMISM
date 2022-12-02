#! /bin/csh -f

./compile_all_mac.csh

rm -rf results_202*

## ISMIP-HOM A

mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_160
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_080
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_040
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_020
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_010
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_005

mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_160
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_080
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_040
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_020
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_010
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_005

mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_160
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_080
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_040
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_020
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_010
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_A_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_005

## ISMIP-HOM C

mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_160
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_080
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_040
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_020
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_010
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_SIASSA_template     benchmarks/ISMIP-HOM/config-files/config_var_005

mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_160
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_080
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_040
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_020
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_010
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_DIVA_template     benchmarks/ISMIP-HOM/config-files/config_var_005

mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_160
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_080
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_040
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_020
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_010
mpiexec  -n 2   UFEMISM_program   benchmarks/ISMIP-HOM/config-files/config_ISMIP_HOM_C_BPA_template     benchmarks/ISMIP-HOM/config-files/config_var_005