#! /bin/csh -f

./compile_all.csh

rm -rf results_2021*

mpiexec -n 2 UFEMISM_program config-files/config_test
