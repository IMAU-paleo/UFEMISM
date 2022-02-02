#! /bin/csh -f

cd ..
./compile_all.csh

rm -rf results_202*

#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA    MISMIP_mod/config-files/config_var_64km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA    MISMIP_mod/config-files/config_var_40km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA    MISMIP_mod/config-files/config_var_32km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA    MISMIP_mod/config-files/config_var_20km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA    MISMIP_mod/config-files/config_var_16km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA    MISMIP_mod/config-files/config_var_10km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA    MISMIP_mod/config-files/config_var_8km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA    MISMIP_mod/config-files/config_var_5km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA    MISMIP_mod/config-files/config_var_4km

#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA_sans    MISMIP_mod/config-files/config_var_64km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA_sans    MISMIP_mod/config-files/config_var_40km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA_sans    MISMIP_mod/config-files/config_var_32km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA_sans    MISMIP_mod/config-files/config_var_20km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA_sans    MISMIP_mod/config-files/config_var_16km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA_sans    MISMIP_mod/config-files/config_var_10km
mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA_sans    MISMIP_mod/config-files/config_var_8km
mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA_sans    MISMIP_mod/config-files/config_var_5km
mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA_sans    MISMIP_mod/config-files/config_var_4km

#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid    MISMIP_mod/config-files/config_var_64km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid    MISMIP_mod/config-files/config_var_40km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid    MISMIP_mod/config-files/config_var_32km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid    MISMIP_mod/config-files/config_var_20km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid    MISMIP_mod/config-files/config_var_16km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid    MISMIP_mod/config-files/config_var_10km
mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid    MISMIP_mod/config-files/config_var_8km
mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid    MISMIP_mod/config-files/config_var_5km
mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid    MISMIP_mod/config-files/config_var_4km

#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid_sans    MISMIP_mod/config-files/config_var_64km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid_sans    MISMIP_mod/config-files/config_var_40km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid_sans    MISMIP_mod/config-files/config_var_32km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid_sans    MISMIP_mod/config-files/config_var_20km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid_sans    MISMIP_mod/config-files/config_var_16km
#mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid_sans    MISMIP_mod/config-files/config_var_10km
mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid_sans    MISMIP_mod/config-files/config_var_8km
mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid_sans    MISMIP_mod/config-files/config_var_5km
mpiexec -n 2 UFEMISM_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid_sans    MISMIP_mod/config-files/config_var_4km