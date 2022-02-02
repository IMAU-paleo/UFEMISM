#! /bin/csh -f

cd ..
./compile_all.csh

rm -rf results_202*

#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA    ISMIP-HOM/config-files/config_var_L005
#
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA_sans    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA_sans    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA_sans    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA_sans    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA_sans    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_DIVA_sans    ISMIP-HOM/config-files/config_var_L005
#
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid    ISMIP-HOM/config-files/config_var_L005
#
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid_sans    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid_sans    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid_sans    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid_sans    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid_sans    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A_hybrid_sans    ISMIP-HOM/config-files/config_var_L005

#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA    ISMIP-HOM/config-files/config_var_L005
#
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA_sans    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA_sans    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA_sans    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA_sans    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA_sans    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_DIVA_sans    ISMIP-HOM/config-files/config_var_L005
#
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid    ISMIP-HOM/config-files/config_var_L005
#
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid_sans    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid_sans    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid_sans    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid_sans    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid_sans    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B_hybrid_sans    ISMIP-HOM/config-files/config_var_L005

#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA    ISMIP-HOM/config-files/config_var_L005
#
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA_sans    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA_sans    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA_sans    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA_sans    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA_sans    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_DIVA_sans    ISMIP-HOM/config-files/config_var_L005
#
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid    ISMIP-HOM/config-files/config_var_L005
#
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid_sans    ISMIP-HOM/config-files/config_var_L160
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid_sans    ISMIP-HOM/config-files/config_var_L080
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid_sans    ISMIP-HOM/config-files/config_var_L040
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid_sans    ISMIP-HOM/config-files/config_var_L020
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid_sans    ISMIP-HOM/config-files/config_var_L010
#mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C_hybrid_sans    ISMIP-HOM/config-files/config_var_L005

mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA    ISMIP-HOM/config-files/config_var_L160
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA    ISMIP-HOM/config-files/config_var_L080
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA    ISMIP-HOM/config-files/config_var_L040
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA    ISMIP-HOM/config-files/config_var_L020
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA    ISMIP-HOM/config-files/config_var_L010
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA    ISMIP-HOM/config-files/config_var_L005

mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA_sans    ISMIP-HOM/config-files/config_var_L160
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA_sans    ISMIP-HOM/config-files/config_var_L080
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA_sans    ISMIP-HOM/config-files/config_var_L040
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA_sans    ISMIP-HOM/config-files/config_var_L020
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA_sans    ISMIP-HOM/config-files/config_var_L010
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_DIVA_sans    ISMIP-HOM/config-files/config_var_L005

mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid    ISMIP-HOM/config-files/config_var_L160
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid    ISMIP-HOM/config-files/config_var_L080
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid    ISMIP-HOM/config-files/config_var_L040
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid    ISMIP-HOM/config-files/config_var_L020
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid    ISMIP-HOM/config-files/config_var_L010
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid    ISMIP-HOM/config-files/config_var_L005

mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid_sans    ISMIP-HOM/config-files/config_var_L160
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid_sans    ISMIP-HOM/config-files/config_var_L080
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid_sans    ISMIP-HOM/config-files/config_var_L040
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid_sans    ISMIP-HOM/config-files/config_var_L020
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid_sans    ISMIP-HOM/config-files/config_var_L010
mpiexec -n 2 UFEMISM_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D_hybrid_sans    ISMIP-HOM/config-files/config_var_L005