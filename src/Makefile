# Makefile for compiling the UFEMISM model
# For "make" syntax search on "make tutorial", or see:
# http://www.gnu.org/software/make/manual/make.html

# ================================
#       Compilation settings
# ================================

# Include settings from (uncommented) file below:

#include Makefile_include_local.txt
include Makefile_include_mac.txt
#include Makefile_include_gemini.txt
#include Makefile_include_snellius.txt

# The OBJ_PATH and the MOD_PATH are defined as:
OBJ_PATH = object-files
MOD_PATH = module-files

# ==========================
#       Building rules
# ==========================

# Pattern rules for creating object files from *.f90 (fortran 90) source files:
$(OBJ_PATH)/%.o : %.f90
	$(F90) $(F90FLAGS) -c -o $@ $<

$(OBJ_PATH)/%.o : %.F90
	$(F90) $(F90FLAGS) -c -o $@ $<

# This rule is used to build an executable from the object files:
% :
	$(F90) $(F90FLAGS)    -o $@ $^ $(LDFLAGS)



# ========================================
#       Building the UFEMISM_program
# ========================================

# Adding this program to the list of all programs:
all_programs += UFEMISM_program

# List of object files
ifeq ($(DO_SELEN),yes)
UFEMISM_program_objects = \
        configuration_module.o \
        data_types_netcdf_module.o \
        data_types_module.o \
        parallel_module.o \
        petsc_module.o \
        parameters_module.o \
        utilities_module.o \
        sparse_matrix_module.o \
        mesh_memory_module.o \
        mesh_help_functions_module.o \
        mesh_Delaunay_module.o \
        mesh_operators_module.o \
        mesh_ArakawaC_module.o \
        mesh_mapping_module.o \
        mesh_creation_module.o \
        mesh_update_module.o \
        zeta_module.o \
        netcdf_module.o \
        netcdf_basic_module.o \
        netcdf_input_module.o \
        netcdf_output_module.o \
        restart_module.o \
        forcing_module.o \
        reference_fields_module.o \
        basal_conditions_and_sliding_module.o \
        SMB_module.o \
        BMB_module.o \
        ocean_module.o \
        climate_module.o \
        general_ice_model_data_module.o \
        general_sea_level_module.o \
        calving_module.o \
        ice_velocity_module.o \
        ice_thickness_module.o \
        thermodynamics_module.o \
        ice_dynamics_module.o \
        isotopes_module.o \
        bedrock_ELRA_module.o \
        SELEN_mapping_module.o \
        SELEN_harmonics_module.o \
        SELEN_fmvals_module.o \
        SELEN_fmzm90.o \
        SELEN_fm.o \
        SELEN_common_hp_module.o \
        SELEN_taboo_hp_module.o \
        SELEN_pmtf_routines.o \
        SELEN_spline_interpolation_module.o \
        SELEN_sealevel_equation_module.o \
        SELEN_main_module.o \
        tests_and_checks_module.o \
        text_output_module.o \
        UFEMISM_main_model.o \
        UFEMISM_program.o
else
UFEMISM_program_objects = \
        configuration_module.o \
        data_types_netcdf_module.o \
        data_types_module.o \
        parallel_module.o \
        petsc_module.o \
        parameters_module.o \
        utilities_module.o \
        sparse_matrix_module.o \
        mesh_memory_module.o \
        mesh_help_functions_module.o \
        mesh_Delaunay_module.o \
        mesh_operators_module.o \
        mesh_ArakawaC_module.o \
        mesh_mapping_module.o \
        mesh_creation_module.o \
        mesh_update_module.o \
        zeta_module.o \
        netcdf_module.o \
        netcdf_basic_module.o \
        netcdf_input_module.o \
        netcdf_output_module.o \
        restart_module.o \
        forcing_module.o \
        reference_fields_module.o \
        basal_conditions_and_sliding_module.o \
        SMB_module.o \
        BMB_module.o \
        ocean_module.o \
        climate_module.o \
        general_ice_model_data_module.o \
        general_sea_level_module.o \
        calving_module.o \
        ice_velocity_module.o \
        ice_thickness_module.o \
        thermodynamics_module.o \
        ice_dynamics_module.o \
        isotopes_module.o \
        bedrock_ELRA_module.o \
        tests_and_checks_module.o \
        text_output_module.o \
        UFEMISM_main_model.o \
        UFEMISM_program.o
endif

# The executable (the target) depends on the object files:
UFEMISM_program : $(UFEMISM_program_objects:%=$(OBJ_PATH)/%)



# ==================================
#       Cleaning up everything
# ==================================

.PHONY : clean

clean :
	rm -fv $(all_programs) $(OBJ_PATH)/*.o $(MOD_PATH)/*.mod *.mod *.lst



# ===============================
#       Building everything
# ===============================

.PHONY : all

all : $(all_programs)