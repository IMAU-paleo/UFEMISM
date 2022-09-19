#    Makefile for compiling the IMAU_ICE model

#    For "make" syntax search on "make tutorial", or see:
#     http://www.gnu.org/software/make/manual/make.html

# For compiling setting, Makefile will point to the uncommented file below:
include Makefile_include_local.txt
#include Makefile_include_mac.txt
#include Makefile_include_snellius.txt

# The OBJ_PATH and the MOD_PATH are defined:
OBJ_PATH = object-files
MOD_PATH = module-files

# Pattern rules for creating object files from *.f90 (fortran 90) source files
$(OBJ_PATH)/%.o : %.f90 Makefile
	$(F90) $(F90FLAGS) -J${MOD_PATH} -fbackslash -c -o $@ $<


# ================================================================
#       Building the UFEMISM_program
# ================================================================

# Adding this program to the list of all programs:
all_programs += UFEMISM_program

# List of object files
UFEMISM_program_objects = \
        configuration_module.o \
        reallocate.o \
        data_types_netcdf_module.o \
        data_types_module.o \
        parallel_module.o \
        mpi_module.o \
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
        mesh_single_module.o \
        zeta_module.o \
        netcdf_module.o \
        forcing_module.o \
        climate_module.o \
        ocean_module.o \
        reference_fields_module.o \
        basal_conditions_and_sliding_module.o \
        general_ice_model_data_module.o \
        general_sea_level_module.o \
        thermodynamics_module.o \
        ice_velocity_module.o \
        ice_thickness_module.o \
        ice_dynamics_module.o \
        BMB_module.o \
        SMB_module.o \
        scalar_data_output_module.o \
        UFEMISM_main_model.o \
        UFEMISM_program.o

# This rule is used to build an executable from the object files:
# The executable (the target) depends on the object files:
UFEMISM_program : $(UFEMISM_program_objects:%=$(OBJ_PATH)/%)
	$(F90) $(F90FLAGS)    -o $@ $^ $(LDFLAGS)


# =============================================================
#       Cleaning up everything
# =============================================================

.PHONY : clean

clean :
	rm -fv $(all_programs) $(OBJ_PATH)/*.o $(MOD_PATH)/*.mod *.mod *.lst .depend



# =============================================================
#       Building everything
# =============================================================

.PHONY : all

all : $(all_programs)



# =============================================================
#       UFEMISM library for unit testing:
# =============================================================

TEST_DIR = ../tests
TEST_OBJ_DIR = ${TEST_DIR}/obj
TEST_MOD_DIR = ${TEST_DIR}/mod
TEST_OBJS = $(filter-out ${TEST_OBJ_DIR}/UFEMISM_program.o,$(UFEMISM_program_objects:%=$(TEST_OBJ_DIR)/%))

$(TEST_OBJ_DIR)/%.o : %.f90 Makefile
	$(F90) $(F90FLAGS) -J${TEST_MOD_DIR} -fbackslash -fprofile-arcs -ftest-coverage -c -o $@ $<

unit_tests_lib: ../tests/UFEMISM.a

../tests/UFEMISM.a : $(TEST_OBJS)
	ar rcs $@ $^


# =============================================================
#      Make dependency graph, better re-compilation + multi
#      threaded compilation (-j flag)
#      Does nothing when makedepf90 is not present
# =============================================================

ifeq ($(shell command -v makedepf90 2>&1 > /dev/null; echo $$?), 0)
include .depend
.depend: Makefile *.f90
	  makedepf90 -b $(OBJ_PATH) *.f90 > .depend;
	  makedepf90 -b $(TEST_OBJ_DIR) *.f90 >> .depend;
else
.NOTPARALLEL:
.depend: 
	$(warning makedepf90 disabled, enabling .NOTPARALLEL: for the whole makefile )
	touch .depend;
endif
