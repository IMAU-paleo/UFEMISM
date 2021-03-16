#    Makefile for compiling the UFEMISM odel

#    For "make" syntax search on "make tutorial", or see:
#     http://www.gnu.org/software/make/manual/make.html

include Makefile.mpif90
#include Makefile.gfortran
#include Makefile.gfortran-mp-4.8-mac
#include Makefile.gfortran-mp-4.9-mac
#include Makefile.gfortran-lisa
#include Makefile.ifort-lisa
#include Makefile.ifort-staff.science

# The OBJ_PATH and the MOD_PATH are defined:
OBJ_PATH = object-files
MOD_PATH = module-files

# Pattern rules for creating object files from *.f90 (fortran 90) source files
$(OBJ_PATH)/%.o : %.f90
	$(F90) $(F90FLAGS) -c -o $@ $<

# This rule is used to build an executable from the object files:
% :
	$(F90) $(F90FLAGS)    -o $@ $^ $(LDFLAGS)



# ================================================================
#       Building the UFEMISM_program
# ================================================================

# Adding this program to the list of all programs:
all_programs += UFEMISM_program

# List of object files for fortran_example_model_using_namelist_program:
UFEMISM_program_objects = \
        configuration_module.o \
        data_types_module.o \
        parallel_module.o \
        parameters_module.o \
        zeta_module.o \
        netcdf_module.o \
        mesh_memory_module.o \
        mesh_help_functions_module.o \
        mesh_Delaunay_module.o \
        mesh_ArakawaC_module.o \
        mesh_derivatives_module.o \
        mesh_mapping_module.o \
        mesh_rectangular_module.o \
        mesh_creation_module.o \
        mesh_update_module.o \
        mesh_smooth_module.o \
        forcing_module.o \
        reference_fields_module.o \
        climate_module.o \
        SMB_module.o \
        BMB_module.o \
        general_ice_model_data_module.o \
        ice_dynamics_module.o \
        thermodynamics_module.o \
        UFEMISM_main_model.o \
        global_text_output_module.o \
        UFEMISM_program.o

# The executable (the target) depends on the object files:
UFEMISM_program : $(UFEMISM_program_objects:%=$(OBJ_PATH)/%)



# =============================================================
#       Cleaning up everything
# =============================================================

.PHONY : clean

clean :
	rm -fv $(all_programs) $(OBJ_PATH)/*.o $(MOD_PATH)/*.mod *.mod *.lst



# =============================================================
#       Building everything
# =============================================================

.PHONY : all

all : $(all_programs)
