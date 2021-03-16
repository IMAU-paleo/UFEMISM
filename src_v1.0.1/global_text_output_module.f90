MODULE global_text_output_module

  USE data_types_module,           ONLY: type_model_region, type_mesh, type_ice_model, type_PD_data_fields, type_init_data_fields, &
                                         type_climate_model, type_climate_matrix, type_SMB_model, type_BMB_model
  USE configuration_module,        ONLY: dp, C

  IMPLICIT NONE

CONTAINS
  
  SUBROUTINE create_text_output_file
    
    IMPLICIT NONE
    
    CHARACTER(LEN=256)                            :: filename
    
    filename = TRIM(C%output_dir) // 'aa_general_output.txt'
    OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
    
    WRITE(UNIT = 1337, FMT = '(A)') 'UFEMISM global output data'
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A)') '     Time     sealevel     CO2      d18O_obs   d18O_mod   d18O_ice   d18O_Tdw     sl_NAM     sl_EAS     sl_GRL     sl_ANT'
    
    CLOSE(UNIT = 1337)
    
  END SUBROUTINE create_text_output_file
  SUBROUTINE write_text_output( time, sl_glob, CO2, d18O_obs, d18O_mod, d18O_ice, d18O_Tdw, sl_NAM, sl_EAS, sl_GRL, sl_ANT)
    ! Write data to global output file
  
    IMPLICIT NONE  
    
    REAL(dp),                   INTENT(IN)        :: time, sl_glob, CO2, d18O_obs, d18O_mod, d18O_ice, d18O_Tdw, sl_NAM, sl_EAS, sl_GRL, sl_ANT
    CHARACTER(LEN=256)                            :: filename
            
    filename = TRIM(C%output_dir) // 'aa_general_output.txt'
    OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
    
    WRITE(UNIT = 1337, FMT = '(F10.1,10F11.2)') time, sl_glob, CO2, d18O_obs, d18O_mod, d18O_ice, d18O_Tdw, sl_NAM, sl_EAS, sl_GRL, sl_ANT
    
    CLOSE(UNIT = 1337)
    
  END SUBROUTINE write_text_output

END MODULE global_text_output_module
