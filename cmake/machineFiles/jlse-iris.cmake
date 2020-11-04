SET (CMAKE_Fortran_COMPILER mpiifx CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpiicx CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpiicpx CACHE FILEPATH "")

SET (NETCDF_DIR /home/azamat/soft/netcdf/4.4.1c-4.2cxx-4.4.4f/ifx CACHE FILEPATH "")
SET (PNETCDF_DIR /home/azamat/soft/pnetcdf/1.12.1/ifx CACHE FILEPATH "")

EXECUTE_PROCESS(COMMAND ${NETCDF_DIR}/bin/nf-config --flibs
  RESULT_VARIABLE NFCONFIG_RESULT
  OUTPUT_VARIABLE NFCONFIG_OUTPUT
  ERROR_VARIABLE  NFCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)  
IF (${NFCONFIG_ERROR})
  MESSAGE(WARNING "${NETCDF_DIR}/bin/nf-config --flibs produced an error. Default linking will be used.")
ELSE ()  
  SET (ADD_LINKER_FLAGS " ${NFCONFIG_OUTPUT} " CACHE STRING "")
  MESSAGE(STATUS "Added ${NFCONFIG_OUTPUT} to link flags: ${ADD_LINKER_FLAGS}")
ENDIF ()
