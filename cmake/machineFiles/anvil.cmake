SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicxx CACHE FILEPATH "")

SET (NetCDF_C_PATH /blues/gpfs/software/centos7/spack-latest/opt/spack/linux-centos7-x86_64/intel-17.0.0/netcdf-4.4.1-tckdgwl CACHE FILEPATH "")
SET (NetCDF_Fortran_PATH /blues/gpfs/software/centos7/spack-latest/opt/spack/linux-centos7-x86_64/intel-17.0.0/netcdf-fortran-4.4.4-urmb6ss CACHE FILEPATH "")
SET (PNETCDF_DIR /blues/gpfs/software/centos7/spack-latest/opt/spack/linux-centos7-x86_64/intel-17.0.0/parallel-netcdf-1.11.0-6qz7skn CACHE FILEPATH "")

EXECUTE_PROCESS(COMMAND ${NetCDF_Fortran_PATH}/bin/nf-config --flibs
  RESULT_VARIABLE NFCONFIG_RESULT
  OUTPUT_VARIABLE NFCONFIG_OUTPUT
  ERROR_VARIABLE  NFCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
IF (${NFCONFIG_ERROR})
  MESSAGE(WARNING "${NetCDF_Fortran_PATH}/bin/nf-config --flibs produced an error. Default linking will be used.")
ELSE ()
  SET (ADD_LINKER_FLAGS " ${NFCONFIG_OUTPUT} " CACHE STRING "")
  MESSAGE(STATUS "Added ${NFCONFIG_OUTPUT} to link flags: ${ADD_LINKER_FLAGS}")
ENDIF ()

EXECUTE_PROCESS(COMMAND ${NetCDF_C_PATH}/bin/nc-config --libs
  RESULT_VARIABLE NFCONFIG_RESULT
  OUTPUT_VARIABLE NFCONFIG_OUTPUT
  ERROR_VARIABLE  NFCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)  
IF (${NFCONFIG_ERROR})
  MESSAGE(WARNING "${NetCDF_C_PATH}/bin/nc-config --libs produced an error. Default linking will be used.")
ELSE ()  
  SET (ADD_LINKER_FLAGS " ${NFCONFIG_OUTPUT} " CACHE STRING "")
  MESSAGE(STATUS "Added ${NFCONFIG_OUTPUT} to link flags: ${ADD_LINKER_FLAGS}")
ENDIF ()

