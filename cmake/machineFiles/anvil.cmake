SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicxx CACHE FILEPATH "")

SET (NETCDF_DIR /soft/spack-0.10.0/opt/spack/linux-centos6-x86_64/intel-17.0.1/netcdf-4.4.1.1-mwaopcxxkftrg4jardmcgtw6e7vmxnk4 CACHE FILEPATH "")
SET (PNETCDF_DIR /soft/spack-0.10.0/opt/spack/linux-centos6-x86_64/intel-17.0.1/parallel-netcdf-1.8.0-p7324i4uontmkrbl3hu6sk7buo2zudpy CACHE FILEPATH "")

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
ENDIF ()

