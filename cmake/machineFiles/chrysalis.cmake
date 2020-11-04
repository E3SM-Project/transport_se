SET (CMAKE_Fortran_COMPILER mpiifort CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpiicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpiicpc CACHE FILEPATH "")

SET (NetCDF_C_PATH /nfs/spack-0.15.4/opt/spack/linux-centos8-x86_64/intel-19.0.5/netcdf-c-4.7.4-xvlcfbf CACHE FILEPATH "")
SET (NetCDF_Fortran_PATH /nfs/spack-0.15.4/opt/spack/linux-centos8-x86_64/intel-19.0.5/netcdf-fortran-4.5.3-25tcwku CACHE FILEPATH "")
SET (PNETCDF_DIR /nfs/spack-0.15.4/opt/spack/linux-centos8-x86_64/intel-19.0.5/parallel-netcdf-1.12.1-2ixlykf CACHE FILEPATH "")
SET (HDF5_DIR /nfs/spack-0.15.4/opt/spack/linux-centos8-x86_64/intel-19.0.5/hdf5-1.10.6-cbztx4u CACHE FILEPATH "")
SET (SZIP_DIR /nfs/spack-0.15.4/opt/spack/linux-centos8-x86_64/intel-19.0.5/libszip-2.1.1-mjrvf4g CACHE FILEPATH "")

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
  SET (ADD_LINKER_FLAGS " ${ADD_LINKER_FLAGS} ${NFCONFIG_OUTPUT} " CACHE STRING "")
  MESSAGE(STATUS "Added ${NFCONFIG_OUTPUT} to link flags: ${ADD_LINKER_FLAGS}")
ENDIF ()

