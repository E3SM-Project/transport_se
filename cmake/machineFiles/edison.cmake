# CMake initial cache file for Edison
SET (CMAKE_Fortran_COMPILER ftn CACHE FILEPATH "")
SET (CMAKE_C_COMPILER cc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER CC CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{NETCDF_DIR} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{HDF5_DIR} CACHE FILEPATH "")
SET (PNETCDF_DIR $ENV{PARALLEL_NETCDF_DIR} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{HDF5_DIR} CACHE FILEPATH "")

SET (CMAKE_SYSTEM_NAME Catamount CACHE FILEPATH "")

# override cmake's intel defaults:
# default cmake options for Intel: 
#      -assume byterecl -fp-model precise -ftz -g -O3 
#SET (FORCE_Fortran_FLAGS "-openmp -fp-model fast -ftz  -O3" CACHE STRING "")
#SET (FORCE_Fortran_FLAGS "-openmp -traceback -fp-model precise -ftz -g -O2" CACHE STRING "")

# recommended:
SET (FORCE_Fortran_FLAGS "-openmp -fast -ftz" CACHE STRING "")

# for debugging:
#SET (FORCE_Fortran_FLAGS "-traceback -openmp -ftz -g -O1" CACHE STRING "")