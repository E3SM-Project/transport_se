# CMake initial cache file for Edison
SET (CMAKE_Fortran_COMPILER ftn CACHE FILEPATH "")
SET (CMAKE_C_COMPILER cc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER CC CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{NETCDF_DIR} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{HDF5_DIR} CACHE FILEPATH "")
SET (PNETCDF_DIR $ENV{PARALLEL_NETCDF_DIR} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{HDF5_DIR} CACHE FILEPATH "")

SET (CMAKE_SYSTEM_NAME Catamount CACHE FILEPATH "")
SET (ADD_Fortran_FLAGS "-fast" CACHE STRING "")

# flags set by cmake script (with OMP enabled):
# intel: Fortran Flags =  -assume byterecl -fp-model precise -ftz -O3 -openmp -fast

# for debugging:
# SET (FORCE_Fortran_FLAGS "-traceback -fbounds-check -openmp -ftz -g -O1" CACHE STRING "")