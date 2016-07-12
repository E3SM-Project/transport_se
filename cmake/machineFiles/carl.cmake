# CMake machine file for Carl
SET (CMAKE_Fortran_COMPILER mpiifort CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpiicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpiicpc CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{NETCDF_DIR} CACHE FILEPATH "")
SET (CURL_DIR $ENV{CURL_DIR} CACHE FILEPATH "")

SET (ADD_Fortran_FLAGS "-xMIC-AVX512 $ENV{CURL}" CACHE STRING "")
SET (ADD_C_FLAGS "-xMIC-AVX512 $ENV{CURL}" CACHE STRING "")
SET (ADD_CXX_FLAGS "-xMIC-AVX512 $ENV{CURL}" CACHE STRING "")

# flags set by cmake script (with OMP enabled):
# intel: Fortran Flags =  -assume byterecl -fp-model precise -ftz -O3 -openmp -fast

# for debugging:
# SET (FORCE_Fortran_FLAGS "-traceback -fbounds-check -openmp -ftz -g -O1" CACHE STRING "")
