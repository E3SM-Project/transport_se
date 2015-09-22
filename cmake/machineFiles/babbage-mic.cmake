# CMake initial cache file for babbage host 
SET (CMAKE_Fortran_COMPILER /global/babbage/nsg/opt/intel/impi/5.0.2.044/intel64/bin/mpiifort CACHE FILEPATH "")
SET (CMAKE_C_COMPILER /global/babbage/nsg/opt/intel/impi/5.0.2.044/intel64/bin/mpiicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER /global/babbage/nsg/opt/intel/impi/5.0.2.044/intel64/bin/mpiicpc CACHE FILEPATH "")

SET (HDF5_DIR /global/homes/a/azamat/babbage/soft/hdf5/1.8.14 CACHE FILEPATH "")
SET (NETCDF_DIR /global/homes/a/azamat/babbage/soft/netcdf/4.4.2 CACHE FILEPATH "")
SET (WITH_PNETCDF FALSE CACHE BOOL "")

SET (ENABLE_INTEL_PHI TRUE CACHE BOOL "")
