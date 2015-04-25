# - Try to find Netcdf
# Once done this will define
#  NETCDF_FOUND - System has Netcdf
#  NETCDF_INCLUDE_DIRS - The Netcdf include directories
# NETCDF_C_LIBRARIES - The C libraries needed to use Netcdf
# NETCDF_Fortran_LIBRARIES - The Fortran libraries needed to use Netcdf
#  NETCDF_LIBRARIES - All the libraries needed to use Netcdf
#  NETCDF_DEFINITIONS - Compiler switches required for using Netcdf

find_path(NETCDF_INCLUDE_DIR netcdf.h
          HINTS ${NETCDF_DIR}/include )

find_path(NETCDF_LIB_DIR NAMES libnetcdf.a libnetcdf.so
          HINTS ${NETCDF_DIR}/lib ${NETCDF_DIR}/lib64 )

find_path(NETCDF_FORTRAN_LIB_DIR libnetcdff.a libnetcdff.so
          HINTS ${NETCDF_DIR}/lib ${NETCDF_DIR}/lib64 )

find_file(NETCDF4_PAR_H netcdf_par.h 
          HINTS ${NETCDF_INCLUDE_DIR}
          NO_DEFAULT_PATH )

#MESSAGE("PAR_H: ${NETCDF4_PAR_H}")
set(NETCDF_C_LIBRARY "-L${NETCDF_LIB_DIR}  -lnetcdf")

if(NOT NETCDF_FORTRAN_LIB_DIR)
  MESSAGE("WARNING: did not find netcdf fortran library")
else()
  set(NETCDF_Fortran_LIBRARY "-L${NETCDF_LIB_DIR} -lnetcdff")
endif()
set(NETCDF_LIBRARIES "${NETCDF_Fortran_LIBRARY} ${NETCDF_C_LIBRARY}")
if(NOT NETCDF4_PAR_H)
  set(NETCDF4_PARALLEL "no")
  MESSAGE("NETCDF built without MPIIO")
else()
  set(NETCDF4_PARALLEL "yes")
  MESSAGE("NETCDF built with hdf5 MPIIO support")
endif()

set(NETCDF_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR} )

# If the above compilation/link step failed then we need to find dependencies
MESSAGE(STATUS "Determining Netcdf dependencies.")

find_path(Netcdf_NC_CONFIG_BIN
          NAMES nc-config
          HINTS ${NETCDF_INCLUDE_DIR}/../bin ${NETCDF_BIN_DIR}
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

IF (NOT ${Netcdf_NC_CONFIG_BIN} STREQUAL Netcdf_NC_CONFIG_BIN-NOTFOUND)

  # Probe nc-config to determine dependencies of Netcdf
  MESSAGE(STATUS "nc-config found at ${Netcdf_NC_CONFIG_BIN}")

  # use nc-confg --has-nc4 to determine if Netcdf depends upon HDF5
  EXECUTE_PROCESS(COMMAND ${Netcdf_NC_CONFIG_BIN}/nc-config --has-nc4
    RESULT_VARIABLE Homme_result
    OUTPUT_VARIABLE Homme_output
    ERROR_VARIABLE Homme_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  IF (${Homme_error})
    MESSAGE(FATAL_ERROR "${Netcdf_NC_CONFIG_BIN}/nc-config --has-nc4 produced an error")
  ELSE ()
    IF (${Homme_output} STREQUAL yes)
      SET (NETCDF_REQUIRE_HDF5 TRUE)
      MESSAGE(STATUS "nc-config: Netcdf depends upon HDF5")
    ELSE ()
      SET (NETCDF_REQUIRE_HDF5 FALSE)
      MESSAGE(STATUS "nc-config: Netcdf does not depend upon HDF5")
    ENDIF ()
  ENDIF ()

  # use nc-confg --has-dap to determine if Netcdf depends upon CURL
  EXECUTE_PROCESS(COMMAND ${Netcdf_NC_CONFIG_BIN}/nc-config --has-dap
    RESULT_VARIABLE Homme_result
    OUTPUT_VARIABLE Homme_output
    ERROR_VARIABLE Homme_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  IF (${Homme_error})
    MESSAGE(FATAL_ERROR "${Netcdf_NC_CONFIG_BIN}/nc-config --has-dap produced an error")
  ELSE ()
    IF (${Homme_output} STREQUAL yes)
      SET (NETCDF_REQUIRE_CURL TRUE)
      MESSAGE(STATUS "nc-config: Netcdf depends upon CURL")
    ELSE ()
      SET (NETCDF_REQUIRE_CURL FALSE)
      MESSAGE(STATUS "nc-config: Netcdf does not depend upon CURL")
    ENDIF ()
  ENDIF ()

ELSE ()
  SET (NETCDF_REQUIRE_HDF5 TRUE)
  SET (NETCDF_REQUIRE_CURL TRUE)
  MESSAGE(STATUS "nc-config not found assuming hdf5 and curl dependencies")
ENDIF ()

IF (${NETCDF_REQUIRE_CURL})

  # For some reasone CURL uses CURL_ROOT rather than CURL_DIR
  #   - change the variable for consistency
  SET(CURL_ROOT ${CURL_DIR})
  find_package(CURL)

  IF (${CURL_FOUND})
    #MESSAGE(STATUS "Found CURL: ${CURL_LIBRARY}")
    set(HDF5_LIBRARIES ${HDF5_LIBRARIES} ${CURL_LIBRARY})
    get_filename_component(CPRNC_CURL_DIR ${CURL_LIBRARY} PATH)
    get_filename_component(CPRNC_CURL_DIR ${CURL_LIBRARY}/.. ABSOLUTE)
    #MESSAGE(STATUS "CPRNC_CURL_DIR=${CPRNC_CURL_DIR}")
  ELSE ()
    MESSAGE(FATAL_ERROR "CURL Not found")
  ENDIF ()
ENDIF ()

IF (${NETCDF_REQUIRE_HDF5})

  find_path(HDF5_INCLUDE_DIR
            hdf5.h
            PATHS ${HDF5_DIR} ${Homme_HDF5_DIR}
            PATH_SUFFIXES include
            NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  find_library(HDF5_LIBRARY
               NAMES libhdf5.a hdf5
               PATHS ${HDF5_DIR} ${Homme_HDF5_DIR}
               PATH_SUFFIXES lib lib64
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  find_library(HDF5hl_LIBRARY
               NAMES libhdf5_hl.a hdf5_hl
               PATHS ${HDF5_DIR} ${Homme_HDF5_DIR}
               PATH_SUFFIXES lib lib64
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)


  if(${HDF5_LIBRARY} STREQUAL "HDF5_LIBRARY-NOTFOUND" OR ${HDF5hl_LIBRARY} STREQUAL "HDF5hl_LIBRARY-NOTFOUND")
    set(HDF5_FOUND OFF)
    MESSAGE(FATAL_ERROR "HDF5 not found, set HDF5_DIR to appropriate installation")
  else()
    set(HDF5_FOUND ON)
    MESSAGE(STATUS "Found HDF5: ${HDF5hl_LIBRARY} ${HDF5_LIBRARY}")
    set(HDF5_LIBRARIES "-L${HDF5_DIR}/lib -lhdf5_hl -lhdf5 ${HDF5_LIBRARIES}")
    MESSAGE(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")

    #set(NETCDF_LIBRARIES ${NETCDF_LIBRARIES} ${HDF5hl_LIBRARY} ${HDF5_LIBRARY})
    get_filename_component(CPRNC_HDF5_DIR ${HDF5_INCLUDE_DIR}/.. ABSOLUTE)
    #MESSAGE(STATUS "CPRNC_HDF5_DIR=${CPRNC_HDF5_DIR}")
  endif()

  # Check to see which dependencies (ZLIB, SZIP) hdf5 has
  INCLUDE(CheckSymbolExists)
  CHECK_SYMBOL_EXISTS(H5_HAVE_SZLIB_H "${HDF5_INCLUDE_DIR}/H5pubconf.h" HDF5_REQUIRE_SZ)
  CHECK_SYMBOL_EXISTS(H5_HAVE_ZLIB_H "${HDF5_INCLUDE_DIR}/H5pubconf.h" HDF5_REQUIRE_ZLIB)

  IF (${HDF5_REQUIRE_ZLIB})

    MESSAGE(STATUS "This HDF5 library requires ZLIB")
    # Find package always finds the shared object
    #find_package(ZLIB REQUIRED)
    find_library(ZLIB_LIBRARY
                 NAMES libz.a z
                 PATHS ${ZLIB_DIR} ${Homme_ZLIB_DIR}
                 PATH_SUFFIXES lib lib64
                 NO_SYSTEM_ENVIRONMENT_PATH)

    IF(${ZLIB_LIBRARY} STREQUAL "ZLIB_LIBRARY-NOTFOUND")
      SET(ZLIB_FOUND OFF)
      MESSAGE(FATAL_ERROR "ZLIB Not found")
    ELSE()
      SET(ZLIB_FOUND ON)
      get_filename_component(CPRNC_ZLIB_DIR ${ZLIB_LIBRARY} PATH)
      get_filename_component(CPRNC_ZLIB_DIR ${CPRNC_ZLIB_DIR}/.. ABSOLUTE)
      #MESSAGE(STATUS "CPRNC_ZLIB_DIR=${CPRNC_ZLIB_DIR}")
    ENDIF ()

    IF (${ZLIB_FOUND})
      MESSAGE(STATUS "Found ZLIB: ${ZLIB_LIBRARY}")
      set(HDF5_LIBRARIES "${HDF5_LIBRARIES} -L${ZLIB_DIR}/lib -lz")
    ELSE ()
      MESSAGE(FATAL_ERROR "ZLIB Not found, set ZLIB_DIR to appropriate installation")
    ENDIF ()
  ENDIF ()

  IF (${HDF5_REQUIRE_SZ})
    MESSAGE(STATUS "This HDF5 library requires SZIP")

    find_library(SZIP_LIBRARY
                 NAMES libsz.a sz
                 PATHS ${SZIP_DIR} ${Homme_SZIP_DIR}
                 PATH_SUFFIXES lib lib64
                 NO_SYSTEM_ENVIRONMENT_PATH)

    IF(${SZIP_LIBRARY} STREQUAL "SZIP_LIBRARY-NOTFOUND")
      SET(SZIP_FOUND OFF)
      MESSAGE(FATAL_ERROR "SZIP Not found, set SZIP_DIR to appropriate installation")
    ELSE()
      SET(SZIP_FOUND ON)
      MESSAGE(STATUS "Found SZIP: ${SZIP_LIBRARY}")
      SET(HDF5_LIBRARIES "${HDF5_LIBRARIES} -L${SZIP_DIR/lib}/lib -lsz")
      get_filename_component(CPRNC_SZIP_DIR ${SZIP_LIBRARY} PATH)
      get_filename_component(CPRNC_SZIP_DIR ${CPRNC_SZIP_DIR}/.. ABSOLUTE)
      #MESSAGE(STATUS "CPRNC_SZIP_DIR=${CPRNC_SZIP_DIR}")
    ENDIF()
  ENDIF ()

ENDIF()

if(${HDF5_FOUND}) 
  MESSAGE(STATUS "Adding hdf5 libraries ${HDF5_LIBRARIES}")
  set(NETCDF_LIBRARIES ${NETCDF_LIBRARIES} ${HDF5_LIBRARIES})
endif()

# Export variables so other projects can use them as well
#  ie. if pio is added with add_subdirectory
SET(NETCDF_INCLUDE_DIR ${NETCDF_INCLUDE_DIR} CACHE STRING "Location of NetCDF include files.")
SET(NETCDF_LIBRARIES ${NETCDF_LIBRARIES} CACHE STRING "Link line for NetCDF.")

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(NETCDF  DEFAULT_MSG NETCDF_LIBRARIES
                                  NETCDF_C_LIBRARY NETCDF_Fortran_LIBRARY NETCDF_INCLUDE_DIR)

mark_as_advanced(NETCDF_INCLUDE_DIR NETCDF_LIBRARIES NETCDF_C_LIBRARY NETCDF_Fortran_LIBRARY NETCDF4_PARALLEL )
