##############################################################################
# Compiler specific options
##############################################################################
MESSAGE(STATUS "CMAKE_Fortran_COMPILER_ID = ${CMAKE_Fortran_COMPILER_ID}")
# Need this for a fix in repro_sum_mod
IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL XL)
  ADD_DEFINITIONS(-DnoI8)
ENDIF ()

IF (DEFINED BASE_FFLAGS)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${BASE_FFLAGS}")
ELSE ()
  IF (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mextend -Kieee -Mflushz")
    # Needed by csm_share
    ADD_DEFINITIONS(-DCPRPGI)
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL PathScale)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -extend-source")
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume byterecl")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model precise -ftz")
    # Needed by csm_share
    ADD_DEFINITIONS(-DCPRINTEL)
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL XL)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -WF,-C! -qstrict -qnosave")
    # Needed by csm_share
    ADD_DEFINITIONS(-DCPRIBM)
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -kind=byte -wmismatch=mpi_send,mpi_recv,mpi_bcast,mpi_allreduce,mpi_reduce,mpi_isend,mpi_irecv,mpi_irsend,mpi_rsend,mpi_gatherv,mpi_gather,mpi_scatterv,mpi_allgather,mpi_alltoallv,mpi_file_read_all,mpi_file_write_all,mpi_file_read_at")
#    SET(OPT_FFLAGS "${OPT_FFLAGS} -ieee=full -O2")
    SET(DEBUG_FFLAGS "${DEBUG_FFLAGS} -g -time -f2003 -ieee=stop")
    ADD_DEFINITIONS(-DHAVE_F2003_PTR_BND_REMAP)
    # Needed by both PIO and csm_share
    ADD_DEFINITIONS(-DCPRNAG)
  ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DHAVE_F2003_PTR_BND_REMAP")
    # Needed by csm_share
    ADD_DEFINITIONS(-DCPRCRAY)
 ENDIF ()
ENDIF ()

##############################################################################
# Optimization flags
# 1) OPT_FLAGS if specified sets the Fortran,C, and CXX optimization flags
# 2) OPT_FFLAGS if specified sets the Fortran optimization flags
# 3) OPT_CFLAGS if specified sets the C optimization flags
# 4) OPT_CXXFLAGS if specified sets the CXX optimization flags
##############################################################################
IF (OPT_FLAGS)
  # Flags for Fortran C and CXX
  SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPT_FLAGS}")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPT_FLAGS}")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_FLAGS}")

ELSE ()

  IF (OPT_FFLAGS)
    # User specified optimization flags
    SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPT_FFLAGS}")
  ELSE ()
    # Defaults
    IF (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL PGI)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL PathScale)
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
      #SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mavx -DTEMP_INTEL_COMPILER_WORKAROUND_001")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL XL)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
    ELSEIF (CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
    ENDIF ()
  ENDIF ()
 
  IF (OPT_CFLAGS)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPT_CFLAGS}")
  ELSE ()
    IF (CMAKE_C_COMPILER_ID STREQUAL GNU)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL PGI)
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL PathScale)
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL Intel)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
      #SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx -DTEMP_INTEL_COMPILER_WORKAROUND_001")
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL XL)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
    ELSEIF (CMAKE_C_COMPILER_ID STREQUAL Cray)
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
    ENDIF ()
  ENDIF ()

  IF (OPT_CXXFLAGS)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_CXXFLAGS}")
  ELSE ()
    IF (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL PGI)
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL PathScale)
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL Intel)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
      #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -DTEMP_INTEL_COMPILER_WORKAROUND_001")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL XL)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
    ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL Cray)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
    ENDIF ()
  ENDIF ()
ENDIF ()

##############################################################################
# DEBUG flags
# 1) DEBUG_FLAGS if specified sets the Fortran,C, and CXX debug flags
# 2) DEBUG_FFLAGS if specified sets the Fortran debugflags
# 3) DEBUG_CFLAGS if specified sets the C debug flags
# 4) DEBUG_CXXFLAGS if specified sets the CXX debug flags
##############################################################################
IF (DEBUG_FLAGS)
  SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${DEBUG_FLAGS}")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${DEBUG_FLAGS}")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEBUG_FLAGS}")
ELSE ()
  IF (DEBUG_FFLAGS)
    SET (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${DEBUG_FFLAGS}")
  ELSE ()
    IF(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray)
      SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    ENDIF ()
  ENDIF ()

  IF (DEBUG_CFLAGS)
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${DEBUG_CFLAGS}")
  ELSE ()
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  ENDIF ()

  IF (DEBUG_CXXFLAGS)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEBUG_CXXFLAGS}")
  ELSE ()
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  ENDIF ()

ENDIF ()


##############################################################################
# OpenMP
# Two flavors:
#   1) HORIZ_OPENMP OpenMP over elements (standard OPENMP)
#   2) COLUMN_OPENMP OpenMP within an element (previously called ELEMENT_OPENMP)
##############################################################################
OPTION(ENABLE_OPENMP "OpenMP across elements" TRUE)
OPTION(ENABLE_HORIZ_OPENMP "OpenMP across elements" TRUE)
OPTION(ENABLE_COLUMN_OPENMP "OpenMP within an element" TRUE)

# If OpenMP is turned off also turn off ENABLE_HORIZ_OPENMP
IF (NOT ${ENABLE_OPENMP}) 
  SET(ENABLE_HORIZ_OPENMP FALSE) 
  SET(ENABLE_COLUMN_OPENMP FALSE) 
ENDIF ()

##############################################################################
IF (ENABLE_HORIZ_OPENMP OR ENABLE_COLUMN_OPENMP)
  IF(NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray) 
    FIND_PACKAGE(OpenMP)
    IF(OPENMP_FOUND)
      MESSAGE(STATUS "Found OpenMP Flags")
      # This file is needed for the timing library - this is currently
      # inaccessible from the timing CMake script
      IF (CMAKE_Fortran_COMPILER_ID STREQUAL XL)
        SET(OpenMP_C_FLAGS "-qsmp=omp")
        IF (ENABLE_COLUMN_OPENMP)
          SET(OpenMP_C_FLAGS "-qsmp=omp:nested_par")
        ENDIF ()
      ENDIF ()
      SET(OpenMP_Fortran_FLAGS "${OpenMP_C_FLAGS}")
      MESSAGE(STATUS "OpenMP_Fortran_FLAGS: ${OpenMP_Fortran_FLAGS}")
      MESSAGE(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
      MESSAGE(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
      MESSAGE(STATUS "OpenMP_EXE_LINKER_FLAGS: ${OpenMP_EXE_LINKER_FLAGS}")
      # The fortran openmp flag should be the same as the C Flag
      SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
      SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    ELSE ()
      MESSAGE(FATAL_ERROR "Unable to find OpenMP")
    ENDIF()
  ENDIF()
 IF (${ENABLE_HORIZ_OPENMP})
   # Set this as global so it can be picked up by all executables
   SET(HORIZ_OPENMP TRUE CACHE BOOL "Threading in the horizontal direction")
   MESSAGE(STATUS "  Using HORIZ_OPENMP")
 ENDIF ()
 
 IF (${ENABLE_COLUMN_OPENMP})
   # Set this as global so it can be picked up by all executables
   SET(COLUMN_OPENMP TRUE CACHE BOOL "Threading in the horizontal direction")
   MESSAGE(STATUS "  Using COLUMN_OPENMP")
 ENDIF ()
ENDIF ()
##############################################################################

##############################################################################
# OpenACC specific flag - only supporting PGI compiler
##############################################################################
OPTION(ENABLE_OPENACC "Whether to build with OpenACC support" FALSE)
IF (${ENABLE_OPENACC})
  IF (${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
    # Need to add -acc to the Fortran FLAGS to see if it will compile 
    # "call acc_init()"
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc")
    TRY_COMPILE(OPENACC_SUCCESS
                ${CMAKE_BINARY_DIR}/tests/compilerTests/
                ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilerTests/openAccTest.f90
                OUTPUT_VARIABLE COMPILE_OUTPUT)
    IF (${OPENACC_SUCCESS})
      MESSAGE(STATUS "Using OpenACC through PGI compiler")
      SET(PREQX_USE_OPENACC TRUE)
    ELSE()
      MESSAGE(FATAL_ERROR "Unable to use OpenACC through the PGI compiler")
    ENDIF()
  ELSEIF (${CMAKE_Fortran_COMPILER_ID} STREQUAL Cray)
    SET(PREQX_USE_OPENACC TRUE)
    MESSAGE(STATUS "Using OpenACC through Cray compiler")
  ELSE ()
    MESSAGE(FATAL_ERROR "OpenACC only supported through the PGI and Cray compilers")
  ENDIF()
ENDIF ()

##############################################################################
# CUDA Fortran specific flags - only supporting PGI compiler
##############################################################################
OPTION(ENABLE_CUDA_FORTRAN "Whether to build with CUDA Fortran support" FALSE)
IF (${ENABLE_CUDA_FORTRAN})
  IF (NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL PGI)
    MESSAGE(FATAL_ERROR "CUDA Fortran only supported through the PGI compiler")
  ELSE ()
    # Set PGI CUDA Fortran flags

    # Set defaults as lowest version of CUDA and device capability allowed
    # To do: determine a way to generalize this
    IF (NOT CUDA_VERSION)
      SET(CUDA_VERSION "4.1")
    ENDIF ()

    # Compute capability: cc2x is for devices with compute capability >= 2.0 
    IF (NOT CUDA_DEVICE_CAPABILITY)
      SET(CUDA_DEVICE_CAPABILITY "cc2x")
    ENDIF ()

    SET(CMAKE_Fortran_FLAGS 
        "${CMAKE_Fortran_FLAGS} -ta=nvidia -Mcuda=${CUDA_VERSION},${CUDA_DEVICE_CAPABILITY},ptxinfo,keepgpu")

    MESSAGE(STATUS "Testing PGI CUDA Fortran Compilation with flags: ${CMAKE_Fortran_FLAGS}")

    TRY_COMPILE(CUDAFOR
                ${CMAKE_BINARY_DIR}/tests/compilerTests/
                ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilerTests/cudaFortranTest.f90
                OUTPUT_VARIABLE COMPILE_OUTPUT)
    IF (${CUDAFOR})
      SET(PREQX_USE_CUDA_FORTRAN TRUE)
      MESSAGE(STATUS "Succeeded. Using CUDA Fortran through PGI compiler")
    ELSE()
      SET(PREQX_USE_CUDA_FORTRAN FALSE)
      MESSAGE(FATAL_ERROR "Unable to use CUDA Fortran through the PGI "
              "compiler. Compilation failed with the following "
              "output.\n${COMPILE_OUTPUT}")
    ENDIF()
  ENDIF()
ENDIF ()

##############################################################################
# Intel Phi (MIC) specific flags - only supporting the Intel compiler
##############################################################################
OPTION(ENABLE_INTEL_PHI "Whether to build with Intel Xeon Phi (MIC) support" FALSE)
IF (ENABLE_INTEL_PHI)
  IF (NOT ${CMAKE_Fortran_COMPILER_ID} STREQUAL Intel)
    MESSAGE(FATAL_ERROR "Intel Phi acceleration only supported through the Intel compiler")
  ELSE ()
    SET(INTEL_PHI_FLAGS "-mmic")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${INTEL_PHI_FLAGS}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  ${INTEL_PHI_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INTEL_PHI_FLAGS}")
    SET(IS_ACCELERATOR TRUE)
    # CMake magic for cross-compilation
    SET(CMAKE_SYSTEM_NAME Linux) 
    SET(CMAKE_SYSTEM_PROCESSOR k1om) 
    SET(CMAKE_SYSTEM_VERSION 1) 
    SET(_CMAKE_TOOLCHAIN_PREFIX  x86_64-k1om-linux-)
    # Specify the location of the target environment
    IF (TARGET_ROOT_PATH)
      SET(CMAKE_FIND_ROOT_PATH ${TARGET_ROOT_PATH})
    ELSE ()
      SET(CMAKE_FIND_ROOT_PATH /usr/linux-k1om-4.7)
    ENDIF ()
  ENDIF ()
ENDIF ()

##############################################################################
# Allow the option to add compiler flags to those provided
##############################################################################
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ADD_Fortran_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ADD_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADD_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ADD_LINKER_FLAGS}")

##############################################################################
# Allow the option to override compiler flags
##############################################################################
IF (FORCE_Fortran_FLAGS)
  SET(CMAKE_Fortran_FLAGS ${FORCE_Fortran_FLAGS})
ENDIF ()
IF (FORCE_C_FLAGS)
  SET(CMAKE_C_FLAGS ${FORCE_C_FLAGS})
ENDIF ()
IF (FORCE_CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${FORCE_CXX_FLAGS})
ENDIF ()
IF (FORCE_LINKER_FLAGS)
  SET(CMAKE_EXE_LINKER_FLAGS ${FORCE_LINKER_FLAGS})
ENDIF ()


