# Relevant Directories within the HOMME repository
# holding source files needed for building the library
# and executable preqx.



SET(UTILS_SHARE_DIR ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR         ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR   ${HOMME_SOURCE_DIR}/src/share)
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/timing)

# Make INCLUDE_DIRS global so the tests can access it
SET (PREQX_INCLUDE_DIRS ${PIO_INCLUDE_DIRS}
                        ${TIMING_INCLUDE_DIRS}
                        ${SRC_SHARE_DIR} CACHE INTERNAL "include directories for preqx")

SET(PREQX_SRCS_F90
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_SHARE_DIR}/dimensions_mod.F90
  ${SRC_SHARE_DIR}/control_mod.F90
  ${SRC_SHARE_DIR}/physical_constants.F90
  ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
  ${SRC_SHARE_DIR}/params_mod.F90
  ${SRC_SHARE_DIR}/gridgraph_mod.F90
  ${SRC_SHARE_DIR}/edge_mod.F90
  ${SRC_SHARE_DIR}/cuda_mod.F90
  ${SRC_SHARE_DIR}/element_mod.F90
  ${SRC_SHARE_DIR}/parallel_mod.F90
  ${SRC_SHARE_DIR}/hybrid_mod.F90
  ${SRC_SHARE_DIR}/reduction_mod.F90
  ${SRC_SHARE_DIR}/quadrature_mod.F90
  ${SRC_SHARE_DIR}/global_norms_mod.F90
  ${SRC_SHARE_DIR}/ll_mod.F90
  ${SRC_SHARE_DIR}/metagraph_mod.F90
  ${SRC_SHARE_DIR}/schedule_mod.F90
  ${SRC_SHARE_DIR}/schedtype_mod.F90 
  ${SRC_SHARE_DIR}/thread_mod.F90
  ${SRC_SHARE_DIR}/bndry_mod.F90
  ${SRC_SHARE_DIR}/derivative_mod.F90
  ${SRC_SHARE_DIR}/viscosity_mod.F90
  ${SRC_SHARE_DIR}/spacecurve_mod.F90
  ${SRC_SHARE_DIR}/cube_mod.F90
  ${SRC_SHARE_DIR}/dof_mod.F90
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_DIR}/repro_sum_mod.F90 
  ${SRC_SHARE_DIR}/time_mod.F90
  ${SRC_SHARE_DIR}/hybvcoord_mod.F90
  ${UTILS_SHARE_DIR}/shr_kind_mod.F90
  ${UTILS_SHARE_DIR}/shr_mpi_mod.F90 
  ${UTILS_SHARE_DIR}/shr_sys_mod.F90 
  ${UTILS_SHARE_DIR}/shr_file_mod.F90 
  ${UTILS_SHARE_DIR}/shr_vmath_mod.F90 
  ${UTILS_SHARE_DIR}/shr_const_mod.F90 
  ${UTILS_SHARE_DIR}/shr_spfn_mod.F90 
  ${UTILS_TIMING_DIR}/perf_mod.F90 
  ${UTILS_TIMING_DIR}/perf_utils.F90 
  ${SRC_SHARE_DIR}/domain_mod.F90 
  ${SRC_DIR}/restart_io_mod.F90 
  ${SRC_SHARE_DIR}/filter_mod.F90
  ${SRC_DIR}/common_movie_mod.F90 
  ${SRC_SHARE_DIR}/namelist_mod.F90 
  ${SRC_SHARE_DIR}/mass_matrix_mod.F90 
  ${SRC_DIR}/checksum_mod.F90 
  ${SRC_DIR}/baroclinic_inst_mod.F90 
  ${SRC_SHARE_DIR}/diffusion_mod.F90 
  ${SRC_DIR}/asp_tests.F90 
  ${SRC_SHARE_DIR}/prim_state_mod.F90 
  ${SRC_DIR}/prim_restart_mod.F90 
  ${SRC_SHARE_DIR}/prim_advance_mod.F90
  ${SRC_SHARE_DIR}/prim_advection_mod.F90
  ${SRC_SHARE_DIR}/metis_mod.F90
  ${SRC_SHARE_DIR}/prim_driver_mod.F90
  ${SRC_DIR}/prim_movie_mod.F90
  ${SRC_DIR}/netcdf_io_mod.F90 
  ${SRC_SHARE_DIR}/mesh_mod.F90 
  ${SRC_DIR}/prim_main.F90 
  ${SRC_SHARE_DIR}/dcmip_123_mod.F90
  ${SRC_SHARE_DIR}/dcmip_wrapper_mod.F90
)

#SET(PREQX_SRCS_C
#  ${UTILS_SHARE_DIR}/shr_vmath_fwrap.c
#  ${SRC_DIR}/jrio.c
#)

# If the user specified a file for custom compiler options use those
IF (DEFINED PREQX_CUSTOM_FLAGS_FILE)
  setCustomCompilerFlags(PREQX_CUSTOM_FLAGS_FILE PREQX_SRCS_F90)
ENDIF ()

# Make SRCS global so the tests can access it
SET(PREQX_SRCS ${PREQX_SRCS} ${PREQX_SRCS_C} ${PREQX_SRCS_F90} CACHE INTERNAL "Sources for preqx")

# Set up defaults
IF (NOT PREQX_NP)
  SET (PREQX_NP 4) 
ENDIF ()
IF (NOT PREQX_PLEV)
  SET (PREQX_PLEV 20) 
ENDIF ()
IF (NOT PREQX_USE_PIO)
  SET (PREQX_USE_PIO FALSE) 
ENDIF ()
IF (NOT PREQX_USE_ENERGY)
  SET (PREQX_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT PREQX_NC)
  SET (PREQX_NC 4)
ENDIF ()
IF (PREQX_QSIZE_D)
  SET (QSIZE_D ${PREQX_QSIZE_D})
ENDIF ()

MESSAGE(STATUS "Building preqx with:")
MESSAGE(STATUS "  NP = ${PREQX_NP}")
MESSAGE(STATUS "  NC = ${PREQX_NC}")
MESSAGE(STATUS "  PLEV = ${PREQX_PLEV}")
MESSAGE(STATUS "  PIO = ${PREQX_USE_PIO}")
MESSAGE(STATUS "  ENERGY = ${PREQX_USE_ENERGY}")

IF (${PREQX_USE_CUDA_FORTRAN})
  SET(USE_CUDA_FORTRAN TRUE)
  MESSAGE(STATUS "  Using CUDA Fortran")
ENDIF ()

IF (${PREQX_USE_OPENACC})
  SET(USE_OPENACC TRUE)
  MESSAGE(STATUS "  Using OpenACC")
ENDIF ()

IF (${ENABLE_HORIZ_OPENMP} AND ${ENABLE_COLUMN_OPENMP})
  SET(NESTED_OPENMP TRUE CACHE BOOL "Both horizontal and column threading")
ENDIF ()

# The following macro handles the compilation with the appropriate variables
############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(preqx preqx ${PREQX_NP} ${PREQX_NC} ${PREQX_PLEV} 
               ${PREQX_USE_PIO}  ${PREQX_USE_ENERGY})


