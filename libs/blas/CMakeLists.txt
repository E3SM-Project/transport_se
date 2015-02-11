
# This is to fix an issue with the IBM compiler wrappers
IF (${CMAKE_Fortran_COMPILER_ID} MATCHES XL)
  ADD_DEFINITIONS(-qfixed)
ENDIF()

IF (${CMAKE_Fortran_COMPILER_ID} MATCHES NAG)
  ADD_DEFINITIONS(-dcfuns)
ENDIF()

SET(HOMME_BLASSRC
  dasum.f
  daxpy.f
  dcopy.f
  ddot.f 
  dgemm.f
  dgemv.f
  dger.f
  dnrm2.f
  drot.f 
  drotg.f
  drotm.f
  dscal.f
  dswap.f
  dtrmm.f
  dtrmv.f
  dsyr2k.f
  dsyr2.f
  dsymv.f
  idamax.f
  lsame.f
  xerbla.f
  dznrm2.f
  zscal.f
  sswap.f
  sdot.f
)

ADD_LIBRARY(blas ${HOMME_BLASSRC})

