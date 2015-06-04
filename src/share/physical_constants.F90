#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
! This module should 'use' only module 'kinds'
!
module physical_constants
  ! ------------------------------
  use kinds, only : real_kind, longdouble_kind
  ! -----------------------------

  implicit none
  private

  real (kind=real_kind), public, parameter :: DD_PI = 3.141592653589793238462643383279_real_kind
  real (kind=longdouble_kind), public, parameter :: QQ_PI = 3.141592653589793238462643383279_longdouble_kind

  ! to run Held-Suarez in JUPITER mode: 
  ! make rearth,g,omega namelist variables and set to:
  ! omega = 2e-4, rearth=7e7, g=23.0
  real (kind=real_kind), public, parameter :: rearth       = 6.376D6         ! m
  real (kind=real_kind), public, parameter :: g            = 9.80616D0         ! m s^-2
  real (kind=real_kind), public            :: omega        = 7.292D-5        ! radians/s
  real (kind=real_kind), public, parameter :: Rgas         = 287.04D0        
  real (kind=real_kind), public, parameter :: Cp           = 1005.0D0
  real (kind=real_kind), public, parameter :: p0           = 100000.0D0        ! surface pressure (mbar)
  real (kind=real_kind), public, parameter :: MWDAIR       = 28.966D0
  real (kind=real_kind), public, parameter :: Rwater_vapor = 461.50D0
  real (kind=real_kind), public, parameter :: Cpwater_vapor= 1870.0D0 !ASC ANDII VERIFY PLS
  real (kind=real_kind), public, parameter :: kappa        = Rgas/Cp
  real (kind=real_kind), public, parameter :: Rd_on_Rv     = Rgas/Rwater_vapor	
  real (kind=real_kind), public, parameter :: Cpd_on_Cpv     = Cp/Cpwater_vapor
  real (kind=real_kind), public, parameter :: rrearth      = 1.0_real_kind/rearth         ! m

end module physical_constants
