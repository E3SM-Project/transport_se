#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_si_mod
  implicit none
  private

  public :: preq_hydrostatic
  public :: preq_pressure

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  CCM3 hydrostatic integral
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine preq_hydrostatic(phi,phis,T_v,p,dp)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev
    use physical_constants, only : rgas
    implicit none

    !------------------------------Arguments---------------------------------------------------------------
    real(kind=real_kind), intent(out) :: phi(np,np,nlev)     
    real(kind=real_kind), intent(in) :: phis(np,np)
    real(kind=real_kind), intent(in) :: T_v(np,np,nlev)
    real(kind=real_kind), intent(in) :: p(np,np,nlev)   
    real(kind=real_kind), intent(in) :: dp(np,np,nlev)  
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=real_kind) Hkk,Hkl          ! diagonal term of energy conversion matrix
    real(kind=real_kind), dimension(np,np,nlev) :: phii       ! Geopotential at interfaces
    !-----------------------------------------------------------------------

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,j,i,hkk,hkl)
#endif
       do j=1,np   !   Loop inversion (AAM)

          do i=1,np
             hkk = dp(i,j,nlev)*0.5d0/p(i,j,nlev)
             hkl = 2*hkk
             phii(i,j,nlev)  = Rgas*T_v(i,j,nlev)*hkl
             phi(i,j,nlev) = phis(i,j) + Rgas*T_v(i,j,nlev)*hkk 
          end do

          do k=nlev-1,2,-1
             do i=1,np
                ! hkk = dp*ckk
                hkk = dp(i,j,k)*0.5d0/p(i,j,k)
                hkl = 2*hkk
                phii(i,j,k) = phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkl
                phi(i,j,k) = phis(i,j) + phii(i,j,k+1) + Rgas*T_v(i,j,k)*hkk
             end do
          end do

          do i=1,np
             ! hkk = dp*ckk
             hkk = 0.5d0*dp(i,j,1)/p(i,j,1)
             phi(i,j,1) = phis(i,j) + phii(i,j,2) + Rgas*T_v(i,j,1)*hkk
          end do

       end do


end subroutine preq_hydrostatic

!----------------------------------------------------------------------- 
! preq_pressure:
!
! Purpose: 
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure. Originally plevs0!
! 
! Method: 
! 
! Author: B. Boville/ Adapted for HOMME by Rich Loft
! 
!-----------------------------------------------------------------------
!
! $Id: prim_si_mod.F90,v 2.10 2005/10/14 20:17:22 jedwards Exp $
! $Author: jedwards $
!
!-----------------------------------------------------------------------

  subroutine preq_pressure (ps0,  ps,               &
       hyai, hybi, hyam, hybm, &
       pint, pmid, pdel)
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nlevp
    implicit none

    !-----------------------------------------------------------------------

    real(kind=real_kind), intent(in)  :: ps0                ! Hybrid coordinate reference pressure (pascals)
    real(kind=real_kind), intent(in)  :: ps(np,np)          ! Surface pressure (pascals)
    real(kind=real_kind), intent(in)  :: hyai(nlevp)        ! Hybrid interface A coefficients
    real(kind=real_kind), intent(in)  :: hybi(nlevp)        ! Hybrid interface B coefficients
    real(kind=real_kind), intent(in)  :: hyam(nlev)         ! Hybrid midpoint  A coefficients
    real(kind=real_kind), intent(in)  :: hybm(nlev)         ! Hybrid midpoint  B coefficients
    real(kind=real_kind), intent(out) :: pint(np,np,nlevp)  ! Pressure at model interfaces
    real(kind=real_kind), intent(out) :: pmid(np,np,nlev)   ! Pressure at model levels
    real(kind=real_kind), intent(out) :: pdel(np,np,nlev)   ! Layer thickness (pint(k+1) - pint(k))
    !-----------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k             ! Horizontal, level indices
    !-----------------------------------------------------------------------
    !
    ! Set interface pressures
    !
    do k=1,nlevp
       do j=1,np
          do i=1,np
             pint(i,j,k) = hyai(k)*ps0 + hybi(k)*ps(i,j)
          end do
       end do
    end do
    !
    ! Set midpoint pressures and layer thicknesses
    !
    do k=1,nlev
       do j=1,np
          do i=1,np
             pmid(i,j,k) = hyam(k)*ps0 + hybm(k)*ps(i,j)
             pdel(i,j,k) = pint(i,j,k+1) - pint(i,j,k)
          end do
       end do
    end do

  end subroutine preq_pressure

end module prim_si_mod
