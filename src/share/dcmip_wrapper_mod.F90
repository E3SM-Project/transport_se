!
!  This module contains code for setting initial-conditions
!  and prescribed-winds for some of the DCMIP 2012 test cases.
!
!  https://www.earthsystemcog.org/projects/dcmip-2012/test_cases
!
!  For questions, please contact: david.hall@colorado.edu
!_______________________________________________________________________

#include "config.h"

module dcmip_wrapper_mod

use kinds,          only: rl=>real_kind
use dimensions_mod, only: np, nlev, qsize, nlevp, qsize_d, nelemd
use control_mod,    only: test_case
use hybrid_mod,     only: hybrid_t
use hybvcoord_mod,  only: hvcoord_t
use time_mod,       only: timelevel_t
use element_mod,    only: element_t, elem_state_t, nt=>timelevels, derived_state_t
use dcmip_123_mod,  only: test1_advection_deformation, test1_advection_hadley
use physical_constants, only: p0, g, Rd=>Rgas

implicit none

! parameters for constant-temperature atmosphere cases

real(rl), parameter :: T0 = 300.d0    ! temperature (K)
real(rl), parameter :: H  = Rd*T0/g   ! scale height (m)

! arrays for accumulation of point-wise values

real(rl) :: v_m(np,np,2,nlev)         ! horizontal velocity at layer midpoints
real(rl) :: T_m(np,np,nlev)           ! temperature at layer midpoints
real(rl) :: q_m(np,np,nlev,qsize_d)   ! tracer mixing ratios at midpoints
real(rl) :: p_i(np,np,nlevp)          ! pressure at layer interfaces
real(rl) :: p_m(np,np,nlev)           ! pressure at layer midpoints
real(rl) :: phi_s(np,np)              ! geopotential at the surface
real(rl) :: z_m(np,np,nlev)           ! z at layer midpoints

real(rl) :: omega_m(np,np,nlev)       ! vertical pressure-velocity at midpoints
real(rl) :: eta_dot_dpdn(np,np,nlevp) ! vertical flux at layer interfaces

! single point-wise values

real(rl) :: T,phis,ps,u,v,w,p,z,rho,q(qsize_d),lon,lat
integer  :: ie,i,j,k                  ! loop indices

CONTAINS

!_______________________________________________________________________
subroutine set_dcmip_1_1_fields(elem, hybrid, hvcoord, nets, nete, n0, tl, time)

  ! Wrapper for DCMIP 1-1: 3D Deformational Flow

  type(element_t),    intent(inout), target :: elem(:)
  type(hybrid_t),     intent(in)            :: hybrid
  type(hvcoord_t),    intent(inout)         :: hvcoord
  integer,            intent(in)            :: nets, nete
  integer,            intent(in)            :: n0
  type (timelevel_t), intent(in)            :: tl                       ! time-level structure
  real(rl),           intent(in)            :: time                     ! simulation time in seconds

  integer,  parameter :: zcoords  = 1                                   ! use z coordinate

  do ie= nets, nete

    ! accumulate layer interfaces values
    do k=1,nlev
      z = H * log(1.0d0/hvcoord%etam(k))
      p = p0*hvcoord%etam(k)  ! for evenly spaced eta

      do j=1,np; do i=1,np
        lon=elem(ie)%spherep(i,j)%lon; lat=elem(ie)%spherep(i,j)%lat
        call test1_advection_deformation(time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2),q(3),q(4))
        call cache_midpoint_values(i,j,k,u,v,p,T,q,z)
      enddo; enddo!i,j
    enddo!k

    ! accumulate layer midpoint values
    do k=1,nlevp
      z = H * log(1.0d0/hvcoord%etai(k))
      p = p0*hvcoord%etam(k)  ! for evenly spaced eta

      do j=1,np; do i=1,np
        lon=elem(ie)%spherep(i,j)%lon; lat=elem(ie)%spherep(i,j)%lat
        call test1_advection_deformation(time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2),q(3),q(4))
        call cache_interface_values(i,j,k,p,phis,rho,w)
      enddo; enddo! i,j
    enddo!k

    ! fill element data-structure
    call set_element_state(elem(ie),hvcoord,n0,time)

  enddo!ie

  !call collect_tracer_statistics(elem,hybrid,tl,time,nets,nete)

end subroutine

!_______________________________________________________________________
subroutine set_dcmip_1_2_fields(elem, hybrid, hvcoord, nets, nete, n0, tl, time)

  ! Wrapper for DCMIP 1-2: Hadley-like Flow

  type(element_t),    intent(inout), target :: elem(:)
  type(hybrid_t),     intent(in)            :: hybrid
  type(hvcoord_t),    intent(inout)         :: hvcoord
  integer,            intent(in)            :: nets, nete
  integer,            intent(in)            :: n0
  type (timelevel_t), intent(in)            :: tl                       ! time-level structure
  real(rl),           intent(in)            :: time                     ! simulation time in seconds

  integer,  parameter :: zcoords  = 1                                   ! use z coordinate

  do ie= nets, nete

    ! accumulate layer interfaces values
    do k=1,nlev
      z = H * log(1.0d0/hvcoord%etam(k))
      p = p0*hvcoord%etam(k)  ! for evenly spaced eta

      do j=1,np; do i=1,np
        lon=elem(ie)%spherep(i,j)%lon; lat=elem(ie)%spherep(i,j)%lat
        call test1_advection_hadley (time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2))
        call cache_midpoint_values(i,j,k,u,v,p,T,q,z)
      enddo; enddo!i,j
    enddo!k

    ! accumulate layer midpoint values
    do k=1,nlevp
      z = H * log(1.0d0/hvcoord%etai(k))
      p = p0*hvcoord%etam(k)  ! for evenly spaced eta

      do j=1,np; do i=1,np
        lon=elem(ie)%spherep(i,j)%lon; lat=elem(ie)%spherep(i,j)%lat
        call test1_advection_hadley (time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2))
        call cache_interface_values(i,j,k,p,phis,rho,w)
      enddo; enddo! i,j
    enddo!k

    ! fill element data-structure
    call set_element_state(elem(ie),hvcoord,n0,time)

  enddo!ie

  !call collect_tracer_statistics(elem,hybrid,tl,time,nets,nete)

end subroutine

!_______________________________________________________________________
subroutine set_element_state(e, hvcoord, tl, time)

  ! set element state from accumulated field values

  type (element_t), target, intent(in)  :: e        ! current element
  type(hvcoord_t),          intent(in)  :: hvcoord  ! vertical coordinate
  integer,                  intent(in)  :: tl       ! time level
  real(rl),                 intent(in)  :: time     ! simulation time in seconds

  integer :: qi,k

  type(elem_state_t),    pointer :: s
  type(derived_state_t), pointer :: d

  real(rl) :: dp(np,np,nlev), dn(nlev)
  real(rl) :: dpdn_m(np,np,nlev)

  ! set state variables

  dp = p_i(:,:,2:nlevp) - p_i(:,:,1:nlev)

  s => e%state
  s%v(:,:,:,:,tl)   = v_m
  s%T(:,:,:,tl)     = T_m
  s%dp3d(:,:,:,tl)  = dp
  s%ps_v(:,:,tl)    = p_i(:,:,nlevp)
  s%lnps(:,:,tl)    = log(s%ps_v(:,:,tl))
  s%phis            = phi_s

  ! set derived variables

  d => e%derived
  d%phi         = z_m * g             ! set geopotential height
  d%omega_p     = omega_m/p_m         ! set omega/p
  d%eta_dot_dpdn= eta_dot_dpdn        ! set vertical flux
  d%dp          = dp

  if(time==0.0d0) then

    ! Set additional initial conditions
    s%Q(:,:,:,:)      = q_m
    do qi=1,qsize_d
      s%Qdp(:,:,:,qi,1)  = q_m(:,:,:,qi)*s%dp3d(:,:,:,tl)
      s%Qdp(:,:,:,qi,2)  = q_m(:,:,:,qi)*s%dp3d(:,:,:,tl)
    enddo

  endif

end subroutine

!_______________________________________________________________________
subroutine cache_midpoint_values(i,j,k,u,v,p,T,q,z)
integer   :: i,j,k
real(rl)  :: T,u,v,w,p,z,rho,q(qsize_d)

v_m(i,j,1,k)  = u;  v_m(i,j,2,k)  = v
p_m(i,j,k)    = p;  T_m(i,j,k)    = T
q_m(i,j,k,:)  = q;  z_m(i,j,k)    = z
omega_m(i,j,k)= -g*rho*w

end subroutine

!_______________________________________________________________________
subroutine cache_interface_values(i,j,k,p,phis,rho,w)
integer   :: i,j,k
real(rl)  :: p,phis,rho,w

p_i(i,j,k) = p
phi_s(i,j) = phis
eta_dot_dpdn(i,j,k) = -g*rho*w

end subroutine

!_______________________________________________________________________
subroutine collect_tracer_statistics(elem,hybrid,tl,time,nets,nete)

  ! write global statistics to file to verify test performance

  use control_mod,      only: statefreq
  use global_norms_mod, only: global_integral
  use reduction_mod,    only: parallelmax,parallelmin

  type(element_t),    intent(in), target :: elem(:)
  type(hybrid_t),     intent(in) :: hybrid
  type (timelevel_t), intent(in) :: tl
  real(rl),           intent(in) :: time
  integer,            intent(in) :: nets, nete

  integer :: ie, qi
  integer :: npts = size(phi_s,1)

  real(rl), dimension(np,np,nets:nete) :: local_Qmass
  real(rl), dimension(nets:nete):: local_qmin, local_qmax
  real(rl), dimension(qsize_d)  :: qmass, qmax, qmin

  ! collect statistics at interval specified by statfreq

  if ( mod(tl%nstep,statefreq)/=0) return

  ! measure global mass, max, and min for each tracer

  do qi = 1,qsize
    do ie=nets,nete
      local_qmass(:,:,ie) = elem(ie)%accum%Qmass(:,:,qi,1)
      local_qmin(ie) = minval(elem(ie)%state%q(:,:,:,qi))
      local_qmax(ie) = maxval(elem(ie)%state%q(:,:,:,qi))
    enddo

    qmass(qi) = global_integral(elem, local_qmass(:,:,nets:nete),hybrid,npts,nets,nete)/g
    qmax (qi) = parallelmax(local_qmax,hybrid)
    qmin (qi) = parallelmin(local_qmin,hybrid)

    if (hybrid%masterthread) then
      print *,"tracer ",qi," Qmass = ",qmass(qi)," Qmax=",qmax(qi)," Qmin=",qmin(qi)
    endif

  enddo

end subroutine

!_______________________________________________________________________
subroutine write_level_files()

  ! write hyai and hybi levels to file for dcmip1-1 test

  real(rl), parameter :: z_top = 12000.d0 ! top of atmosphere (m)
  real(rl), parameter :: c     = 2.0d0

  real(rl) :: zi(nlevp), eta_i(nlevp)
  real(rl) :: Am(nlev), Bm(nlev), Ai(nlevp), Bi(nlevp)

  integer :: k

  ! get evenly spaced z
  forall(k=1:nlevp) zi(k) = z_top - z_top*(k-1)/nlev

  ! get eta coord from z for constant temperature atmosphere
  eta_i = exp(-zi/H)

  ! get hybrid coordinates at interfaces
  Bi = ((eta_i - eta_i(1))/(1.0-eta_i(1)))**c
  Ai = eta_i - Bi

  ! get hybrid coordinate at midpoints by averaging interfaces
  Bm = 0.5d0*(Bi(2:nlevp) + Bi(1:nlev))
  Am = 0.5d0*(Ai(2:nlevp) + Ai(1:nlev))

  ! write interface-level file

  open(UNIT=12, FILE="dcmip11i-64.ascii", ACTION="write", STATUS="replace")
  write(12,*) nlevp, " ! hyai";
  do k=1,nlevp; write(12,*) Ai(k); enddo
  write(12,*) nlevp, " ! hybi";
  do k=1,nlevp; write(12,*) Bi(k); enddo
  close(12)

  ! write midpoint-level file

  open(UNIT=12, FILE="dcmip11m-64.ascii", ACTION="write", STATUS="replace")
  write(12,*) nlev, " ! hyam";
  do k=1,nlev; write(12,*) Am(k); enddo
  write(12,*) nlev, " ! hybm";
  do k=1,nlev; write(12,*) Bm(k); enddo
  close(12)

end subroutine

end module
