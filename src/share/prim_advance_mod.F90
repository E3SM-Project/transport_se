#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
!#define _DBG_ !DBG
!
!
module prim_advance_mod
  use edge_mod, only : EdgeBuffer_t
  use kinds, only : real_kind, iulog
  use perf_mod, only: t_startf, t_stopf, t_barrierf, t_adj_detailf ! _EXTERNAL
  use parallel_mod, only : abortmp, parallel_t

  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_init

  type (EdgeBuffer_t) :: edge1
  type (EdgeBuffer_t) :: edge2
  type (EdgeBuffer_t) :: edge3p1

  real (kind=real_kind) :: initialized_for_dt   = 0

  real (kind=real_kind), allocatable :: ur_weights(:)

contains

  subroutine prim_advance_init(par,integration)
    use edge_mod, only : initEdgeBuffer
    use dimensions_mod, only : nlev
    use control_mod, only : qsplit,rsplit
    type(parallel_t) :: par
    character(len=*)    , intent(in) :: integration
    integer :: i

    if (rsplit==0) then
       call initEdgeBuffer(par,edge3p1,3*nlev+1)
    else
       ! need extra buffer space for dp3d
       call initEdgeBuffer(par,edge3p1,4*nlev+1)
    endif

    if(integration == 'semi_imp') then
       call initEdgeBuffer(par,edge1,nlev)
       call initEdgeBuffer(par,edge2,2*nlev)
    end if

    ! compute averaging weights for RK+LF (tstep_type=1) timestepping:
    allocate(ur_weights(qsplit))
    ur_weights(:)=0.0d0

    if(mod(qsplit,2).NE.0)then
       ur_weights(1)=1.0d0/qsplit
       do i=3,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    else
       do i=2,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    endif

  end subroutine prim_advance_init


  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,&
       dt, tl,  nets, nete, compute_diagnostics)
    use bndry_mod, only : bndry_exchangev
    use control_mod, only : prescribed_wind, qsplit, tstep_type, rsplit, qsplit, moisture, integration, test_case
    use derivative_mod, only : derivative_t, vorticity, divergence, gradient, gradient_wk
    use dimensions_mod, only : np, nlev, nlevp, nvar, nc
!    use prim_state_mod, only : prim_printstate
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, initEdgeBuffer
    use element_mod, only : element_t
    use hybvcoord_mod, only : hvcoord_t
    use hybrid_mod, only : hybrid_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use time_mod, only : TimeLevel_t,  timelevel_qdp, tevolve
    use diffusion_mod, only :  prim_diffusion
    use dcmip_wrapper_mod, only: set_dcmip_1_1_fields, set_dcmip_1_2_fields

#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
#endif

#ifndef CAM
    use asp_tests, only : asp_advection_vertical
#else
    use control_mod, only : prescribed_vertwind
#endif

    implicit none

    type (element_t), intent(inout), target   :: elem(:)
    type (derivative_t), intent(in)   :: deriv
    type (hvcoord_t)                  :: hvcoord

    type (hybrid_t)    , intent(in):: hybrid

    real (kind=real_kind), intent(in) :: dt
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete
    logical, intent(in)               :: compute_diagnostics

    ! =================
    ! Local
    ! =================
    real (kind=real_kind) ::  dt2, time, dt_vis, x, eta_ave_w
    real (kind=real_kind) ::  eta_dot_dpdn(np,np,nlevp)
    real (kind=real_kind) ::  dp(np,np)
    real (kind=real_kind) ::  tempdp3d(np,np)
    real (kind=real_kind) ::  tempmass(nc,nc)
    real (kind=real_kind) ::  tempflux(nc,nc,4)
    real (kind=real_kind) ::  deta
    integer :: ie,nm1,n0,np1,nstep,method,qsplit_stage,k, qn0
    integer :: n,i,j,lx,lenx


!JMD    call t_barrierf('sync_prim_advance_exp', hybrid%par%comm)
    call t_startf('prim_advance_exp')
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! timelevel to use for accessing Qdp() to compute virtual temperature
    qn0 = -1    ! -1 = disabled (assume dry dynamics)
    if ( moisture /= "dry") then
       call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel
    endif


! integration = "explicit"
!
!   tstep_type=0  pure leapfrog except for very first timestep   CFL=1
!                    typically requires qsplit=4 or 5
!   tstep_type=1  RK2 followed by qsplit-1 leapfrog steps        CFL=close to qsplit
!                    typically requires qsplit=4 or 5
!   tstep_type=2  RK2-SSP 3 stage (as used by tracers)           CFL=.58
!                    optimal in terms of SSP CFL, but not        CFLSSP=2
!                    optimal in terms of CFL
!                    typically requires qsplit=3
!                    but if windspeed > 340m/s, could use this
!                    with qsplit=1
!   tstep_type=3  classic RK3                                    CFL=1.73 (sqrt(3))
!
!   tstep_type=4  Kinnmark&Gray RK4 4 stage                      CFL=sqrt(8)=2.8
!                 should we replace by standard RK4 (CFL=sqrt(8))?
!                 (K&G 1st order method has CFL=3)
!   tstep_type=5  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
!                 From Paul Ullrich.  3rd order for nonlinear terms also
!                 K&G method is only 3rd order for linear
!                 optimal: for windspeeds ~120m/s,gravity: 340m/2
!                 run with qsplit=1
!                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
!
! integration = "full_imp"
!
!   tstep_type=1  Backward Euler or BDF2 implicit dynamics
!

! default weights for computing mean dynamics fluxes
    eta_ave_w = 1d0/qsplit

    if(tstep_type==0)then
       method=0                ! pure leapfrog
       if (nstep==0) method=1  ! but use RK2 on first step
    else if (tstep_type==1) then
       method=0                           ! LF
       qsplit_stage = mod(nstep,qsplit)
       if (qsplit_stage==0) method=1      ! RK2 on first of qsplit steps
       ! RK2 + LF scheme has tricky weights:
       eta_ave_w=ur_weights(qsplit_stage+1)
    else
       method = tstep_type                ! other RK variants
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! fix dynamical variables, skip dynamics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    if (1==prescribed_wind) then
      time=tl%nstep*dt

      if(test_case(1:8)=="dcmip1-1") then
        call set_dcmip_1_1_fields(elem, hybrid,hvcoord,nets,nete,tl%np1,tl,time)

      else if(test_case(1:8)=="dcmip1-2") then
        call set_dcmip_1_2_fields(elem, hybrid,hvcoord,nets,nete,tl%np1,tl,time)

      else
          ! Apply constant temperature and velocity fields
         do ie=nets,nete
            ! velocity and temperature are held constant in time
            elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,n0)
            elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,n0)
         end do
      endif

      do ie=nets,nete

          ! accumulate mean velocity
          do k=1,nlev
          elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+&
            eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*elem(ie)%derived%dp(:,:,k)
          elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+&
            eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*elem(ie)%derived%dp(:,:,k)
          enddo
        end do

    endif


    call t_stopf('prim_advance_exp')
    end subroutine prim_advance_exp










end module prim_advance_mod

