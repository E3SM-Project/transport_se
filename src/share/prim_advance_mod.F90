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
  use control_mod, only : se_prescribed_wind_2d

  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_init, &
       applyCAMforcing_dynamics, applyCAMforcing, overwrite_SEdensity

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
    use control_mod, only : prescribed_wind, qsplit, tstep_type, rsplit, qsplit, moisture, integration
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

#ifdef TRILINOS
    real (c_double) ,allocatable, dimension(:) :: xstate(:)

!    type (element_t) :: pc_elem(size(elem))
!    type (element_t) :: jac_elem(size(elem))

! state_object is a derived data type passed thru noxinit as a pointer
    type(derived_type) ,target         :: state_object
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    type(derived_type) ,target         :: pre_object
    type(derived_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre
    type(derived_type) ,target         :: jac_object
    type(derived_type) ,pointer         :: jptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_jac

    integer(c_int) :: ierr = 0

  interface

   subroutine noxsolve(vectorSize,vector,v_container,p_container,j_container,ierr) &
!   subroutine noxsolve(vectorSize,vector,v_container) &
     bind(C,name='noxsolve')
    use ,intrinsic :: iso_c_binding
      integer(c_int)                :: vectorSize
      real(c_double)  ,dimension(*) :: vector
      type(c_ptr)                   :: v_container
      type(c_ptr)                   :: p_container  !precon ptr
      type(c_ptr)                   :: j_container  !analytic jacobian ptr
      integer(c_int)                :: ierr         !error flag
    end subroutine noxsolve

  end interface

#endif

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
    if (1==prescribed_wind .and. .not.se_prescribed_wind_2d) then
       time=tl%nstep*dt
       do ie=nets,nete

          ! assume most fields are constant in time
          elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,n0)
          elem(ie)%state%lnps(:,:,np1) = elem(ie)%state%lnps(:,:,n0)
          elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,n0)
          elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,n0)
          elem(ie)%derived%div = 0
	  !get eta_dot_dpdn at n0, ps_v is not used or calculated in this routine
          call asp_advection_vertical(time,hvcoord,elem(ie)%state%ps_v(:,:,n0),&
              eta_dot_dpdn)

          ! accumulate mean fluxes for advection
          if (rsplit==0) then
             elem(ie)%derived%eta_dot_dpdn(:,:,:) = &
                  elem(ie)%derived%eta_dot_dpdn(:,:,:) + eta_dot_dpdn(:,:,:)*eta_ave_w
          else
             ! lagrangian case.  mean vertical velocity = 0
             ! compute dp3d on floating levels
             elem(ie)%derived%eta_dot_dpdn(:,:,:) = 0

             do k=1,nlev
                elem(ie)%state%dp3d(:,:,k,np1) = elem(ie)%state%dp3d(:,:,k,n0)  &
                     + dt*(eta_dot_dpdn(:,:,k+1) - eta_dot_dpdn(:,:,k))
             enddo

          end if
          ! subcycling code uses a mean flux to advect tracers
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(k,dp)
#endif
          do k=1,nlev
             if (rsplit==0) then
                dp(:,:) =&
                     ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                     ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0) 
             else
                dp(:,:) = elem(ie)%state%dp3d(:,:,k,tl%n0)
             end if

             elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+&
                  eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*dp(:,:)
             elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+&
                  eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*dp(:,:)
          enddo
       end do
    endif


    call t_stopf('prim_advance_exp')
    end subroutine prim_advance_exp





  subroutine applyCAMforcing(elem,fvm,hvcoord,np1,np1_qdp,dt_q,nets,nete)
  use dimensions_mod, only : np, nc, nlev, qsize, ntrac
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t
  use control_mod, only : moisture, tracer_grid_type
  use control_mod, only : TRACER_GRIDTYPE_GLL, TRACER_GRIDTYPE_FVM
  use physical_constants, only: Cp
  use fvm_control_volume_mod, only : fvm_struct, n0_fvm
  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  type(fvm_struct)     , intent(inout) :: fvm(:)
  real (kind=real_kind), intent(in) :: dt_q
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: np1,nets,nete,np1_qdp

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  real (kind=real_kind) :: beta(np,np),E0(np,np),ED(np,np),dp0m1(np,np),dpsum(np,np)
  logical :: wet

  wet = (moisture /= "dry")

  do ie=nets,nete
     ! apply forcing to Qdp
     elem(ie)%derived%FQps(:,:,1)=0
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,i,j,v1)
#endif
     !
     ! even when running fvm tracers we need to updates forcing on ps and qv on GLL grid
     !
     ! for fvm tracer qsize is usually 1 (qv)
     !
     do q=1,qsize
        do k=1,nlev
           do j=1,np
              do i=1,np
                 v1 = dt_q*elem(ie)%derived%FQ(i,j,k,q,1)
                 !if (elem(ie)%state%Qdp(i,j,k,q,np1) + v1 < 0 .and. v1<0) then
                 if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) + v1 < 0 .and. v1<0) then
                    !if (elem(ie)%state%Qdp(i,j,k,q,np1) < 0 ) then
                    if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) < 0 ) then
                       v1=0  ! Q already negative, dont make it more so
                    else
                       !v1 = -elem(ie)%state%Qdp(i,j,k,q,np1)
                       v1 = -elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
                    endif
                 endif
                 !elem(ie)%state%Qdp(i,j,k,q,np1) = elem(ie)%state%Qdp(i,j,k,q,np1)+v1
                 elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)+v1
                 if (q==1) then
                    elem(ie)%derived%FQps(i,j,1)=elem(ie)%derived%FQps(i,j,1)+v1/dt_q
                 endif
              enddo
           enddo
        enddo
     enddo
     ! Repeat for the fvm tracers
     do q = 1, ntrac
        do k = 1, nlev
           do j = 1, nc
              do i = 1, nc
                 v1 = fvm(ie)%fc(i,j,k,q)
                 if (fvm(ie)%c(i,j,k,q,n0_fvm) + v1 < 0 .and. v1<0) then
                    if (fvm(ie)%c(i,j,k,q,n0_fvm) < 0 ) then
                       v1 = 0  ! C already negative, dont make it more so
                    else
                       v1 = -fvm(ie)%c(i,j,k,q,n0_fvm)
                    end if
                 end if
                 fvm(ie)%c(i,j,k,q,np1_qdp) = fvm(ie)%c(i,j,k,q,n0_fvm) + v1
                 !                    if (q == 1) then
                 !!XXgoldyXX: Should update the pressure forcing here??!!??
                 !                    elem(ie)%derived%FQps(i,j,1)=elem(ie)%derived%FQps(i,j,1)+v1/dt_q
                 !                  end if
              end do
           end do
        end do
     end do

     if (wet .and. qsize>0) then
        ! to conserve dry mass in the precese of Q1 forcing:
        elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,np1) + &
             dt_q*elem(ie)%derived%FQps(:,:,1)
     endif

#if 0
     ! disabled - energy fixers will be moving into CAM physics
     ! energy fixer for FQps term
     ! dp1 = dp0 + d(FQps)
     ! dp0-dp1 = -d(FQps)
     ! E0-E1 = sum( dp0*ED) - sum( dp1*ED) = sum( dp0-dp1) * ED )
     ! compute E0-E1
     E0=0
     do k=1,nlev
        ED(:,:) = ( 0.5d0* &
             (elem(ie)%state%v(:,:,1,k,np1)**2 + elem(ie)%state%v(:,:,2,k,np1)**2)&
             + cp*elem(ie)%state%T(:,:,k,np1)  &
             + elem(ie)%state%phis(:,:) )

        dp0m1(:,:) = -dt_q*( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%derived%FQps(:,:,1)

        E0(:,:) = E0(:,:) + dp0m1(:,:)*ED(:,:)
     enddo
     ! energy fixer:
     ! Tnew = T + beta
     ! cp*dp*beta  = E0-E1   beta = (E0-E1)/(cp*sum(dp))

     dpsum(:,:) = ( hvcoord%hyai(nlev+1) - hvcoord%hyai(1) )*hvcoord%ps0 + &
          ( hvcoord%hybi(nlev+1) - hvcoord%hybi(1) )*elem(ie)%state%ps_v(:,:,np1)

     beta(:,:)=E0(:,:)/(dpsum(:,:)*cp)
     do k=1,nlev
        elem(ie)%state%T(:,:,k,np1)=elem(ie)%state%T(:,:,k,np1)+beta(:,:)
     enddo
#endif




     ! Qdp(np1) and ps_v(np1) were updated by forcing - update Q(np1)
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,k,i,j,dp)
#endif
     do q=1,qsize
        do k=1,nlev
           do j=1,np
              do i=1,np
                 dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,np1)
                 elem(ie)%state%Q(i,j,k,q) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)/dp
              enddo
           enddo
        enddo
     enddo

     elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
          dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
          dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing



  subroutine applyCAMforcing_dynamics(elem,hvcoord,np1,dt_q,nets,nete)
  use dimensions_mod, only : np, nlev, qsize
  use element_mod, only : element_t
  use hybvcoord_mod, only : hvcoord_t

  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  real (kind=real_kind), intent(in) :: dt_q
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: np1,nets,nete

  ! local
  integer :: i,j,k,ie,q
  real (kind=real_kind) :: v1,dp
  logical :: wet

  do ie=nets,nete
     elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
          dt_q*elem(ie)%derived%FT(:,:,:,1)
     elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
          dt_q*elem(ie)%derived%FM(:,:,:,:,1)
  enddo
  end subroutine applyCAMforcing_dynamics




  subroutine smooth_phis(phis,elem,hybrid,deriv,nets,nete,minf,numcycle)
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : smooth_phis_nudt, hypervis_scaling
  use hybrid_mod, only : hybrid_t
  use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
  use bndry_mod, only : bndry_exchangev
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t , laplace_sphere_wk
  use time_mod, only : TimeLevel_t
  implicit none

  integer :: nets,nete
  real (kind=real_kind), dimension(np,np,nets:nete), intent(inout)   :: phis
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind), intent(in)   :: minf
  integer,               intent(in) :: numcycle

  ! local
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(nets:nete) :: pmin,pmax
  real (kind=real_kind) :: mx,mn
  integer :: nt,ie,ic,i,j,order,order_max, iuse
  logical :: use_var_coef


  ! compute local element neighbor min/max
  do ie=nets,nete
     pstens(:,:,ie)=minval(phis(:,:,ie))
     call edgeVpack(edge3p1,pstens(:,:,ie),1,0,elem(ie)%desc)
  enddo
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     call edgeVunpackMin(edge3p1, pstens(:,:,ie), 1, 0, elem(ie)%desc)
     pmin(ie)=minval(pstens(:,:,ie))
  enddo
  do ie=nets,nete
     pstens(:,:,ie)=maxval(phis(:,:,ie))
     call edgeVpack(edge3p1,pstens(:,:,ie),1,0,elem(ie)%desc)
  enddo
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     call edgeVunpackMax(edge3p1, pstens(:,:,ie), 1, 0, elem(ie)%desc)
     pmax(ie)=maxval(pstens(:,:,ie))
  enddo

  ! order = 1   grad^2 laplacian
  ! order = 2   grad^4 (need to add a negative sign)
  ! order = 3   grad^6
  ! order = 4   grad^8 (need to add a negative sign)
  order_max = 1


  use_var_coef=.true.
  if (hypervis_scaling/=0) then
     ! for tensorHV option, we turn off the tensor except for *last* laplace operator
     use_var_coef=.false.
     if (hypervis_scaling>=3) then
        ! with a 3.2 or 4 scaling, assume hyperviscosity
        order_max = 2
     endif
  endif


  do ic=1,numcycle
     pstens=phis

     do order=1,order_max-1

        do ie=nets,nete
           pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),var_coef=use_var_coef)
           call edgeVpack(edge3p1,pstens(:,:,ie),1,0,elem(ie)%desc)
        enddo
        call bndry_exchangeV(hybrid,edge3p1)
        do ie=nets,nete
           call edgeVunpack(edge3p1, pstens(:,:,ie), 1, 0, elem(ie)%desc)
           pstens(:,:,ie)=pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
        enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
     enddo
     do ie=nets,nete
        pstens(:,:,ie)=laplace_sphere_wk(pstens(:,:,ie),deriv,elem(ie),var_coef=.true.)
     enddo
     if (mod(order_max,2)==0) pstens=-pstens

     do ie=nets,nete
        !  ps(t+1) = ps(t) + Minv * DSS * M * RHS
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  RHS ]
        ! but output of biharminc_wk is of the form M*RHS.  rewrite as:
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  M*RHS/M ]
        ! so we can apply limiter to ps(t) +  (M*RHS)/M
#if 1
        mn=pmin(ie)
        mx=pmax(ie)
        iuse = numcycle+1  ! always apply min/max limiter
#endif
        phis(:,:,ie)=phis(:,:,ie) + &
           smooth_phis_nudt*pstens(:,:,ie)/elem(ie)%spheremp(:,:)


        ! remove new extrema.  could use conservative reconstruction from advection
        ! but no reason to conserve mean PHI.
        if (ic < iuse) then
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)>mx) phis(i,j,ie)=mx
           if (phis(i,j,ie)<mn) phis(i,j,ie)=mn
        enddo
        enddo
        endif


        ! user specified minimum
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)<minf) phis(i,j,ie)=minf
        enddo
        enddo

        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%spheremp(:,:)
        call edgeVpack(edge3p1,phis(:,:,ie),1,0,elem(ie)%desc)
     enddo
     call bndry_exchangeV(hybrid,edge3p1)
     do ie=nets,nete
        call edgeVunpack(edge3p1, phis(:,:,ie), 1, 0, elem(ie)%desc)
        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%rspheremp(:,:)
     enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
  enddo
  end subroutine smooth_phis

  subroutine overwrite_SEdensity(elem, fvm, dt_q, hybrid,nets,nete, np1)
    use fvm_reconstruction_mod, only: reconstruction
    use fvm_filter_mod, only: monotonic_gradient_cart, recons_val_cart
    use dimensions_mod, only : np, nlev, nc,nhe
    use hybrid_mod, only : hybrid_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
    use bndry_mod, only : bndry_exchangev
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t , laplace_sphere_wk
    use time_mod, only : TimeLevel_t
    use fvm_control_volume_mod, only : fvm_struct
    use spelt_mod, only : spelt_struct


    type (element_t) , intent(inout)        :: elem(:)

#if defined(_SPELT)
      type(spelt_struct), intent(inout) :: fvm(:)
#else
      type(fvm_struct), intent(inout) :: fvm(:)
#endif
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    integer, intent(in)                     :: np1
    integer :: ie, k

    real (kind=real_kind)             :: xp,yp, tmpval, dt_q
    integer                           :: i, j,ix, jy, starti,endi,tmpi

    real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons

    if ((nc .ne. 4) .or. (np .ne. 4)) then
      if(hybrid%masterthread) then
        print *,"You are in OVERWRITE SE AIR DENSITY MODE"
        print *,"This only works for nc=4 and np=4"
        print *,"Write a new search algorithm or pay $10000!"
      endif
      stop
    endif
#if defined(_FVM)
    do ie=nets,nete
      call reconstruction(fvm(ie)%psc, fvm(ie),recons)
      call monotonic_gradient_cart(fvm(ie)%psc, fvm(ie),recons, elem(ie)%desc)
      do j=1,np
        do i=1,np
          xp=tan(elem(ie)%cartp(i,j)%x)
          yp=tan(elem(ie)%cartp(i,j)%y)
          ix=i
          jy=j
          ! Search index along "x"  (bisection method)
!           starti = 1
!           endi = nc+1
!           do
!              if  ((endi-starti) <=  1)  exit
!              tmpi = (endi + starti)/2
!              if (xp  >  fvm%acartx(tmpi)) then
!                 starti = tmpi
!              else
!                 endi = tmpi
!              endif
!           enddo
!           ix = starti
!
!         ! Search index along "y"
!           starti = 1
!           endi = nc+1
!           do
!              if  ((endi-starti) <=  1)  exit
!              tmpi = (endi + starti)/2
!              if (yp  >  fvm%acarty(tmpi)) then
!                 starti = tmpi
!              else
!                 endi = tmpi
!              endif
!           enddo
!           jy = starti

          call recons_val_cart(fvm(ie)%psc, xp,yp,fvm(ie)%spherecentroid,recons,ix,jy,tmpval)
          elem(ie)%state%ps_v(i,j,np1)= elem(ie)%state%ps_v(i,j,np1) +&
               dt_q*(tmpval - elem(ie)%state%ps_v(i,j,np1) )/(7*24*60*60)
        end do
      end do
      elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,np1)*elem(ie)%spheremp(:,:)
     call edgeVpack(edge3p1,elem(ie)%state%ps_v(:,:,np1),1,0,elem(ie)%desc)
  enddo
  call bndry_exchangeV(hybrid,edge3p1)
  do ie=nets,nete
     call edgeVunpack(edge3p1, elem(ie)%state%ps_v(:,:,np1), 1, 0, elem(ie)%desc)
     elem(ie)%state%ps_v(:,:,np1)=elem(ie)%state%ps_v(:,:,np1)*elem(ie)%rspheremp(:,:)
  enddo
#endif
  end subroutine overwrite_SEdensity


end module prim_advance_mod

