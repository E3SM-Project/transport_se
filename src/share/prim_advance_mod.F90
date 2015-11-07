#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
!#define _DBG_ !DBG

module prim_advance_mod

  use edge_mod,     only: EdgeBuffer_t
  use kinds,        only: real_kind
  use perf_mod,     only: t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only: parallel_t

  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_init

  type (EdgeBuffer_t) :: edge1
  type (EdgeBuffer_t) :: edge2
  type (EdgeBuffer_t) :: edge3p1

  real (kind=real_kind), allocatable :: ur_weights(:)

contains

  !_____________________________________________________________________
  subroutine prim_advance_init(par,integration)

    ! Initialize prim advance

    use edge_mod,       only: initEdgeBuffer
    use dimensions_mod, only: nlev
    use control_mod,    only: qsplit,rsplit

    type(parallel_t), intent(in) :: par
    character(len=*), intent(in) :: integration

    integer :: i

    ! Initialize edge buffer

    if (rsplit==0) then
       call initEdgeBuffer(par,edge3p1,3*nlev+1)
    else
       ! need extra buffer space for dp3d
       call initEdgeBuffer(par,edge3p1,4*nlev+1)
    endif

    ! Compute averaging weights for RK+LF (tstep_type=1) timestepping

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

  !_____________________________________________________________________
  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid, dt, tl, nets, nete, compute_diagnostics)

    ! Apply winds presecribed by test cases

    use control_mod,    only: prescribed_wind, rsplit, qsplit, test_case
    use derivative_mod, only: derivative_t
    use dimensions_mod, only: np, nlev, nlevp, nvar, nc
    use edge_mod,       only: EdgeBuffer_t, initEdgeBuffer
    use element_mod,    only: element_t
    use hybvcoord_mod,  only: hvcoord_t
    use hybrid_mod,     only: hybrid_t
    use time_mod,       only: TimeLevel_t,  timelevel_qdp
    use dcmip_wrapper_mod, only: set_dcmip_1_1_fields, set_dcmip_1_2_fields

    implicit none

    type (element_t),       intent(inout), target   :: elem(:)
    type (derivative_t),    intent(in) :: deriv
    type (hvcoord_t)                   :: hvcoord
    type (hybrid_t),        intent(in) :: hybrid
    real (kind=real_kind),  intent(in) :: dt
    type (TimeLevel_t),     intent(in) :: tl
    integer,                intent(in) :: nets
    integer,                intent(in) :: nete
    logical,                intent(in) :: compute_diagnostics

    real (kind=real_kind) ::  dt2, time, dt_vis, x, eta_ave_w
    real (kind=real_kind) ::  eta_dot_dpdn(np,np,nlevp)
    real (kind=real_kind) ::  dp(np,np)
    real (kind=real_kind) ::  tempdp3d(np,np)
    real (kind=real_kind) ::  tempmass(nc,nc)
    real (kind=real_kind) ::  tempflux(nc,nc,4)
    real (kind=real_kind) ::  deta

    integer :: ie,nm1,n0,np1,nstep,qsplit_stage,k, qn0
    integer :: n,i,j,lx,lenx

    call t_startf('prim_advance_exp')

    ! Extract time indices from time-level stucture

    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep
    time  = tl%nstep*dt

    ! Set weights for mean dynamics fluxes, for timestep type = 1

    qsplit_stage = mod(nstep,qsplit)
    eta_ave_w    = ur_weights(qsplit_stage+1)

    ! Set prescribed fields for DCMIP test cases

    if(test_case(1:8)=="dcmip1-1") then
      call set_dcmip_1_1_fields(elem, hybrid,hvcoord,nets,nete,tl%np1,tl,time)

    else if(test_case(1:8)=="dcmip1-2") then
      call set_dcmip_1_2_fields(elem, hybrid,hvcoord,nets,nete,tl%np1,tl,time)

    else

      ! Apply constant T and v fields for ASP tests

      do ie=nets,nete
        elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,n0)
        elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,n0)
      end do
    endif

    ! Accumulate mean velocity

    do ie=nets,nete
        do k=1,nlev
        elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+&
          eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*elem(ie)%derived%dp(:,:,k)
        elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+&
          eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*elem(ie)%derived%dp(:,:,k)
        enddo
      end do

    call t_stopf('prim_advance_exp')
    end subroutine prim_advance_exp

end module prim_advance_mod

