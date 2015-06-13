#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#if 0
SUBROUTINES:
   prim_advec_tracers_remap_rk2()
      SEM 2D RK2 + monotone remap + hyper viscosity
      SEM 2D RK2 can use sign-preserving or monotone reconstruction

Notes on Lagrange+REMAP advection
dynamics will compute mean fluxes, so that (i.e. for qsplit=3)

    dp(t+3)-dp(t) = -3dt div(Udp_sum/3) - 3dt d(eta_dot_dpdn_sum/3)  + 3dt D(dpdiss_sum/3)

Where the floating lagrangian component:
    dp_star(t+3) = dp(t)  -3dt div(Udp_sum/3)  + 3dt D(dpdiss_sum/3)
OR:
    dp_star(t+3) = dp(t+1) + 3dt d( eta_dot_dpdn_ave(t) )


For RK2 advection of Q:  (example of 2 stage RK for tracers):   dtq = qsplit*dt
For consistency, if Q=1
  dp1  = dp(t)- dtq div[ U1 dp(t)]
  dp2  = dp1  - dtq div[ U2 dp1  ]  + 2*dtq D( dpdiss_ave )
  dp*  = (dp(t) + dp2 )/2
       =  dp(t) - dtq  div[ U1 dp(t) + U2 dp1 ]/2   + dtq D( dpdiss_ave )

so we require:
  U1 = Udp_ave / dp(t)
  U2 = Udp_ave / dp1

For tracer advection:
  Qdp1  = Qdp(t)- dtq div[ U1 Qdp(t)]
  Qdp2  = Qdp1  - dtq div[ U2 Qdp1  ]  + 2*dtq D( Q dpdiss_ave )
  Qdp*  = (Qdp(t) + Qdp2 )/2
       =  Qdp(t) - dtq  div[ U1 Qdp(t) + U2 Qdp1 ]   + dtq D( Q dpdiss_ave )

Qdp1:  limit Q, with Q = Qdp1-before-DSS/(dp1-before-DSS)      with dp1 as computed above
Qdp2:  limit Q, with Q = Qdp2-before-DSS/(dp2-before-DSS)      with dp2 as computed above

For dissipation: Q = Qdp1-after-DSS / dp1-after-DSS


last step:
  remap Qdp* to Qdp(t+1)   [ dp_star(t+1) -> dp(t+1) ]



#endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Begin GPU remap module  !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vertremap_mod

  !**************************************************************************************
  !
  !  Purpose:
  !        Construct sub-grid-scale polynomials using piecewise spline method with
  !        monotone filters.
  !
  !  References: PCM - Zerroukat et al., Q.J.R. Meteorol. Soc., 2005. (ZWS2005QJR)
  !              PSM - Zerroukat et al., Int. J. Numer. Meth. Fluids, 2005. (ZWS2005IJMF)
  !
  !**************************************************************************************

  use kinds, only                  : real_kind,int_kind
  use dimensions_mod, only         : np,nlev,qsize,nlevp,npsq,ntrac,nc
  use hybvcoord_mod, only          : hvcoord_t
  use element_mod, only            : element_t
  use fvm_control_volume_mod, only : fvm_struct
  use spelt_mod, only              : spelt_struct
  use perf_mod, only               : t_startf, t_stopf  ! _EXTERNAL
  use parallel_mod, only           : abortmp, parallel_t
  use control_mod, only : vert_remap_q_alg

  public remap1                  ! remap any field, splines, monotone
  public remap1_nofilter         ! remap any field, splines, no filter
! todo: tweak interface to match remap1 above, rename remap1_ppm:
  public remap_q_ppm             ! remap state%Q, PPM, monotone

  contains

!=======================================================================================================!






!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
!using PPM instead of splines.
subroutine remap_Q_ppm(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  use control_mod, only        : vert_remap_q_alg
  implicit none
  integer,intent(in) :: nx,qsize
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! Local Variables
  integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
  real(kind=real_kind), dimension(       nlev+2 ) :: pio    !Pressure at interfaces for old grid
  real(kind=real_kind), dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
  real(kind=real_kind), dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
  real(kind=real_kind), dimension(       nlev   ) :: z1, z2
  real(kind=real_kind) :: ppmdx(10,0:nlev+1)  !grid spacings
  real(kind=real_kind) :: mymass, massn1, massn2
  integer :: i, j, k, q, kk, kid(nlev)

  call t_startf('remap_Q_ppm')
  do j = 1 , nx
    do i = 1 , nx

      pin(1)=0
      pio(1)=0
      do k=1,nlev
         dpn(k)=dp2(i,j,k)
         dpo(k)=dp1(i,j,k)
         pin(k+1)=pin(k)+dpn(k)
         pio(k+1)=pio(k)+dpo(k)
      enddo



      pio(nlev+2) = pio(nlev+1) + 1.  !This is here to allow an entire block of k threads to run in the remapping phase.
                                      !It makes sure there's an old interface value below the domain that is larger.
      pin(nlev+1) = pio(nlev+1)       !The total mass in a column does not change.
                                      !Therefore, the pressure of that mass cannot either.
      !Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence.
      do k = 1 , gs
        dpo(1   -k) = dpo(       k)
        dpo(nlev+k) = dpo(nlev+1-k)
      enddo

      !Compute remapping intervals once for all tracers. Find the old grid cell index in which the
      !k-th new cell interface resides. Then integrate from the bottom of that old cell to the new
      !interface location. In practice, the grid never deforms past one cell, so the search can be
      !simplified by this. Also, the interval of integration is usually of magnitude close to zero
      !or close to dpo because of minimial deformation.
      !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so
      !I set them equal to each other.
      do k = 1 , nlev
        kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
        !Find the index of the old grid cell in which this new cell's bottom interface resides.
        do while ( pio(kk) <= pin(k+1) )
          kk = kk + 1
        enddo
        kk = kk - 1                   !kk is now the cell index we're integrating over.
        if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds.
                                      !Top bounds match anyway, so doesn't matter what coefficients are used
        kid(k) = kk                   !Save for reuse
        z1(k) = -0.5D0                !This remapping assumes we're starting from the left interface of an old grid cell
                                      !In fact, we're usually integrating very little or almost all of the cell in question
        z2(k) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5 ) / dpo(kk)  !PPM interpolants are normalized to an independent
                                                                        !coordinate domain [-0.5,0.5].
      enddo

      !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the
      !limiting. So anything that depends only on the grid is pre-computed outside the tracer loop.
      ppmdx(:,:) = compute_ppm_grids( dpo )

      !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
      !mass accumulation
      do q = 1 , qsize
        !Accumulate the old mass up to old grid cell interface locations to simplify integration
        !during remapping. Also, divide out the grid spacing so we're working with actual tracer
        !values and can conserve mass. The option for ifndef ZEROHORZ I believe is there to ensure
        !tracer consistency for an initially uniform field. I copied it from the old remap routine.
        masso(1) = 0.
        do k = 1 , nlev
          ao(k) = Qdp(i,j,k,q)
          masso(k+1) = masso(k) + ao(k) !Accumulate the old mass. This will simplify the remapping
          ao(k) = ao(k) / dpo(k)        !Divide out the old grid spacing because we want the tracer mixing ratio, not mass.
        enddo
        !Fill in ghost values. Ignored if vert_remap_q_alg == 2
        do k = 1 , gs
          ao(1   -k) = ao(       k)
          ao(nlev+k) = ao(nlev+1-k)
        enddo
        !Compute monotonic and conservative PPM reconstruction over every cell
        coefs(:,:) = compute_ppm( ao , ppmdx )
        !Compute tracer values on the new grid by integrating from the old cell bottom to the new
        !cell interface to form a new grid mass accumulation. Taking the difference between
        !accumulation at successive interfaces gives the mass inside each cell. Since Qdp is
        !supposed to hold the full mass this needs no normalization.
        massn1 = 0.
        do k = 1 , nlev
          kk = kid(k)
          massn2 = masso(kk) + integrate_parabola( coefs(:,kk) , z1(k) , z2(k) ) * dpo(kk)
          Qdp(i,j,k,q) = massn2 - massn1
          massn1 = massn2
        enddo
      enddo
    enddo
  enddo
  call t_stopf('remap_Q_ppm')
end subroutine remap_Q_ppm


!=======================================================================================================!


!This computes grid-based coefficients from Collela & Woodward 1984.
function compute_ppm_grids( dx )   result(rslt)
  use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=real_kind), intent(in) :: dx(-1:nlev+2)  !grid spacings
  real(kind=real_kind)             :: rslt(10,0:nlev+1)  !grid spacings
  real(kind=real_kind) :: l, c, r, rr
  integer :: j, indB, indE

  call t_startf('compute_ppm_grids')
  !Calculate grid-based coefficients for stage 1 (1--3) and stage 2 (4--10) of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    c = dx(j)
    l = dx(j-1)
    r = dx(j+1)
    rr = dx(j+2)
    rslt( 1,j) = c / ( l + c + r )
    rslt( 2,j) = ( 2.*l + c ) / ( c + r )
    rslt( 3,j) = ( 2.*r + c ) / ( c + l )
    rslt( 4,j) = c / ( c + r )
    rslt( 5,j) = 1. / ( l + c + r + rr )
    rslt( 6,j) = ( 2. * c * r ) / ( c + r )
    rslt( 7,j) = ( c  + l ) / ( 2.*c + r )
    rslt( 8,j) = ( rr + r ) / ( 2.*c + r )
    rslt( 9,j) = c * ( l + c ) / ( 2.*c + r )
    rslt(10,j) = r * ( r + rr ) / ( c + 2.*r )
  enddo
  !Calculate grid-based coefficients for stage 1 of compute_ppm
  rslt( 1,j+1) = dx(j+1) / ( dx(j) + dx(j+1) + dx(j+2) )
  rslt( 2,j+1) = ( 2.*dx(j  ) + dx(j+1) ) / ( dx(j+2) + dx(j+1) )
  rslt( 3,j+1) = ( 2.*dx(j+2) + dx(j+1) ) / ( dx(j  ) + dx(j+1) )
  call t_stopf('compute_ppm_grids')
end function compute_ppm_grids

!=======================================================================================================!



!This computes a limited parabolic interpolant using a net 5-cell stencil, but the stages of computation are broken up into 3 stages
function compute_ppm( a , dx )    result(coefs)
  use control_mod, only: vert_remap_q_alg
  implicit none
  real(kind=real_kind), intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
  real(kind=real_kind), intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
  real(kind=real_kind) ::             coefs(0:2,   nlev  )  !PPM coefficients (for parabola)
  real(kind=real_kind) :: ai (0:nlev  )                     !fourth-order accurate, then limited interface values
  real(kind=real_kind) :: dma(0:nlev+1)                     !An expression from Collela's '84 publication
  real(kind=real_kind) :: da                                !Ditto
  ! Hold expressions based on the grid (which are cumbersome).
  real(kind=real_kind) :: al, ar                            !Left and right interface values for cell-local limiting
  real(kind=real_kind) :: c, rc, cl
  integer :: j, indB, indE

  call t_startf('compute_ppm')
  ! Stage 1: Compute dma for each cell, allowing a 1-cell ghost stencil below and above the domain
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    rc = a(j+1) - a(j)
    cl = a(j  ) - a(j-1)
    da = dx(1,j) * ( dx(2,j) * rc + dx(3,j) * cl )
    dma(j) = minval( (/ abs(da) , 2. * abs( cl ) , 2. * abs( rc ) /) ) * sign(1.D0,da)
    if ( rc * cl <= 0. ) dma(j) = 0.
  enddo

  ! Stage 2: Compute ai for each cell interface in the physical domain (dimension nlev+1)
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    rc = a(j+1) - a(j)
    ai(j) = a(j) + dx(4,j) * rc + dx(5,j) * ( dx(6,j) * ( dx(7,j) - dx(8,j) ) * rc &
         - dx(9,j) * dma(j+1) + dx(10,j) * dma(j) )
  enddo

  ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain
  ! (dimension nlev) using ai on either side and ao within the cell.
  if (vert_remap_q_alg == 2) then
    indB = 3
    indE = nlev-2
  else
    indB = 1
    indE = nlev
  endif
  do j = indB , indE
    al = ai(j-1)
    ar = ai(j  )
    c  = a(j)
    if ( (ar - c) * (c - al) <= 0. ) then
      al = c
      ar = c
    endif
    if ( (ar - al) * (c - (al + ar)/2.) >  (ar - al)**2/6. ) al = 3.*c - 2. * ar
    if ( (ar - al) * (c - (al + ar)/2.) < -(ar - al)**2/6. ) ar = 3.*c - 2. * al
    !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
    coefs(0,j) = 1.5 * c - ( al + ar ) / 4.
    coefs(1,j) = ar - al
    coefs(2,j) = -6. * a(j) + 3. * ( al + ar )
  enddo

  !If we're not using a mirrored boundary condition, then make the two cells bordering the top and bottom
  !material boundaries piecewise constant. Zeroing out the first and second moments, and setting the zeroth
  !moment to the cell mean is sufficient to maintain conservation.
  if (vert_remap_q_alg == 2) then
    coefs(0,1:2) = a(1:2)
    coefs(1:2,1:2) = 0.
    coefs(0,nlev-1:nlev) = a(nlev-1:nlev)
    coefs(1:2,nlev-1:nlev) = 0.D0
  endif
  call t_stopf('compute_ppm')
end function compute_ppm

!=======================================================================================================!


!Simple function computes the definite integral of a parabola in normalized coordinates, xi=(x-x0)/dx,
!given two bounds. Make sure this gets inlined during compilation.
function integrate_parabola( a , x1 , x2 )    result(mass)
  implicit none
  real(kind=real_kind), intent(in) :: a(0:2)  !Coefficients of the parabola
  real(kind=real_kind), intent(in) :: x1      !lower domain bound for integration
  real(kind=real_kind), intent(in) :: x2      !upper domain bound for integration
  real(kind=real_kind)             :: mass
  mass = a(0) * (x2 - x1) + a(1) * (x2*x2 - x1*x1) / 0.2D1 + a(2) * (x2*x2*x2 - x1*x1*x1) / 0.3D1
end function integrate_parabola


!=============================================================================================!



end module vertremap_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End GPU remap module    !!
!! by Rick Archibald, 2010  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!=======================================================================================================!




module prim_advection_mod
!
! two formulations.  both are conservative
! u grad Q formulation:
!
!    d/dt[ Q] +  U grad Q  +  eta_dot dp/dn dQ/dp  = 0
!                            ( eta_dot dQ/dn )
!
!    d/dt[ dp/dn ] = div( dp/dn U ) + d/dn ( eta_dot dp/dn )
!
! total divergence formulation:
!    d/dt[dp/dn Q] +  div( U dp/dn Q ) + d/dn ( eta_dot dp/dn Q ) = 0
!
! for convience, rewrite this as dp Q:  (since dn does not depend on time or the horizonal):
! equation is now:
!    d/dt[dp Q] +  div( U dp Q ) + d( eta_dot_dpdn Q ) = 0
!
!
  use kinds, only              : real_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize, ntrac, nc, nep
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, cp
  use derivative_mod, only     : gradient, vorticity, gradient_wk, derivative_t, divergence, &
                                 gradient_sphere, divergence_sphere
  use element_mod, only        : element_t
  use fvm_control_volume_mod, only        : fvm_struct
  use spelt_mod, only          : spelt_struct
  use filter_mod, only         : filter_t, filter_P
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, smooth, TimeLevel_Qdp
  use prim_si_mod, only        : preq_pressure
  use diffusion_mod, only      : scalar_diffusion, diffusion_init
  use control_mod, only        : integration, test_case, filter_freq_advection,  hypervis_order, &
        statefreq, moisture, TRACERADV_TOTAL_DIVERGENCE, TRACERADV_UGRADQ, &
        nu_q, limiter_option, hypervis_subcycle_q
  use edge_mod, only           : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack, initedgebuffer, &
       edgevunpackmin, ghostbuffer3D_t, initghostbuffer3D
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use viscosity_mod, only      : biharmonic_wk_scalar, biharmonic_wk_scalar_minmax, neighbor_minmax
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only   : abortmp

  implicit none

  private
  save

  public :: Prim_Advec_Init1, Prim_Advec_Init2
  public :: Prim_Advec_Tracers_remap, Prim_Advec_Tracers_remap_rk2
  public :: vertical_remap

  type (EdgeBuffer_t) :: edgeAdv, edgeAdvQ3, edgeAdv_p1, edgeAdvQ2, edgeAdv1,  edgeveloc
  type (ghostBuffer3D_t)   :: ghostbuf_tr

  integer,parameter :: DSSeta = 1
  integer,parameter :: DSSomega = 2
  integer,parameter :: DSSdiv_vdp_ave = 3
  integer,parameter :: DSSno_var = -1

  real(kind=real_kind), allocatable :: qmin(:,:,:), qmax(:,:,:)

  type (derivative_t), public, allocatable   :: deriv(:) ! derivative struct (nthreads)

contains

  subroutine Prim_Advec_Init1(par, n_domains)
    use dimensions_mod, only : nlev, qsize, nelemd
    use parallel_mod, only : parallel_t
    use control_mod, only : use_semi_lagrange_transport
    use interpolate_mod,        only : interpolate_tracers_init
    type(parallel_t) :: par
    integer, intent(in) :: n_domains

    ! Shared buffer pointers.
    ! Using "=> null()" in a subroutine is usually bad, because it makes
    ! the variable have an implicit "save", and therefore shared between
    ! threads. But in this case we want shared pointers.
    real(kind=real_kind), pointer :: buf_ptr(:) => null()
    real(kind=real_kind), pointer :: receive_ptr(:) => null()

    ! this might be called with qsize=0
    ! allocate largest one first
    ! Currently this is never freed. If it was, only this first one should
    ! be freed, as only it knows the true size of the buffer.
    call initEdgeBuffer(par,edgeAdvQ3,max(nlev,qsize*nlev*3), buf_ptr, receive_ptr)  ! Qtens,Qmin, Qmax

    ! remaining edge buffers can share %buf and %receive with edgeAdvQ3
    ! (This is done through the optional 1D pointer arguments.)
    call initEdgeBuffer(par,edgeAdv1,nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(par,edgeAdv,qsize*nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(par,edgeAdv_p1,qsize*nlev + nlev,buf_ptr,receive_ptr)
    call initEdgeBuffer(par,edgeAdvQ2,qsize*nlev*2,buf_ptr,receive_ptr)  ! Qtens,Qmin, Qmax
    call initEdgeBuffer(par,edgeveloc,2*nlev)

    ! Don't actually want these saved, if this is ever called twice.
    nullify(buf_ptr)
    nullify(receive_ptr)

    allocate(deriv(0:n_domains-1))

    ! this static array is shared by all threads, so dimension for all threads (nelemd), not nets:nete:
    allocate (qmin(nlev,qsize,nelemd))
    allocate (qmax(nlev,qsize,nelemd))

    if  (use_semi_lagrange_transport) then
       call initghostbuffer3D(ghostbuf_tr,nlev*qsize,np)
       call interpolate_tracers_init()
    endif

  end subroutine Prim_Advec_Init1

  subroutine Prim_Advec_Init2(hybrid,fvm_corners, fvm_points, spelt_refnep)
    use kinds,          only : longdouble_kind
    use dimensions_mod, only : nc, nep
    use derivative_mod, only : derivinit

    type (hybrid_t), intent(in) :: hybrid
    real(kind=longdouble_kind), intent(in) :: fvm_corners(nc+1)
    real(kind=longdouble_kind), intent(in) :: fvm_points(nc)
    real(kind=longdouble_kind), intent(in) :: spelt_refnep(1:nep)

    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call derivinit(deriv(hybrid%ithr),fvm_corners, fvm_points, spelt_refnep)
  end subroutine Prim_Advec_Init2


!=================================================================================================!

  subroutine Prim_Advec_Tracers_remap( elem , deriv , hvcoord , flt , hybrid , dt , tl , nets , nete )
    use control_mod   , only : use_semi_lagrange_transport
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (filter_t)      , intent(in   ) :: flt
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete


    call Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , flt , hybrid , dt , tl , nets , nete )
  end subroutine Prim_Advec_Tracers_remap




!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! forward-in-time 2 level vertically lagrangian step
!  this code takes a lagrangian step in the horizontal
! (complete with DSS), and then applies a vertical remap
!
! This routine may use dynamics fields at timelevel np1
! In addition, other fields are required, which have to be
! explicitly saved by the dynamics:  (in elem(ie)%derived struct)
!
! Fields required from dynamics: (in
!    omega_p   it will be DSS'd here, for later use by CAM physics
!              we DSS omega here because it can be done for "free"
!    Consistent mass/tracer-mass advection (used if subcycling turned on)
!       dp()   dp at timelevel n0
!       vn0()  mean flux  < U dp > going from n0 to np1
!
! 3 stage
!    Euler step from t     -> t+.5
!    Euler step from t+.5  -> t+1.0
!    Euler step from t+1.0 -> t+1.5
!    u(t) = u(t)/3 + u(t+2)*2/3
!
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
  subroutine Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , flt , hybrid , dt , tl , nets , nete )
    use perf_mod      , only : t_startf, t_stopf            ! _EXTERNAL
    use derivative_mod, only : divergence_sphere
    use control_mod   , only : vert_remap_q_alg, qsplit
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (filter_t)      , intent(in   ) :: flt
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    real (kind=real_kind), dimension(np,np,2     ) :: gradQ
    real (kind=real_kind), dimension(np,np  ,nlev) :: dp_star
    real (kind=real_kind), dimension(np,np  ,nlev) :: dp_np1
    integer :: i,j,k,l,ie,q,nmin
    integer :: nfilt,rkstage,rhs_multiplier
    integer :: n0_qdp, np1_qdp

    call t_barrierf('sync_prim_advec_tracers_remap_k2', hybrid%par%comm)
    call t_startf('prim_advec_tracers_remap_rk2')
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp) !time levels for qdp are not the same
    rkstage = 3 !   3 stage RKSSP scheme, with optimal SSP CFL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RK2 2D advection step
    ! note: stage 3 we take the oppertunity to DSS omega
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! use these for consistent advection (preserve Q=1)
    ! derived%vdp_ave        =  mean horiz. flux:   U*dp
    ! derived%eta_dot_dpdn    =  mean vertical velocity (used for remap)
    ! derived%omega_p         =  advection code will DSS this for the physics, but otherwise
    !                            it is not needed
    ! Also: save a copy of div(U dp) in derived%div(:,:,:,1), which will be DSS'd
    !       and a DSS'ed version stored in derived%div(:,:,:,2)
    do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k, gradQ)
#endif
      do k=1,nlev
        ! div( U dp Q),
        gradQ(:,:,1)=elem(ie)%derived%vn0(:,:,1,k)
        gradQ(:,:,2)=elem(ie)%derived%vn0(:,:,2,k)
        elem(ie)%derived%divdp(:,:,k) = divergence_sphere(gradQ,deriv,elem(ie))
      enddo
      elem(ie)%derived%divdp_proj(:,:,:) = elem(ie)%derived%divdp(:,:,:)
    enddo

    !rhs_multiplier is for obtaining dp_tracers at each stage:
    !dp_tracers(stage) = dp - rhs_multiplier*dt*divdp_proj
    rhs_multiplier = 0
    call euler_step( np1_qdp , n0_qdp  , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSdiv_vdp_ave , rhs_multiplier )

    rhs_multiplier = 1
    call euler_step( np1_qdp , np1_qdp , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSeta         , rhs_multiplier )

    rhs_multiplier = 2
    call euler_step( np1_qdp , np1_qdp , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSomega       , rhs_multiplier )

    !to finish the 2D advection step, we need to average the t and t+2 results to get a second order estimate for t+1.
    call qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nets , nete )

    call t_stopf('prim_advec_tracers_remap_rk2')
  end subroutine prim_advec_tracers_remap_rk2

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nets , nete )
#if USE_CUDA_FORTRAN
    use cuda_mod, only: qdp_time_avg_cuda
#endif
    implicit none
    type(element_t)     , intent(inout) :: elem(:)
    integer             , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete , limiter_option
    integer :: ie
#if USE_CUDA_FORTRAN
    call qdp_time_avg_cuda( elem , rkstage , n0_qdp , np1_qdp , limiter_option , 0d0 , nets , nete )
    return
#endif
    do ie=nets,nete
      elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp) =               &
                   ( elem(ie)%state%Qdp(:,:,:,1:qsize,n0_qdp) + &
                     (rkstage-1)*elem(ie)%state%Qdp(:,:,:,1:qsize,np1_qdp) ) / rkstage
    enddo
  end subroutine qdp_time_avg

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds          , only : real_kind
  use dimensions_mod , only : np, npdg, nlev
  use hybrid_mod     , only : hybrid_t
  use element_mod    , only : element_t
  use derivative_mod , only : derivative_t, divergence_sphere, gradient_sphere, vorticity_sphere
  use edge_mod       , only : edgevpack, edgevunpack
  use bndry_mod      , only : bndry_exchangev
  use hybvcoord_mod  , only : hvcoord_t
#if USE_CUDA_FORTRAN
  use cuda_mod, only: euler_step_cuda
#endif
  implicit none
  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=real_kind), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier

  ! local
  real(kind=real_kind), dimension(np,np                       ) :: divdp, dpdiss
  real(kind=real_kind), dimension(np,np,2                     ) :: gradQ
  real(kind=real_kind), dimension(np,np,2,nlev                ) :: Vstar
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: Qtens
  real(kind=real_kind), dimension(np,np  ,nlev                ) :: dp,dp_star
  real(kind=real_kind), dimension(np,np  ,nlev,qsize,nets:nete) :: Qtens_biharmonic
  real(kind=real_kind), pointer, dimension(:,:,:)               :: DSSvar
  real(kind=real_kind) :: dp0
  real(kind=real_kind) :: tmp1(nlev),tmp2(nlev)
  integer :: ie,q,i,j,k
  integer :: rhs_viss = 0

#if USE_CUDA_FORTRAN
  call euler_step_cuda( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  return
#endif
! call t_barrierf('sync_euler_step', hybrid%par%comm)
!   call t_startf('euler_step')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss = 0

    ! when running lim8, we also need to limit the biharmonic, so that term needs
    ! to be included in each euler step.  three possible algorithms here:
    ! 1) most expensive:
    !     compute biharmonic (which also computes qmin/qmax) during all 3 stages
    !     be sure to set rhs_viss=1
    !     cost:  3 biharmonic steps with 3 DSS
    !
    ! 2) cheapest:
    !     compute biharmonic (which also computes qmin/qmax) only on first stage
    !     be sure to set rhs_viss=3
    !     reuse qmin/qmax for all following stages (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps with 1 DSS
    !     main concern:  viscosity
    !
    ! 3)  compromise:
    !     compute biharmonic (which also computes qmin/qmax) only on last stage
    !     be sure to set rhs_viss=3
    !     compute qmin/qmax directly on first stage
    !     reuse qmin/qmax for 2nd stage stage (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps, 2 DSS
    !
    ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
    do ie = nets , nete
      ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k, q)
#endif
      do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
        dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier*dt*elem(ie)%derived%divdp_proj(:,:,k)
        do q = 1 , qsize
          Qtens_biharmonic(:,:,k,q,ie) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp)/dp(:,:,k)
        enddo
      enddo
    enddo

    ! compute element qmin/qmax
    if ( rhs_multiplier == 0 ) then
      do ie = nets , nete
        do k = 1 , nlev
          do q = 1 , qsize
            qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
            qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
          enddo
        enddo
      enddo
      ! update qmin/qmax based on neighbor data for lim8
      call neighbor_minmax(elem,hybrid,edgeAdvQ2,nets,nete,qmin(:,:,nets:nete),qmax(:,:,nets:nete))
    endif

    ! lets just reuse the old neighbor min/max, but update based on local data
    if ( rhs_multiplier == 1 ) then
      do ie = nets , nete
        do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
          do q = 1 , qsize
            qmin(k,q,ie)=min(qmin(k,q,ie),minval(Qtens_biharmonic(:,:,k,q,ie)))
            qmax(k,q,ie)=max(qmax(k,q,ie),maxval(Qtens_biharmonic(:,:,k,q,ie)))
          enddo
        enddo
      enddo
    endif

    ! get niew min/max values, and also compute biharmonic mixing term
    if ( rhs_multiplier == 2 ) then
      rhs_viss = 3
      ! compute element qmin/qmax
      do ie = nets , nete
        do k = 1  ,nlev
          do q = 1 , qsize
            qmin(k,q,ie)=minval(Qtens_biharmonic(:,:,k,q,ie))
            qmax(k,q,ie)=maxval(Qtens_biharmonic(:,:,k,q,ie))
          enddo
        enddo
      enddo

      call biharmonic_wk_scalar_minmax( elem , qtens_biharmonic , deriv , edgeAdvQ3 , hybrid , &
           nets , nete , qmin(:,:,nets:nete) , qmax(:,:,nets:nete) )
      do ie = nets , nete
#if (defined COLUMN_OPENMP)
        !$omp parallel do private(k, q, dp0)
#endif
        do k = 1 , nlev    !  Loop inversion (AAM)
          dp0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
          do q = 1 , qsize
            ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
            qtens_biharmonic(:,:,k,q,ie) = &
                     -rhs_viss*dt*nu_q*dp0*Qtens_biharmonic(:,:,k,q,ie) / elem(ie)%spheremp(:,:)
          enddo
        enddo
      enddo
    endif



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ie = nets , nete
    ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
    ! all zero so we only have to DSS 1:nlev
    if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
    if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
    if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

    ! Compute velocity used to advance Qdp
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
    do k = 1 , nlev    !  Loop index added (AAM)
      dp(:,:,k) = elem(ie)%derived%dp(:,:,k) - rhs_multiplier * dt * elem(ie)%derived%divdp_proj(:,:,k)
      ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
      ! but that's ok because rhs_multiplier=0 on the first stage:
      Vstar(:,:,1,k) = elem(ie)%derived%vn0(:,:,1,k) / dp(:,:,k)
      Vstar(:,:,2,k) = elem(ie)%derived%vn0(:,:,2,k) / dp(:,:,k)
    enddo

    ! advance Qdp
#if (defined COLUMN_OPENMP)
! commenting this out because of errors at lines 857--860 on Mira
! !$omp parallel do private(q,k,gradQ,dp_star,qtens,dpdiss)
#endif
    do q = 1 , qsize
      do k = 1 , nlev  !  dp_star used as temporary instead of divdp (AAM)
        ! div( U dp Q),
        gradQ(:,:,1) = Vstar(:,:,1,k) * elem(ie)%state%Qdp(:,:,k,q,n0_qdp)
        gradQ(:,:,2) = Vstar(:,:,2,k) * elem(ie)%state%Qdp(:,:,k,q,n0_qdp)
        dp_star(:,:,k) = divergence_sphere( gradQ , deriv , elem(ie) )
        Qtens(:,:,k) = elem(ie)%state%Qdp(:,:,k,q,n0_qdp) - dt * dp_star(:,:,k)
        ! optionally add in hyperviscosity computed above:
        if ( rhs_viss /= 0 ) Qtens(:,:,k) = Qtens(:,:,k) + Qtens_biharmonic(:,:,k,q,ie)
      enddo

      if ( limiter_option == 8) then
        do k = 1 , nlev  ! Loop index added (AAM)
          ! UN-DSS'ed dp at timelevel n0+1:
          dp_star(:,:,k) = dp(:,:,k) - dt * elem(ie)%derived%divdp(:,:,k)
        enddo

        do k=1,nlev
           tmp1(k) = sum(Qtens(:,:,k)*elem(ie)%spheremp(:,:))
        enddo

        ! apply limiter to Q = Qtens / dp_star
        call limiter_optim_iter_full( Qtens(:,:,:) , elem(ie)%spheremp(:,:) , qmin(:,q,ie) , &
                                      qmax(:,q,ie) , dp_star(:,:,:) )

        do k=1,nlev
           tmp2(k) = sum(Qtens(:,:,k)*elem(ie)%spheremp(:,:))
           if (abs(tmp1(k)-tmp2(k)).gt. 1e-5 ) then
              print *,'ie,k=',ie,k,' diff=',abs(tmp1(k)-tmp2(k))
              print *,'sums=',tmp1(k),tmp2(k)
           endif
        enddo

      endif




      ! apply mass matrix, overwrite np1 with solution:
      ! dont do this earlier, since we allow np1_qdp == n0_qdp
      ! and we dont want to overwrite n0_qdp until we are done using it
      do k = 1 , nlev
        elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%spheremp(:,:) * Qtens(:,:,k)
      enddo

    enddo

    if ( DSSopt == DSSno_var ) then
      call edgeVpack(edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
    else
      call edgeVpack(edgeAdv_p1 , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
      ! also DSS extra field
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k)
#endif
      do k = 1 , nlev
        DSSvar(:,:,k) = elem(ie)%spheremp(:,:) * DSSvar(:,:,k)
      enddo
      call edgeVpack( edgeAdv_p1 , DSSvar(:,:,1:nlev) , nlev , nlev*qsize , elem(ie)%desc )
    endif
  enddo

  if ( DSSopt == DSSno_var ) then
    call bndry_exchangeV( hybrid , edgeAdv    )
  else
    call bndry_exchangeV( hybrid , edgeAdv_p1 )
  endif

  do ie = nets , nete
    if ( DSSopt == DSSeta         ) DSSvar => elem(ie)%derived%eta_dot_dpdn(:,:,:)
    if ( DSSopt == DSSomega       ) DSSvar => elem(ie)%derived%omega_p(:,:,:)
    if ( DSSopt == DSSdiv_vdp_ave ) DSSvar => elem(ie)%derived%divdp_proj(:,:,:)

    if ( DSSopt == DSSno_var ) then
      call edgeVunpack( edgeAdv    , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,q)
#endif
      do q = 1 , qsize
        do k = 1 , nlev    !  Potential loop inversion (AAM)
          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
        enddo
      enddo
    else
      call edgeVunpack( edgeAdv_p1 , elem(ie)%state%Qdp(:,:,:,:,np1_qdp) , nlev*qsize , 0 , elem(ie)%desc )
#if (defined COLUMN_OPENMP)
      !$omp parallel do private(q,k)
#endif
      do q = 1 , qsize
        do k = 1 , nlev    !  Potential loop inversion (AAM)
          elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Qdp(:,:,k,q,np1_qdp)
        enddo
      enddo
      call edgeVunpack( edgeAdv_p1 , DSSvar(:,:,1:nlev) , nlev , qsize*nlev , elem(ie)%desc )

      do k = 1 , nlev
        DSSvar(:,:,k) = DSSvar(:,:,k) * elem(ie)%rspheremp(:,:)
      enddo
    endif
  enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
!$OMP BARRIER
#endif
#endif
!   call t_stopf('euler_step')
  end subroutine euler_step


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter_optim_iter_full(ptens,sphweights,minp,maxp,dpmass)
    !THIS IS A NEW VERSION OF LIM8, POTENTIALLY FASTER BECAUSE INCORPORATES KNOWLEDGE FROM
    !PREVIOUS ITERATIONS

    !The idea here is the following: We need to find a grid field which is closest
    !to the initial field (in terms of weighted sum), but satisfies the min/max constraints.
    !So, first we find values which do not satisfy constraints and bring these values
    !to a closest constraint. This way we introduce some mass change (addmass),
    !so, we redistribute addmass in the way that l2 error is smallest.
    !This redistribution might violate constraints thus, we do a few iterations.
    use kinds         , only : real_kind
    use dimensions_mod, only : np, np, nlev
    real (kind=real_kind), dimension(np*np,nlev), intent(inout)            :: ptens
    real (kind=real_kind), dimension(np*np     ), intent(in   )            :: sphweights
    real (kind=real_kind), dimension(      nlev), intent(inout)            :: minp
    real (kind=real_kind), dimension(      nlev), intent(inout)            :: maxp
    real (kind=real_kind), dimension(np*np,nlev), intent(in   ), optional  :: dpmass

    real (kind=real_kind), dimension(np*np,nlev) :: weights
    integer  k1, k, i, j, iter, i1, i2
    integer :: whois_neg(np*np), whois_pos(np*np), neg_counter, pos_counter
    real (kind=real_kind) :: addmass, weightssum, mass
    real (kind=real_kind) :: x(np*np),c(np*np)
    real (kind=real_kind) :: al_neg(np*np), al_pos(np*np), howmuch
    real (kind=real_kind) :: tol_limiter = 1e-15
    integer, parameter :: maxiter = 5

    do k = 1 , nlev
      weights(:,k) = sphweights(:) * dpmass(:,k)
      ptens(:,k) = ptens(:,k) / dpmass(:,k)
    enddo

    do k = 1 , nlev
      c = weights(:,k)
      x = ptens(:,k)

      mass = sum(c*x)

      ! relax constraints to ensure limiter has a solution:
      ! This is only needed if runnign with the SSP CFL>1 or
      ! due to roundoff errors
      if( (mass / sum(c)) < minp(k) ) then
        minp(k) = mass / sum(c)
      endif
      if( (mass / sum(c)) > maxp(k) ) then
        maxp(k) = mass / sum(c)
      endif


!!!      do i2 = 1 , maxIter

      addmass = 0.0d0
      pos_counter = 0
      neg_counter = 0

      ! apply constraints, compute change in mass caused by constraints
      do k1 = 1 , np*np
        if ( ( x(k1) >= maxp(k) ) ) then
          addmass = addmass + ( x(k1) - maxp(k) ) * c(k1)
          x(k1) = maxp(k)
          whois_pos(k1) = -1
        else
          pos_counter = pos_counter+1
          whois_pos(pos_counter) = k1
        endif
        if ( ( x(k1) <= minp(k) ) ) then
          addmass = addmass - ( minp(k) - x(k1) ) * c(k1)
          x(k1) = minp(k)
          whois_neg(k1) = -1
        else
          neg_counter = neg_counter+1
          whois_neg(neg_counter) = k1
        endif
      enddo

      ! iterate to find field that satifies constraints and is l2-norm closest to original
      weightssum = 0.0d0
      if ( addmass > 0 ) then
        do i2 = 1 , maxIter
          weightssum = 0.0
          do k1 = 1 , pos_counter
            i1 = whois_pos(k1)
            weightssum = weightssum + c(i1)
            al_pos(i1) = maxp(k) - x(i1)
          enddo

          if( ( pos_counter > 0 ) .and. ( addmass > tol_limiter * abs(mass) ) ) then
            do k1 = 1 , pos_counter
              i1 = whois_pos(k1)
              howmuch = addmass / weightssum
              if ( howmuch > al_pos(i1) ) then
                howmuch = al_pos(i1)
                whois_pos(k1) = -1
              endif
              addmass = addmass - howmuch * c(i1)
              weightssum = weightssum - c(i1)
              x(i1) = x(i1) + howmuch
            enddo
            !now sort whois_pos and get a new number for pos_counter
            !here neg_counter and whois_neg serve as temp vars
            neg_counter = pos_counter
            whois_neg = whois_pos
            whois_pos = -1
            pos_counter = 0
            do k1 = 1 , neg_counter
              if ( whois_neg(k1) .ne. -1 ) then
                pos_counter = pos_counter+1
                whois_pos(pos_counter) = whois_neg(k1)
              endif
            enddo
          else
            exit
          endif
        enddo
      else
        do i2 = 1 , maxIter
           weightssum = 0.0
           do k1 = 1 , neg_counter
             i1 = whois_neg(k1)
             weightssum = weightssum + c(i1)
             al_neg(i1) = x(i1) - minp(k)
           enddo

           if ( ( neg_counter > 0 ) .and. ( (-addmass) > tol_limiter * abs(mass) ) ) then
             do k1 = 1 , neg_counter
               i1 = whois_neg(k1)
               howmuch = -addmass / weightssum
               if ( howmuch > al_neg(i1) ) then
                 howmuch = al_neg(i1)
                 whois_neg(k1) = -1
               endif
               addmass = addmass + howmuch * c(i1)
               weightssum = weightssum - c(i1)
               x(i1) = x(i1) - howmuch
             enddo
             !now sort whois_pos and get a new number for pos_counter
             !here pos_counter and whois_pos serve as temp vars
             pos_counter = neg_counter
             whois_pos = whois_neg
             whois_neg = -1
             neg_counter = 0
             do k1 = 1 , pos_counter
               if ( whois_pos(k1) .ne. -1 ) then
                 neg_counter = neg_counter+1
                 whois_neg(neg_counter) = whois_pos(k1)
               endif
             enddo
           else
             exit
           endif
         enddo
      endif
!!!      enddo
      ptens(:,k) = x
    enddo

    do k = 1 , nlev
      ptens(:,k) = ptens(:,k) * dpmass(:,k)
    enddo
  end subroutine limiter_optim_iter_full


!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter2d_minmax(Qdp,dp,spheremp,qmin,qmax)
!
! mass conserving limiter (2D only).  to be called just before DSS
!
! in pure 2D advection, the element mass will not be negative before DSS
! this routine will redistribute to remove negative values (conservative)
!
! call with Qdp and assocated dp
!
!
  implicit none
  real (kind=real_kind), intent(inout) :: Qdp(np,np,nlev)
  real (kind=real_kind), intent(in   ) :: spheremp(np,np)
  real (kind=real_kind), intent(in   ) ::  dp(np,np,nlev)
  real (kind=real_kind), intent(in   ) :: qmin(nlev),qmax(nlev)

  ! local
  integer i,j,k
  real (kind=real_kind) :: mass,mass_new,area,mass2
  real (kind=real_kind) :: Q(np,np,nlev)

  Q(:,:,:)=Qdp(:,:,:)/dp(:,:,:)  ! convert to concentration


#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,mass,area,mass2,mass_new,i,j)
#endif
  do k = 1 , nlev
    mass = sum( Qdp(:,:,k)*spheremp(:,:) )
    area = sum( dp(:,:,k)*spheremp(:,:) )

!   if (mass>0) print *,k,mass/area,qmin(k),qmax(k)
    ! max limiter
    if ( maxval(Q(:,:,k)) > qmax(k) ) then
      Q(:,:,k) = qmax(k) - Q(:,:,k)      ! some of these will be negative
      mass2 = area * qmax(k) - mass

      if (mass2 < 0) Q(:,:,k) = -Q(:,:,k)
      mass_new = 0
      do j = 1 , np
        do i = 1 , np
          if ( Q(i,j,k) < 0 ) then
            Q(i,j,k) = 0
          else
            mass_new = mass_new + Q(i,j,k) * dp(i,j,k) * spheremp(i,j)
          endif
        enddo
      enddo

      ! now scale the all positive values to restore mass
      if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass2) / mass_new
      if ( mass2    < 0 ) Q(:,:,k) = -Q(:,:,k)

      Q(:,:,k) = qmax(k) - Q(:,:,k)
    endif

    ! min limiter
    if ( minval(Q(:,:,k)) < qmin(k) ) then
      Q(:,:,k) = Q(:,:,k) - qmin(k)
      mass2 = mass - area * qmin(k)
      ! negative mass.  so reduce all postive values to zero
      ! then increase negative values as much as possible
      if ( mass2 < 0 ) Q(:,:,k) = -Q(:,:,k)
      mass_new = 0
      do j = 1 , np
        do i = 1 , np
          if ( Q(i,j,k) < 0 ) then
            Q(i,j,k) = 0
          else
            mass_new = mass_new + Q(i,j,k) * dp(i,j,k) * spheremp(i,j)
          endif
        enddo
      enddo

      ! now scale the all positive values to restore mass
      if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass2) / mass_new
      if ( mass2    < 0 ) Q(:,:,k) = -Q(:,:,k)

      Q(:,:,k) = Q(:,:,k) + qmin(k)
    endif

    Qdp(:,:,k) = Q(:,:,k) * dp(:,:,k)
  enddo
  end subroutine limiter2d_minmax

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

  subroutine limiter2d_zero(Q)
  ! mass conserving zero limiter (2D only).  to be called just before DSS
  !
  ! this routine is called inside a DSS loop, and so Q had already
  ! been multiplied by the mass matrix.  Thus dont include the mass
  ! matrix when computing the mass = integral of Q over the element
  !
  ! ps is only used when advecting Q instead of Qdp
  ! so ps should be at one timelevel behind Q
  implicit none
  real (kind=real_kind), intent(inout) :: Q(np,np,nlev)

  ! local
  real (kind=real_kind) :: dp(np,np)
  real (kind=real_kind) :: mass,mass_new,ml
  integer i,j,k

  do k = nlev , 1 , -1
    mass = 0
    do j = 1 , np
      do i = 1 , np
        !ml = Q(i,j,k)*dp(i,j)*spheremp(i,j)  ! see above
        ml = Q(i,j,k)
        mass = mass + ml
      enddo
    enddo

    ! negative mass.  so reduce all postive values to zero
    ! then increase negative values as much as possible
    if ( mass < 0 ) Q(:,:,k) = -Q(:,:,k)
    mass_new = 0
    do j = 1 , np
      do i = 1 , np
        if ( Q(i,j,k) < 0 ) then
          Q(i,j,k) = 0
        else
          ml = Q(i,j,k)
          mass_new = mass_new + ml
        endif
      enddo
    enddo

    ! now scale the all positive values to restore mass
    if ( mass_new > 0 ) Q(:,:,k) = Q(:,:,k) * abs(mass) / mass_new
    if ( mass     < 0 ) Q(:,:,k) = -Q(:,:,k)
  enddo
  end subroutine limiter2d_zero

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------




  subroutine vertical_remap(hybrid,elem,fvm,hvcoord,dt,np1,np1_qdp,np1_fvm,nets,nete)
  ! This routine is called at the end of the vertically Lagrangian
  ! dynamics step to compute the vertical flux needed to get back
  ! to reference eta levels
  !
  ! input:
  !     derived%dp()  delta p on levels at beginning of timestep
  !     state%dp3d(np1)  delta p on levels at end of timestep
  ! output:
  !     state%ps_v(np1)          surface pressure at time np1
  !     derived%eta_dot_dpdn()   vertical flux from final Lagrangian
  !                              levels to reference eta levels
  !
  use kinds, only : real_kind
  use hybvcoord_mod, only : hvcoord_t
  use vertremap_mod, only : remap_q_ppm ! _EXTERNAL (actually INTERNAL)
  use control_mod, only : rsplit, tracer_transport_type 
  use parallel_mod, only : abortmp
  use hybrid_mod     , only : hybrid_t
  use derivative_mod, only : interpolate_gll2fvm_points
#if USE_CUDA_FORTRAN
  use cuda_mod, only: vertical_remap_cuda
#endif
  use control_mod, only : se_prescribed_wind_2d

  type (hybrid_t), intent(in) :: hybrid  ! distributed parallel structure (shared)
#if defined(_SPELT)
  type(spelt_struct), intent(inout) :: fvm(:)
  real (kind=real_kind) :: cdp(1:nep,1:nep,nlev,ntrac-1)
  real (kind=real_kind)  :: psc(nep,nep), dpc(nep,nep,nlev),dpc_star(nep,nep,nlev)
#else
  type(fvm_struct), intent(inout) :: fvm(:)
  real (kind=real_kind) :: cdp(1:nc,1:nc,nlev,ntrac)
  real (kind=real_kind)  :: psc(nc,nc), dpc(nc,nc,nlev),dpc_star(nc,nc,nlev)
#endif

  !    type (hybrid_t), intent(in)       :: hybrid  ! distributed parallel structure (shared)
  type (element_t), intent(inout)   :: elem(:)
  type (hvcoord_t)                  :: hvcoord
  real (kind=real_kind)             :: dt

  integer :: ie,i,j,k,np1,nets,nete,np1_qdp,np1_fvm
  integer :: q
  real (kind=real_kind), dimension(np,np,nlev)  :: dp,dp_star
  real (kind=real_kind), dimension(np,np,nlev,2)  :: ttmp

#if USE_CUDA_FORTRAN
  call vertical_remap_cuda(elem,fvm,hvcoord,dt,np1,np1_qdp,nets,nete)
  return
#endif

  call t_startf('vertical_remap')

  ! reference levels:
  !   dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v(i,j)
  !   hybi(1)=0          pure pressure at top of atmosphere
  !   hyai(1)=ptop
  !   hyai(nlev+1) = 0   pure sigma at bottom
  !   hybi(nlev+1) = 1
  !
  ! sum over k=1,nlev
  !  sum(dp(k)) = (hyai(nlev+1)-hyai(1))*ps0 + (hybi(nlev+1)-hybi(1))*ps_v
  !             = -ps0 + ps_v
  !  ps_v =  ps0+sum(dp(k))
  !
  ! reference levels:
  !    dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v
  ! floating levels:
  !    dp_star(k) = dp(k) + dt_q*(eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  ! hence:
  !    (dp_star(k)-dp(k))/dt_q = (eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  !

  do ie=nets,nete

     ! prescribed wind case: density was not updated by dynamics.  
     ! lets use the density from the tracer scheme:
!$omp parallel
!$omp do
     do k = 1 , nlev    
        elem(ie)%state%dp3d(:,:,k,np1)=&
        elem(ie)%derived%dp(:,:,k) - dt*elem(ie)%derived%divdp_proj(:,:,k)
        dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
      enddo
!$omp end do
!$omp single
      if (minval(dp_star)<0) call abortmp('negative layer thickness.  timestep or remap time too large')
      ! compute surface pressure for new levels
      elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
           sum(elem(ie)%state%dp3d(:,:,:,np1),3)
!$omp end single
!$omp do
      do k=1,nlev
         dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                     ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
      enddo
!$omp end do nowait
!$omp end parallel
      
      ! remap the gll tracers from lagrangian levels (dp_star)  to REF levels dp
      call remap_q_ppm(elem(ie)%state%Qdp(:,:,:,:,np1_qdp),np,qsize,dp_star,dp)

  enddo
  call t_stopf('vertical_remap')
  end subroutine vertical_remap




end module prim_advection_mod
