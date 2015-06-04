#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!#define _DBG_ print *,"file: ",__FILE__," line: ",__LINE__," ithr: ",hybrid%ithr
#define _DBG_
module prim_driver_mod
  use kinds, only : real_kind, iulog, longdouble_kind
  use dimensions_mod, only : np, nlev, nlevp, nelem, nelemd, nelemdmax, GlobalUniqueCols, qsize, nc
  use hybrid_mod, only : hybrid_t
  use quadrature_mod, only : quadrature_t, test_gauss, test_gausslobatto, gausslobatto
  use prim_restart_mod, only : initrestartfile
  use restart_io_mod , only : RestFile,readrestart
  use filter_mod, only : filter_t
  use derivative_mod, only : derivative_t
  use reduction_mod, only : reductionbuffer_ordered_1d_t, red_min, red_max, red_max_int, &
         red_sum, red_sum_int, red_flops, initreductionbuffer

  use element_mod, only : element_t, timelevels,  allocate_element_desc
  use thread_mod, only : omp_get_num_threads
  implicit none
  private
  public :: prim_init1, prim_init2 , prim_run_subcycle, prim_finalize

  type (quadrature_t)   :: gp               ! element GLL points
  type (filter_t)       :: flt              ! Filter struct for v and p grid
  type (filter_t)       :: flt_advection    ! Filter struct for v grid for advection only
  real*8  :: tot_iter
  type (ReductionBuffer_ordered_1d_t), save :: red   ! reduction buffer               (shared)

contains

  subroutine prim_init1(elem, par, dom_mt, Tl)

    ! --------------------------------
    use thread_mod, only : nthreads, omp_get_thread_num, omp_set_num_threads, &
                           vert_num_threads
    ! --------------------------------
    use control_mod, only : runtype, restartfreq, filter_counter, integration, topology, &
         partmethod, while_iter
    ! --------------------------------
    use prim_state_mod, only : prim_printstate_init
    ! --------------------------------
    use namelist_mod, only : readnl
    ! --------------------------------
    use mesh_mod, only : MeshUseMeshFile
    ! --------------------------------
    use time_mod, only : nmax, time_at, timelevel_init, timelevel_t
    ! --------------------------------
    ! --------------------------------
    use mass_matrix_mod, only : mass_matrix
    ! --------------------------------
    use cube_mod,  only : cubeedgecount , cubeelemcount, cubetopology
    ! --------------------------------
    use mesh_mod, only : MeshSetCoordinates, MeshUseMeshFile, MeshCubeTopology, &
         MeshCubeElemCount, MeshCubeEdgeCount
    use cube_mod, only : cube_init_atomic, rotation_init_atomic, set_corner_coordinates, assign_node_numbers_to_elem
    ! --------------------------------
    use metagraph_mod, only : metavertex_t, metaedge_t, localelemcount, initmetagraph
    ! --------------------------------
    use derivative_mod, only : allocate_subcell_integration_matrix
    ! --------------------------------
    use gridgraph_mod, only : gridvertex_t, gridedge_t, allocate_gridvertex_nbrs, deallocate_gridvertex_nbrs
    ! --------------------------------
    use schedtype_mod, only : schedule
    ! --------------------------------
    use schedule_mod, only : genEdgeSched,  PrintSchedule
    ! --------------------------------
    use prim_advection_mod, only: prim_advec_init1
    ! --------------------------------
    use prim_advance_mod, only: prim_advance_init
    ! --------------------------------
    use diffusion_mod, only      : diffusion_init
    ! --------------------------------
    use parallel_mod, only : iam, parallel_t, syncmp, abortmp, global_shared_buf, nrepro_vars
#ifdef _MPI
    use parallel_mod, only : mpiinteger_t, mpireal_t, mpi_max, mpi_sum, haltmp
#endif
    ! --------------------------------
    use metis_mod, only : genmetispart
    ! --------------------------------
    use spacecurve_mod, only : genspacepart
    ! --------------------------------
    use dof_mod, only : global_dof, CreateUniqueIndex, SetElemOffset
    ! --------------------------------
    use params_mod, only : SFCURVE
    ! --------------------------------
    use domain_mod, only : domain1d_t, decompose
    ! --------------------------------
    use physical_constants, only : dd_pi
    ! --------------------------------
    use repro_sum_mod, only: repro_sum, repro_sum_defaultopts, &
         repro_sum_setopts

    implicit none

    type (element_t), pointer :: elem(:)
    type (parallel_t), intent(in) :: par
    type (domain1d_t), pointer :: dom_mt(:)
    type (timelevel_t), intent(out) :: Tl
    ! Local Variables

    type (GridVertex_t), target,allocatable :: GridVertex(:)
    type (GridEdge_t),   target,allocatable :: Gridedge(:)
    type (MetaVertex_t), target,allocatable :: MetaVertex(:)

    integer :: ii,ie, ith
    integer :: nets, nete
    integer :: nelem_edge,nedge
    integer :: nstep
    integer :: nlyr
    integer :: iMv
    integer :: err, ierr, l, j

    real(kind=real_kind), allocatable :: aratio(:,:)
    real(kind=real_kind) :: area(1),xtmp
    character(len=80) rot_type   ! cube edge rotation type

    integer  :: i
    integer,allocatable :: TailPartition(:)
    integer,allocatable :: HeadPartition(:)

    integer total_nelem
    real(kind=real_kind) :: approx_elements_per_task
    integer :: n_domains

    logical :: repro_sum_use_ddpdd, repro_sum_recompute
    real(kind=real_kind) :: repro_sum_rel_diff_max

    ! =====================================
    ! Read in model control information
    ! =====================================

    call readnl(par)
    if (MeshUseMeshFile) then
       total_nelem = MeshCubeElemCount()
    else
       total_nelem = CubeElemCount()
    end if

    approx_elements_per_task = dble(total_nelem)/dble(par%nprocs)
    if  (approx_elements_per_task < 1.0D0) then
       if(par%masterproc) print *,"number of elements=", total_nelem
       if(par%masterproc) print *,"number of procs=", par%nprocs
       call abortmp('There is not enough parallelism in the job, that is, there is less than one elements per task.')
    end if
    ! ====================================
    ! initialize reproducible sum module
    ! ====================================
    call repro_sum_defaultopts(                           &
       repro_sum_use_ddpdd_out=repro_sum_use_ddpdd,       &
       repro_sum_rel_diff_max_out=repro_sum_rel_diff_max, &
       repro_sum_recompute_out=repro_sum_recompute        )
    call repro_sum_setopts(                              &
       repro_sum_use_ddpdd_in=repro_sum_use_ddpdd,       &
       repro_sum_rel_diff_max_in=repro_sum_rel_diff_max, &
       repro_sum_recompute_in=repro_sum_recompute,       &
       repro_sum_master=par%masterproc,                      &
       repro_sum_logunit=6                           )
       if(par%masterproc) print *, "Initialized repro_sum"

    ! ====================================
    ! Set cube edge rotation type for model
    ! unnecessary complication here: all should
    ! be on the same footing. RDL
    ! =====================================
    rot_type="contravariant"

    if (par%masterproc) then
       ! =============================================
       ! Compute total simulated time...
       ! =============================================
       write(iulog,*)""
       write(iulog,*)" total simulated time = ",Time_at(nmax)
       write(iulog,*)""
       ! =============================================
       ! Perform Gauss/Gauss Lobatto tests...
       ! =============================================
       call test_gauss(np)
       call test_gausslobatto(np)
    end if

    ! ===============================================================
    ! Allocate and initialize the graph (array of GridVertex_t types)
    ! ===============================================================

    if (topology=="cube") then

       if (par%masterproc) then
          write(iulog,*)"creating cube topology..."
       end if

       if (MeshUseMeshFile) then
           nelem = MeshCubeElemCount()
           nelem_edge = MeshCubeEdgeCount()
       else
           nelem      = CubeElemCount()
           nelem_edge = CubeEdgeCount()
       end if

       allocate(GridVertex(nelem))
       allocate(GridEdge(nelem_edge))

       do j =1,nelem
          call allocate_gridvertex_nbrs(GridVertex(j))
       end do

       if (MeshUseMeshFile) then
           if (par%masterproc) then
               write(iulog,*) "Set up grid vertex from mesh..."
           end if
           call MeshCubeTopology(GridEdge, GridVertex)
       else
           call CubeTopology(GridEdge,GridVertex)
        end if

       if(par%masterproc)       write(iulog,*)"...done."
    end if
    if(par%masterproc) write(iulog,*)"total number of elements nelem = ",nelem

    !debug  call PrintGridVertex(GridVertex)


    if(partmethod .eq. SFCURVE) then
       if(par%masterproc) write(iulog,*)"partitioning graph using SF Curve..."
       call genspacepart(GridEdge,GridVertex)
    else
        if(par%masterproc) write(iulog,*)"partitioning graph using Metis..."
       call genmetispart(GridEdge,GridVertex)
    endif

    ! ===========================================================
    ! given partition, count number of local element descriptors
    ! ===========================================================
    allocate(MetaVertex(1))
    allocate(Schedule(1))


    nelem_edge=SIZE(GridEdge)

    allocate(TailPartition(nelem_edge))
    allocate(HeadPartition(nelem_edge))
    do i=1,nelem_edge
       TailPartition(i)=GridEdge(i)%tail%processor_number
       HeadPartition(i)=GridEdge(i)%head%processor_number
    enddo

    ! ====================================================
    !  Generate the communication graph
    ! ====================================================
    call initMetaGraph(iam,MetaVertex(1),GridVertex,GridEdge)


    nelemd = LocalElemCount(MetaVertex(1))

    if(nelemd .le. 0) then
       call abortmp('Not yet ready to handle nelemd = 0 yet' )
       stop
    endif
#ifdef _MPI
    call mpi_allreduce(nelemd,nelemdmax,1,MPIinteger_t,MPI_MAX,par%comm,ierr)
#else
    nelemdmax=nelemd
#endif


    if (nelemd>0) then
       allocate(elem(nelemd))
       call allocate_element_desc(elem)
    endif

    ! ====================================================
    !  Generate the communication schedule
    ! ====================================================

    call genEdgeSched(elem,iam,Schedule(1),MetaVertex(1))


    allocate(global_shared_buf(nelemd,nrepro_vars))
    global_shared_buf=0.0_real_kind
    !  nlyr=edge3p1%nlyr
    !  call MessageStats(nlyr)
    !  call testchecksum(par,GridEdge)

    ! ========================================================
    ! load graph information into local element descriptors
    ! ========================================================

    !  do ii=1,nelemd
    !     elem(ii)%vertex = MetaVertex(iam)%members(ii)
    !  enddo

    call syncmp(par)

    ! =================================================================
    ! Set number of domains (for 'decompose') equal to number of threads
    !  for OpenMP across elements, equal to 1 for OpenMP within element
    ! =================================================================
    n_domains = min(Nthreads,nelemd)

    ! =================================================================
    ! Initialize shared boundary_exchange and reduction buffers
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'init shared boundary_exchange buffers'
    call InitReductionBuffer(red,3*nlev,n_domains)
    call InitReductionBuffer(red_sum,5)
    call InitReductionBuffer(red_sum_int,1)
    call InitReductionBuffer(red_max,1)
    call InitReductionBuffer(red_max_int,1)
    call InitReductionBuffer(red_min,1)
    call initReductionBuffer(red_flops,1)


    gp=gausslobatto(np)  ! GLL points

    if (topology=="cube") then
       if(par%masterproc) write(iulog,*) "initializing cube elements..."
       if (MeshUseMeshFile) then
           call MeshSetCoordinates(elem)
       else
           do ie=1,nelemd
               call set_corner_coordinates(elem(ie))
           end do
           call assign_node_numbers_to_elem(elem, GridVertex)
       end if
       do ie=1,nelemd
          call cube_init_atomic(elem(ie),gp%points)
       enddo
    end if

    ! =================================================================
    ! Initialize mass_matrix
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running mass_matrix'
    call mass_matrix(par,elem)
    allocate(aratio(nelemd,1))

    if (topology=="cube") then
       area = 0
       do ie=1,nelemd
          aratio(ie,1) = sum(elem(ie)%mp(:,:)*elem(ie)%metdet(:,:))
       enddo
       call repro_sum(aratio, area, nelemd, nelemd, 1, commid=par%comm)
       area(1) = 4*dd_pi/area(1)  ! ratio correction
       deallocate(aratio)
       if (par%masterproc) &
            write(iulog,'(a,f20.17)') " re-initializing cube elements: area correction=",area(1)

       do ie=1,nelemd
          call cube_init_atomic(elem(ie),gp%points,area(1))
          call rotation_init_atomic(elem(ie),rot_type)
       enddo
    end if


    if(par%masterproc) write(iulog,*) 're-running mass_matrix'
    call mass_matrix(par,elem)


    ! =================================================================
    ! Determine the global degree of freedome for each gridpoint
    ! =================================================================
    if(par%masterproc) write(iulog,*) 'running global_dof'
    call global_dof(par,elem)

    ! =================================================================
    ! Create Unique Indices
    ! =================================================================

    do ie=1,nelemd
       call CreateUniqueIndex(elem(ie)%GlobalId,elem(ie)%gdofP,elem(ie)%idxP)
    enddo

    call SetElemOffset(par,elem, GlobalUniqueCols)

    do ie=1,nelemd
       elem(ie)%idxV=>elem(ie)%idxP
    end do

    !JMD call PrintDofP(elem)
    !JMD call PrintDofV(elem)



    call prim_printstate_init(par)
    ! Initialize output fields for plotting...


    while_iter = 0
    filter_counter = 0

    ! initialize flux terms to 0

    do ie=1,nelemd
       elem(ie)%derived%FM=0.0
       elem(ie)%derived%FQ=0.0
       elem(ie)%derived%FQps=0.0
       elem(ie)%derived%FT=0.0
       elem(ie)%derived%pecnd=0.0

       elem(ie)%accum%Qvar=0
       elem(ie)%accum%Qmass=0
       elem(ie)%accum%Q1mass=0

       elem(ie)%derived%Omega_p=0
       elem(ie)%state%dp3d=0

    enddo


    ! ==========================================================
    !  This routines initalizes a Restart file.  This involves:
    !      I)  Setting up the MPI datastructures
    ! ==========================================================

    if(restartfreq > 0 .or. runtype>=1)  then
       call initRestartFile(elem(1)%state,par,RestFile)
    endif

    !DBG  write(iulog,*) 'prim_init: after call to initRestartFile'

    deallocate(GridEdge)
    do j =1,nelem
       call deallocate_gridvertex_nbrs(GridVertex(j))
    end do
    deallocate(GridVertex)

    do j = 1, MetaVertex(1)%nmembers
       call deallocate_gridvertex_nbrs(MetaVertex(1)%members(j))
    end do
    deallocate(MetaVertex)
    deallocate(TailPartition)
    deallocate(HeadPartition)

    n_domains = min(Nthreads,nelemd)
    call omp_set_num_threads(n_domains)

    ! =====================================
    ! Set number of threads...
    ! =====================================
    if(par%masterproc) then
       write(iulog,*) "Main:NThreads=",NThreads
       write(iulog,*) "Main:n_domains = ",n_domains
    endif

    allocate(dom_mt(0:n_domains-1))
    do ith=0,n_domains-1
       dom_mt(ith)=decompose(1,nelemd,n_domains,ith)
    end do
    ith=0
    nets=1
    nete=nelemd

    call prim_advance_init(par,integration)
    call Prim_Advec_Init1(par, n_domains)
    call diffusion_init(par)

    ! =======================================================
    ! Allocate memory for subcell flux calculations.
    ! =======================================================
    call allocate_subcell_integration_matrix(np, nc)

    call TimeLevel_init(tl)
    if(par%masterproc) write(iulog,*) 'end of prim_init'
  end subroutine prim_init1
!=======================================================================================================!

  subroutine prim_init2(elem, hybrid, nets, nete, tl, hvcoord)

    use parallel_mod, only : parallel_t, haltmp, syncmp, abortmp
    use time_mod, only : timelevel_t, tstep, phys_tscale, timelevel_init, nendstep, nsplit, TimeLevel_Qdp
    use prim_state_mod, only : prim_printstate, prim_diag_scalars
    use filter_mod, only : filter_t, fm_filter_create, taylor_filter_create, &
         fm_transfer, bv_transfer
    use control_mod, only : runtype, integration, filter_mu, filter_mu_advection, test_case, &
        vfile_int, filter_freq, filter_freq_advection, &
         transfer_type, vform, vfile_mid, filter_type, kcut_fm, wght_fm, p_bv, &
         s_bv, topology, moisture, rsplit, qsplit, rk_stage_user,&
         sub_case, &
         limiter_option, nu, nu_q, nu_div, tstep_type, hypervis_subcycle, &
         hypervis_subcycle_q
    use control_mod, only : tracer_transport_type
    use control_mod, only : pertlim                     !used for homme temperature perturbations
    use prim_si_ref_mod, only:  prim_set_mass
#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
#endif
    use thread_mod, only : nthreads
    use derivative_mod, only : derivinit, v2pinit
    use global_norms_mod, only : test_global_integral, print_cfl
    use hybvcoord_mod, only : hvcoord_t
    use prim_advection_mod, only: prim_advec_init2, deriv
    use baroclinic_inst_mod, only : binst_init_state, jw_baroclinic
    use asp_tests, only : asp_tracer, asp_baroclinic, asp_rossby, asp_mountain, asp_gravity_wave, dcmip2_schar
    use dcmip_wrapper_mod, only: set_dcmip_1_1_fields, set_dcmip_1_2_fields

#if USE_CUDA_FORTRAN
    use cuda_mod, only: cuda_mod_init
#endif

    type (element_t), intent(inout) :: elem(:)
    type (hybrid_t), intent(in) :: hybrid

    type (TimeLevel_t), intent(inout)    :: tl              ! time level struct
    type (hvcoord_t), intent(inout)      :: hvcoord         ! hybrid vertical coordinate struct

     integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)


    ! ==================================
    ! Local variables
    ! ==================================

    real (kind=real_kind) :: dt              ! "timestep dependent" timestep
!   variables used to calculate CFL
    real (kind=real_kind) :: dtnu            ! timestep*viscosity parameter
    real (kind=real_kind) :: dt_dyn_vis      ! viscosity timestep used in dynamics
    real (kind=real_kind) :: dt_tracer_vis      ! viscosity timestep used in tracers

    real (kind=real_kind) :: dp


    real (kind=real_kind) :: ps(np,np)       ! surface pressure

    character(len=80)     :: fname
    character(len=8)      :: njusn
    character(len=4)      :: charnum

    real (kind=real_kind) :: Tp(np)     ! transfer function

    integer :: simday
    integer :: i,j,k,ie,iptr,t,q
    integer :: ierr
    integer :: nfrc
    integer :: n0_qdp


    ! ==========================
    ! begin executable code
    ! ==========================
    if (topology == "cube") then
       call test_global_integral(elem, hybrid,nets,nete)
    end if


    ! compute most restrictive dt*nu for use by variable res viscosity:
    if (tstep_type == 0) then
       ! LF case: no tracers, timestep seen by viscosity is 2*tstep
       dt_tracer_vis = 0
       dt_dyn_vis = 2*tstep
       dtnu = 2.0d0*tstep*max(nu,nu_div)
    else
       ! compute timestep seen by viscosity operator:
       dt_dyn_vis = tstep
       if (qsplit>1 .and. tstep_type == 1) then
          ! tstep_type==1: RK2 followed by LF.  internal LF stages apply viscosity at 2*dt
          dt_dyn_vis = 2*tstep
       endif
       dt_tracer_vis=tstep*qsplit

       ! compute most restrictive condition:
       ! note: dtnu ignores subcycling
       dtnu=max(dt_dyn_vis*max(nu,nu_div), dt_tracer_vis*nu_q)
       ! compute actual viscosity timesteps with subcycling
       dt_tracer_vis = dt_tracer_vis/hypervis_subcycle_q
       dt_dyn_vis = dt_dyn_vis/hypervis_subcycle
    endif


    ! ==================================
    ! Initialize derivative structure
    ! ==================================
    call Prim_Advec_Init2(hybrid)

    ! ====================================
    ! In the semi-implicit case:
    ! initialize vertical structure and
    ! related matrices..
    ! ====================================
    ! ==========================================
    ! Initialize pressure and velocity grid
    ! filter matrix...
    ! ==========================================
    if (transfer_type == "bv") then
       Tp    = bv_transfer(p_bv,s_bv,np)
    else if (transfer_type == "fm") then
       Tp    = fm_transfer(kcut_fm,wght_fm,np)
    end if
    if (filter_type == "taylor") then
       flt           = taylor_filter_create(Tp, filter_mu,gp)
       flt_advection = taylor_filter_create(Tp, filter_mu_advection,gp)
    else if (filter_type == "fischer") then
       flt           = fm_filter_create(Tp, filter_mu, gp)
       flt_advection = fm_filter_create(Tp, filter_mu_advection, gp)
    end if



    if (hybrid%masterthread) then
       if (filter_freq>0 .or. filter_freq_advection>0) then
          write(iulog,*) "transfer function type in preq=",transfer_type
          write(iulog,*) "filter type            in preq=",filter_type
          write(*,'(a,99f10.6)') "dynamics: I-mu + mu*Tp(:) = ",&
               (1-filter_mu)+filter_mu*Tp(:)
          write(*,'(a,99f10.6)') "advection: I-mu + mu*Tp(:) = ",&
               (1-filter_mu_advection)+filter_mu_advection*Tp(:)
       endif
    endif

#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
    if (hybrid%ithr==0) then
       call syncmp(hybrid%par)
    end if
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif

    if (topology /= "cube") then
       call abortmp('Error: only cube topology supported for primaitve equations')
    endif

    ! =================================
    ! HOMME stand alone initialization
    ! =================================

    if(runtype >= 1) then
       ! ===========================================================
       ! runtype==1   Exact Restart
       ! runtype==2   Initial run, but take inital condition from Restart file
       ! ===========================================================
       if (hybrid%masterthread) then
          write(iulog,*) 'runtype: RESTART of primitive equations'
       end if

       call ReadRestart(elem,hybrid%ithr,nets,nete,tl)

       ! scale PS to achieve prescribed dry mass
       if (runtype /= 1) &
            call prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)

       if (runtype==2) then
          ! copy prognostic variables:  tl%n0 into tl%nm1
          do ie=nets,nete
             elem(ie)%state%v(:,:,:,:,tl%nm1)=elem(ie)%state%v(:,:,:,:,tl%n0)
             elem(ie)%state%T(:,:,:,tl%nm1)=elem(ie)%state%T(:,:,:,tl%n0) 
             elem(ie)%state%ps_v(:,:,tl%nm1)=elem(ie)%state%ps_v(:,:,tl%n0)
             elem(ie)%state%lnps(:,:,tl%nm1)=elem(ie)%state%lnps(:,:,tl%n0)
          enddo
       endif ! runtype==2
    else  ! initial run  RUNTYPE=0
       ! ===========================================================
       ! Initial Run  - compute initial condition
       ! ===========================================================
       if (hybrid%masterthread) then
          write(iulog,*) ' runtype: INITIAL primitive equations'
       endif
       ! ========================================================
       ! Initialize the test cases
       ! ========================================================

       if(test_case(1:8)=="dcmip1-1") then
          if (hybrid%masterthread) then
            print *,"initializing DCMIP test 1-1: 3d deformational flow"
          endif
          call set_dcmip_1_1_fields(elem, hybrid,hvcoord,nets,nete,tl%n0,tl,time=0.0d0)

       else if(test_case(1:8)=="dcmip1-2") then
          if (hybrid%masterthread) then
            print *,"initializing DCMIP test 1-2: Hadley-like circulation"
          endif
          call set_dcmip_1_2_fields(elem, hybrid,hvcoord,nets,nete,tl%n0,tl,time=0.0d0)

       else if (test_case(1:10) == "baroclinic") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Polvani-Scott-Thomas baroclinic instability test'
          end if

          call binst_init_state(elem, hybrid,nets,nete,hvcoord)
       else if (test_case(1:16) == "asp_gravity_wave") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing ASP gravity wave test'
          end if
          call asp_gravity_wave(elem, hybrid,hvcoord,nets,nete, sub_case)
       else if (test_case(1:12) == "asp_mountain") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing ASP mountain Rossby test'
          end if
          call asp_mountain(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:10) == "asp_rossby") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing ASP Rossby Haurwitz test'
          end if
          call asp_rossby(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:10) == "asp_tracer") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing pure tracer advection tests'
          end if
          call asp_tracer(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:14) == "asp_baroclinic") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Jablonowski and Williamson ASP baroclinic instability test'
          end if
          call asp_baroclinic(elem, hybrid,hvcoord,nets,nete)
       else if (test_case(1:13) == "jw_baroclinic") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing Jablonowski and Williamson baroclinic instability test V1'
          end if
          call jw_baroclinic(elem, hybrid,hvcoord,nets,nete)
       
       else if (test_case(1:12) == "dcmip2_schar") then
          if (hybrid%masterthread) then
             write(iulog,*) 'initializing DCMIP2 test 2-0'
          end if
          call dcmip2_schar(elem, hybrid,hvcoord,nets,nete)
       else
          call abortmp('Error: unrecognized test case')
       endif

       if (hybrid%masterthread) then
          write(iulog,*) '...done'
       end if

       ! scale PS to achieve prescribed dry mass
       call prim_set_mass(elem, tl,hybrid,hvcoord,nets,nete)

       do ie=nets,nete

          elem(ie)%state%T=elem(ie)%state%T &
                * (1.0_real_kind + pertlim)  ! set perlim in ctl_nl namelist for
                                             ! temperature field initial
                                             ! perterbation
       enddo
 
       ! ========================================
       ! Print state and movie output
       ! ========================================
    end if  ! runtype

!$OMP MASTER
    tl%nstep0=2   ! This will be the first full leapfrog step
    if (runtype==1) then
       tl%nstep0=tl%nstep+1            ! restart run: first step = first first full leapfrog step
    endif
    if (runtype==2) then
       ! branch run
       ! reset time counters to zero since timestep may have changed
       nEndStep = nEndStep-tl%nstep ! restart set this to nmax + tl%nstep
       tl%nstep=0
    endif
!$OMP END MASTER
!$OMP BARRIER

    ! For new runs, and branch runs, convert state variable to (Qdp)
    ! because initial conditon reads in Q, not Qdp
    ! restart runs will read dpQ from restart file
    ! need to check what CAM does on a branch run
    if (runtype==0 .or. runtype==2) then
       do ie=nets,nete
          elem(ie)%derived%omega_p(:,:,:) = 0D0
       end do
       do ie=nets,nete
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k, t, q, i, j, dp)
#endif
          do k=1,nlev    !  Loop inversion (AAM)
             do q=1,qsize
                do i=1,np
                   do j=1,np
                      dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                           ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(i,j,tl%n0)
                      
                      elem(ie)%state%Qdp(i,j,k,q,1)=elem(ie)%state%Q(i,j,k,q)*dp
                      elem(ie)%state%Qdp(i,j,k,q,2)=elem(ie)%state%Q(i,j,k,q)*dp
                      
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif




    ! timesteps to use for advective stability:  tstep*qsplit and tstep
    call print_cfl(elem,hybrid,nets,nete,dtnu)

    if (hybrid%masterthread) then
       ! CAM has set tstep based on dtime before calling prim_init2(),
       ! so only now does HOMME learn the timstep.  print them out:
       write(iulog,'(a,2f9.2)') "dt_remap: (0=disabled)   ",tstep*qsplit*rsplit

       if (qsize>0) then
          write(iulog,'(a,2f9.2)') "dt_tracer (SE), per RK stage: ",tstep*qsplit,(tstep*qsplit)/(rk_stage_user-1)
       end if
       write(iulog,'(a,2f9.2)')    "dt_dyn:                  ",tstep
       write(iulog,'(a,2f9.2)')    "dt_dyn (viscosity):      ",dt_dyn_vis
       write(iulog,'(a,2f9.2)')    "dt_tracer (viscosity):   ",dt_tracer_vis

    end if


#if USE_CUDA_FORTRAN
    !Inside this routine, we enforce an OMP BARRIER and an OMP MASTER. It's left out of here because it's ugly
    call cuda_mod_init(elem,hybrid,deriv(hybrid%ithr),hvcoord)
#endif
    if (hybrid%masterthread) write(iulog,*) "initial state:"
    call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)

  end subroutine prim_init2

!=======================================================================================================!




  subroutine prim_run_subcycle(elem, hybrid,nets,nete, dt, tl, hvcoord,nsubstep)
!
!   advance all variables (u,v,T,ps,Q,C) from time t to t + dt_q
!
!   input:
!       tl%nm1   not used
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!
!
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, timelevel_qdp, nsplit
    use control_mod, only: statefreq,&
           energy_fixer, ftype, qsplit, rsplit, test_cfldep, disable_diagnostics
    use prim_state_mod, only : prim_printstate, prim_diag_scalars, prim_energy_halftimes
    use parallel_mod, only : abortmp
    use reduction_mod, only : parallelmax
    use prim_advection_mod, only : vertical_remap
#if USE_CUDA_FORTRAN
    use cuda_mod, only: copy_qdp_h2d, copy_qdp_d2h
#endif


    type (element_t) , intent(inout)        :: elem(:)

    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    integer, intent(in)                     :: nsubstep  ! nsubstep = 1 .. nsplit
    real(kind=real_kind) :: st, st1, dp, dt_q, dt_remap
    integer :: ie, t, q,k,i,j,n
    integer :: n0_qdp,np1_qdp,r, nstep_end

    real (kind=real_kind)                          :: maxcflx, maxcfly
    real (kind=real_kind) :: dp_np1(np,np)
    logical :: compute_diagnostics, compute_energy


    ! ===================================
    ! Main timestepping loop
    ! ===================================
    dt_q = dt*qsplit
    dt_remap = dt_q
    nstep_end = tl%nstep + qsplit



    ! compute diagnostics and energy for STDOUT
    ! compute energy if we are using an energy fixer
    compute_diagnostics=.false.
    compute_energy=energy_fixer > 0

    if (MODULO(nstep_end,statefreq)==0 .or. nstep_end==tl%nstep0 .or. tl%nstep==0) then
       compute_diagnostics=.true.
       compute_energy = .true.
    endif

    if(disable_diagnostics) compute_diagnostics=.false.

    if (compute_diagnostics) &
       call prim_diag_scalars(elem,hvcoord,tl,4,.true.,nets,nete)



    ! E(1) Energy after CAM forcing
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,1,.true.,nets,nete)

    ! qmass and variance, using Q(n0),Qdp(n0)
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,1,.true.,nets,nete)


#if USE_CUDA_FORTRAN
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call copy_qdp_h2d( elem , n0_qdp )
#endif

    ! take a timestep of trancers and dynamis
    call prim_step(elem, hybrid,nets,nete, dt, tl, hvcoord,compute_diagnostics,1)
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)

    ! vertical remap
    call vertical_remap(hybrid,elem,hvcoord,dt_remap,tl%np1,np1_qdp,nets,nete)

#if USE_CUDA_FORTRAN
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)
    call copy_qdp_d2h( elem , np1_qdp )
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! time step is complete.  update some diagnostic variables:
    ! lnps (we should get rid of this)
    ! Q    (mixing ratio)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       elem(ie)%state%lnps(:,:,tl%np1)= LOG(elem(ie)%state%ps_v(:,:,tl%np1))
#if (defined COLUMN_OPENMP)
       !$omp parallel do default(shared), private(k,q,dp_np1)
#endif
       do k=1,nlev    !  Loop inversion (AAM)
          !if (k == 1) then
           !write(*,*) "In prim run there are ", omp_get_num_threads(), " in the current team in parallel region"
          !endif
          dp_np1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
               ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%np1)
          do q=1,qsize
             elem(ie)%state%Q(:,:,k,q)=elem(ie)%state%Qdp(:,:,k,q,np1_qdp)/dp_np1(:,:)
          enddo
       enddo
    enddo





    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - 2*dt
    !   u(n0)    dynamics at  t+dt_remap - dt
    !   u(np1)   dynamics at  t+dt_remap
    !
    !   Q(1)   Q at t+dt_remap
    if (compute_diagnostics) call prim_diag_scalars(elem,hvcoord,tl,2,.false.,nets,nete)
    if (compute_energy) call prim_energy_halftimes(elem,hvcoord,tl,2,.false.,nets,nete)

    if (compute_diagnostics) then
       call prim_diag_scalars(elem,hvcoord,tl,3,.false.,nets,nete)
       call prim_energy_halftimes(elem,hvcoord,tl,3,.false.,nets,nete)
     endif

    ! =================================
    ! update dynamics time level pointers
    ! =================================
    call TimeLevel_update(tl,"leapfrog")

    ! now we have:
    !   u(nm1)   dynamics at  t+dt_remap - dt       (Robert-filtered)
    !   u(n0)    dynamics at  t+dt_remap
    !   u(np1)   undefined


    ! ============================================================
    ! Print some diagnostic information
    ! ============================================================
    if (compute_diagnostics) then
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if
  end subroutine prim_run_subcycle






  subroutine prim_step(elem, hybrid,nets,nete, dt, tl, hvcoord, compute_diagnostics,rstep)
!
!   Take qsplit dynamics steps and one tracer step
!   for vertically lagrangian option, this subroutine does only the horizontal step
!
!   input:
!       tl%nm1   not used
!       tl%n0    data at time t
!       tl%np1   new values at t+dt_q
!
!   then we update timelevel pointers:
!       tl%nm1 = tl%n0
!       tl%n0  = tl%np1
!   so that:
!       tl%nm1   tracers:  t    dynamics:  t+(qsplit-1)*dt
!       tl%n0    time t + dt_q
!
!
    use hybvcoord_mod, only : hvcoord_t
    use time_mod, only : TimeLevel_t, timelevel_update, nsplit
    use control_mod, only: statefreq, integration, ftype, qsplit, nu_p, test_cfldep, rsplit
    use control_mod, only : tracer_transport_type
    use control_mod, only : tracer_grid_type, TRACER_GRIDTYPE_GLL
    use prim_advance_mod, only : prim_advance_exp
    use prim_advection_mod, only : prim_advec_tracers_remap, deriv
    use derivative_mod, only : subcell_integration
    use parallel_mod, only : abortmp
    use reduction_mod, only : parallelmax
    use time_mod,    only : time_at

    type (element_t) , intent(inout)        :: elem(:)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)
    type (hvcoord_t), intent(in)      :: hvcoord         ! hybrid vertical coordinate struct

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    real(kind=real_kind), intent(in)        :: dt  ! "timestep dependent" timestep
    type (TimeLevel_t), intent(inout)       :: tl
    integer, intent(in)                     :: rstep ! vertical remap subcycling step

    real(kind=real_kind) :: st, st1, dp, dt_q
    integer :: ie, t, q,k,i,j,n

    real (kind=real_kind)                          :: maxcflx, maxcfly

    real (kind=real_kind) ::  tempdp3d(np,np), x
    real (kind=real_kind) ::  tempmass(nc,nc)
    real (kind=real_kind) ::  tempflux(nc,nc,4)

    real (kind=real_kind) :: dp_np1(np,np)
    logical :: compute_diagnostics

    dt_q = dt*qsplit
 
    ! ===============
    ! initialize mean flux accumulation variables and save some variables at n0
    ! for use by advection
    ! ===============
    do ie=nets,nete
      elem(ie)%derived%eta_dot_dpdn=0     ! mean vertical mass flux
      elem(ie)%derived%vn0=0              ! mean horizontal mass flux
      elem(ie)%derived%omega_p=0
      if (nu_p>0) then
         elem(ie)%derived%dpdiss_ave=0
         elem(ie)%derived%dpdiss_biharmonic=0
      endif

        ! save dp at time t for use in tracers
#if (defined COLUMN_OPENMP)
!$omp parallel do default(shared), private(k)
#endif
         do k=1,nlev
            elem(ie)%derived%dp(:,:,k)=&
                 ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                 ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,tl%n0)
         enddo
    enddo

    ! ===============
    ! compute analytic winds
    ! ===============
    call prim_advance_exp(elem, deriv(hybrid%ithr), hvcoord,   &
         hybrid, dt, tl, nets, nete, compute_diagnostics)

    ! ===============
    ! advance tracers
    ! ===============
    call Prim_Advec_Tracers_remap(elem, deriv(hybrid%ithr),hvcoord,flt_advection,hybrid,&
         dt_q,tl,nets,nete)


  end subroutine prim_step


!=======================================================================================================!


  subroutine prim_finalize(hybrid)
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)


    ! ==========================
    ! end of the hybrid program
    ! ==========================
  end subroutine prim_finalize



end module prim_driver_mod



