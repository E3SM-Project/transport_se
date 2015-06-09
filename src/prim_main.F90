#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

program prim_main

  use hybvcoord_mod,    only: hvcoord_t, hvcoord_init
  use parallel_mod,     only: parallel_t, initmp, haltmp
  use hybrid_mod,       only: hybrid_t
  use time_mod,         only: tstep, nendstep, timelevel_t, TimeLevel_init
  use dimensions_mod,   only: nelemd
  use domain_mod,       only: domain1d_t
  use element_mod,      only: element_t
  use common_io_mod,    only: output_dir
  use restart_io_mod,   only: writerestart
  use hybrid_mod,       only: hybrid_create
  use common_movie_mod, only: nextoutputstep
  use perf_mod,         only: t_initf, t_prf, t_finalizef, t_startf, t_stopf
  use control_mod,      only: restartfreq, vfile_mid, vfile_int, runtype
  use prim_driver_mod,  only: prim_init1, prim_init2, prim_finalize, prim_run_subcycle
  use thread_mod,       only: nthreads, vert_num_threads, omp_get_thread_num, &
                              omp_set_num_threads, omp_get_nested, &
                              omp_get_num_threads, omp_get_max_threads
#ifdef PIO_INTERP
  use interp_movie_mod,       only : interp_movie_output, interp_movie_finish, interp_movie_init
  use interpolate_driver_mod, only : interpolate_driver
#else
  use prim_movie_mod,   only : prim_movie_output, prim_movie_finish,prim_movie_init
#endif

  implicit none

  type (element_t),  pointer :: elem(:)   ! array of elements
  type (domain1d_t), pointer :: dom_mt(:) ! element start and end indices
  type (parallel_t)   :: par              ! parallel structure for MPI
  type (hybrid_t)     :: hybrid           ! parallel structure OpenMP+MPI
  type (TimeLevel_t)  :: tl               ! time level struct
  type (hvcoord_t)    :: hvcoord          ! hybrid vertical coordinates

  real*8  :: timeit, et, st
  integer :: nets,nete
  integer :: ithr
  integer :: ierr
  integer :: nstep

  ! Initialize MPI

  par=initmp()

  ! Initialize performance timers

  call t_initf('input.nl',LogPrint=par%masterproc, Mpicom=par%comm, MasterTask=par%masterproc)
  call t_startf('Total')

  ! Read namelist, generate mesh, allocate memory, initialize mass-matrix (prim_init1)

  call t_startf('prim_init1')
  call prim_init1(elem, par,dom_mt,tl)
  call t_stopf('prim_init1')

  ! Verify nested OpenMP, if needed

#if (defined HORIZ_OPENMP && defined COLUMN_OPENMP)
   call omp_set_nested(.true.)
   if (omp_get_nested() == 0) then
     call haltmp("Nested threading required but not available. Set OMP_NESTED=true")
   endif
#endif

  ! Print process and thread decomposition

#if (defined HORIZ_OPENMP)
  !$OMP PARALLEL NUM_THREADS(nthreads), DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid)
  call omp_set_num_threads(vert_num_threads)
#endif
  ithr = omp_get_thread_num()
  nets = dom_mt(ithr)%start
  nete = dom_mt(ithr)%end
#if (defined HORIZ_OPENMP)
  !$OMP CRITICAL
#endif
  if (par%rank<100) then 
     write(6,9) par%rank,ithr,nets,nete 
  endif
9 format("process: ",i2,1x,"thread: ",i2,1x,"element limits: ",i5," - ",i5)
#if (defined HORIZ_OPENMP)
  !$OMP END CRITICAL
  !$OMP END PARALLEL
#endif

  ! Read vertical coordinates from file

  ithr   = omp_get_thread_num()
  hybrid = hybrid_create(par,ithr,1)
  nets   = 1
  nete   = nelemd
  hvcoord = hvcoord_init(vfile_mid, vfile_int, .true., hybrid%masterthread, ierr)
  if (ierr /= 0) then
     call haltmp("error in hvcoord_init")
  end if

  ! Initialize: derivatives, restart runs, test_cases, and filters (prim_init2)

  if(par%masterproc) print *,"Primitive Equation Initialization..."
#if (defined HORIZ_OPENMP)
  !$OMP PARALLEL NUM_THREADS(nthreads), DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid)
  call omp_set_num_threads(vert_num_threads)
#endif
  ithr    = omp_get_thread_num()
  hybrid  = hybrid_create(par,ithr,nthreads)
  nets    = dom_mt(ithr)%start
  nete    = dom_mt(ithr)%end
  call t_startf('prim_init2')
  call prim_init2(elem,  hybrid,nets,nete,tl, hvcoord)
  call t_stopf('prim_init2')
#if (defined HORIZ_OPENMP)
  !$OMP END PARALLEL
#endif

  ! Ensure output directory exists

  ithr    = omp_get_thread_num()
  hybrid  = hybrid_create(par,ithr,1)
  if (hybrid%masterthread) then 
     open(unit=447,file=trim(output_dir) // "/output_dir_test",iostat=ierr)
     if ( ierr==0 ) then
        print *,'Directory ',output_dir, ' does exist: initialing IO'
        close(447)
     else
        print *,'Error creating file in directory ',output_dir
        call haltmp("Please be sure the directory exist or specify 'output_dir' in the namelist.")
     end if
  endif
  
  ! Initialize history files

#ifdef PIO_INTERP
  call interp_movie_init( elem, hybrid, 1, nelemd, hvcoord, tl )
#else
  call prim_movie_init( elem, par, hvcoord, tl )
#endif

  ! Write initial state for new runs

  if (runtype == 0 ) then
#ifdef PIO_INTERP
     call interp_movie_output(elem, tl, hybrid, 0d0, 1, nelemd, hvcoord=hvcoord)
#else
     call prim_movie_output(elem, tl, hvcoord, hybrid, 1,nelemd)
#endif
  endif

  ! Perform main timestepping Loop

  if(par%masterproc) print *,"Entering main timestepping loop"
  do while(tl%nstep < nEndStep)
#if (defined HORIZ_OPENMP)
     !$OMP PARALLEL NUM_THREADS(nthreads), DEFAULT(SHARED), PRIVATE(ithr,nets,nete,hybrid)
     call omp_set_num_threads(vert_num_threads)
#endif
     ithr   = omp_get_thread_num()
     hybrid = hybrid_create(par,ithr,nthreads)
     nets   = dom_mt(ithr)%start
     nete   = dom_mt(ithr)%end
     nstep  = nextoutputstep(tl)

     ! Integrate PDEs until next output-step

     do while(tl%nstep<nstep)
        call t_startf('prim_run')
        call prim_run_subcycle(elem, hybrid,nets,nete, tstep, tl, hvcoord,1)
        call t_stopf('prim_run')
     end do
#if (defined HORIZ_OPENMP)
     !$OMP END PARALLEL
#endif

     ! Write history file

     ithr   = omp_get_thread_num()
     hybrid = hybrid_create(par,ithr,1)
#ifdef PIO_INTERP
     call interp_movie_output(elem, tl, hybrid, 0d0, 1, nelemd, hvcoord=hvcoord)
#else
     call prim_movie_output(elem, tl, hvcoord, hybrid, 1,nelemd)
#endif

     ! Write restart files if needed

      if((restartfreq > 0) .and. (MODULO(tl%nstep,restartfreq) ==0)) then
        call WriteRestart(elem, ithr,1,nelemd,tl)
     endif
  end do

  ! Write final history file

  if(par%masterproc) print *,"Finished main timestepping loop",tl%nstep
  call prim_finalize(hybrid)
  if(par%masterproc) print *,"closing history files"
#ifdef PIO_INTERP
  call interp_movie_finish
#else
  call prim_movie_finish
#endif

  ! Stop performance timers and write timing data

  call t_stopf('Total')
  if(par%masterproc) print *,"writing timing data"
  call t_prf('HommeTime', par%comm)
  if(par%masterproc) print *,"calling t_finalizef"
  call t_finalizef()

  ! Halt MPI and exit

  call haltmp("exiting program...")
end program prim_main








