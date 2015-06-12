#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module namelist_mod
  !-----------------
  use kinds, only : real_kind, iulog
  !-----------------
  use params_mod, only : recursive, sfcurve
  !-----------------
  use cube_mod, only : rotate_grid
  use control_mod, only : &
       MAX_STRING_LEN,&
       MAX_FILE_LEN,  &
       partmethod,    &       ! Mesh partitioning method (METIS)
       topology,      &       ! Mesh topology
       test_case,     &       ! test case
       uselapi,       &
       multilevel,    &
       numnodes,      &
       sub_case,      &
       tasknum,       &       ! used dg model in AIX machine
       remapfreq,     &       ! number of steps per remapping call
       remap_type,    &       ! selected remapping option
       statefreq,     &       ! number of steps per printstate call
       accumfreq,     &       ! frequency in steps of accumulation
       accumstart,    &       ! model day to start accumulating state variables
       accumstop,     &       ! model day to stop  accumulating state variables
       restartfreq,   &
       restartfile,   &       ! name of the restart file for INPUT
       restartdir,    &       ! name of the restart directory for OUTPUT
       runtype,       &
       integration,   &       ! integration method
       tracer_advection_formulation, &   ! conservation or non-conservation formulaton
       tstep_type, &
       cubed_sphere_map, &
       compute_mean_flux, &
       qsplit, &
       rsplit, &
       physics, &
       rk_stage_user, &
       LFTfreq,       &
       TRACERADV_TOTAL_DIVERGENCE, &
       TRACERADV_UGRADQ, &
       prescribed_wind, &
       ftype,        &
       energy_fixer,        &
       limiter_option, &
       fine_ne,       &
       max_hypervis_courant, &
       nu,            &
       nu_s,          &
       nu_q,          &
       nu_div,          &
       nu_p,          &
       nu_top,        &
       hypervis_scaling,   & ! use tensor HV instead of scalar coefficient
       disable_diagnostics, & ! Use to disable diagnostics for timing reasons
       psurf_vis,    &
       hypervis_order,    &
       hypervis_power,    &
       hypervis_subcycle, &
       hypervis_subcycle_q, &
       smooth_phis_numcycle, &
       smooth_sgh_numcycle, &
       smooth_phis_nudt, &
       initial_total_mass, &  ! set > 0 to set the initial_total_mass
       u_perturb,     &          ! J&W bareclinic test perturbation size
       moisture,      &
       filter_type,   &
       transfer_type, &
       filter_freq,   &
       filter_mu,     &
       filter_freq_advection,   &
       filter_mu_advection,     &
       p_bv,          &
       s_bv,          &
       wght_fm,       &
       kcut_fm,       &
       vform,           &
       vfile_mid,       &
       vfile_int,       &
       vert_remap_q_alg, &
       pertlim,      &
       tracer_transport_type,           &
       TRACERTRANSPORT_SE_GLL,          &
       tracer_grid_type,                &
       TRACER_GRIDTYPE_GLL,             &
       test_cfldep, &
       se_prescribed_wind_2d

  !-----------------
  use thread_mod, only : nthreads, nthreads_accel, omp_get_max_threads, vert_num_threads
  !-----------------
  use dimensions_mod, only : ne, np, npdg, nnodes, nmpi_per_node, npart, qsize, qsize_d, set_mesh_dimensions
  !-----------------
  use time_mod, only : tstep, ndays,nmax, nendstep,secpday, smooth, secphr, nsplit, phys_tscale
  !-----------------
  use parallel_mod, only : parallel_t,  iam, abortmp, &
       partitionfornodes, useframes, mpireal_t, mpilogical_t, mpiinteger_t, mpichar_t
  !-----------------

  use common_io_mod, only : &
       output_prefix,       &
       output_type,         &
       output_timeunits,    &
       output_start_time,   &
       output_end_time,     &
       output_frequency,    &
       output_dir,          &
       output_varnames1,    &
       output_varnames2,    &
       output_varnames3,    &
       output_varnames4,    &
       output_varnames5,    &
       max_output_variables,&
       max_output_streams,  &
       num_io_procs,       &
       io_stride,          &
       varname_len,         &
       infilenames,         &
       MAX_INFILES

  use physical_constants, only: omega

  use common_movie_mod, only : setvarnames

  implicit none
  private
!
! This module should contain no global data and should only be 'use'd where readnl is called
!
  public :: readnl

 contains
! ============================================
! readnl:
!
!  Read in the namelists...
!
! ============================================

  subroutine readnl(par)
    use mesh_mod, only : MeshOpen
    type (parallel_t), intent(in) ::  par
    character(len=MAX_FILE_LEN) :: mesh_file
    integer :: se_ftype, se_limiter_option
    integer :: se_phys_tscale, se_nsplit
    integer :: i, ii, j
    integer  :: ierr
    character(len=80) :: errstr, arg
    real(kind=real_kind) :: dt_max
    character(len=32) :: tracer_transport_method = 'se_gll'

    ! ============================================
    ! Namelists
    ! ============================================

    namelist /ctl_nl/ PARTMETHOD,    &       ! Mesh partitioning method (METIS)
                      TOPOLOGY,      &       ! Mesh topology
                     qsize,         &       ! number of SE tracers
                     nthreads,      &       ! Number of threads per process
                     vert_num_threads,      &       ! Number of threads per process
                     nthreads_accel,      &       ! Number of threads per an accelerator process
                     limiter_option, &
                     smooth,        &        ! Timestep Filter
                     omega,         &
                     pertlim,        &        !temperature initial perturbation
                     tracer_transport_method, &
                     npart,         &
                     uselapi,       &
                     multilevel,    &
                     useframes,     &
                     numnodes,      &
                     ne,            &       ! element resolution factor
                     tasknum,       &
                     remapfreq,     &       ! number of steps per remapping call
                     remap_type,    &       ! selected remapping option
                     statefreq,     &       ! number of steps per printstate call
                     accumfreq,     &       ! frequency in steps of accumulation
                     accumstart,    &       ! model day to start accumulating state variables
                     accumstop,     &       ! model day to stop  accumulating state variables
                     integration,   &       ! integration method
                     tracer_advection_formulation, &
                     tstep_type, &
                     npdg, &
                     compute_mean_flux, &
                     cubed_sphere_map, &
                     qsplit, &
                     rsplit, &
                     physics, &             ! The type of physics, 0=none, 1=multicloud or 2= emanuel.
                     rk_stage_user, &
                     LFTfreq,       &
                     disable_diagnostics, &
                     prescribed_wind, &
                     se_ftype,        &       ! forcing type
                     energy_fixer,        &       ! forcing type
                     fine_ne,       &
                     max_hypervis_courant, &
                     nu,            &
                     nu_s,          &
                     nu_q,          &
                     nu_div,          &
                     nu_p,          &
                     nu_top,        &
                     psurf_vis,    &
                     hypervis_order,    &
                     hypervis_power,    &
                     hypervis_subcycle, &
                     hypervis_subcycle_q, &
                     hypervis_scaling, &
                     smooth_phis_numcycle, &
                     smooth_sgh_numcycle, &
                     smooth_phis_nudt, &
                     initial_total_mass, &
                     u_perturb,     &
                     rotate_grid,   &
                     mesh_file,     &               ! Name of mesh file
                     vert_remap_q_alg, &
                     se_prescribed_wind_2d

    namelist /ctl_nl/test_case,       &       ! test case
                     sub_case,        &       ! generic test case parameter
                     nmax,            &       ! number of steps
                     ndays,           &       ! number of days to simulate
                     restartfreq,     &
                     restartfile,     &       ! name of the restart file for INPUT
                     restartdir,      &       ! name of the restart directory for OUTPUT
                     runtype,         &
                     tstep,           &       ! tracer time step
                     moisture

    namelist /filter_nl/filter_type,   &
                        transfer_type, &
                        filter_freq,   &
                        filter_mu,     &
                        filter_freq_advection,   &
                        filter_mu_advection,     &
                        p_bv,          &
                        s_bv,          &
                        wght_fm,       &
                        kcut_fm       

    namelist /vert_nl/vform,           &
                      vfile_mid,       &
                      vfile_int

    namelist /analysis_nl/    &
         output_prefix,       &
         output_timeunits,    &
         output_start_time,   &
         output_end_time,     &
         output_frequency,    &
         output_dir,          &
         output_type,         &
         io_stride,           &
         num_io_procs,        &
         infilenames,         &
         output_varnames1,    &
         output_varnames2,    &
         output_varnames3,    &
         output_varnames4,    &
         output_varnames5

    ! ==========================
    ! Set the default partmethod
    ! ==========================

    PARTMETHOD    = RECURSIVE
    npart         = 1
    useframes     = 0
    multilevel    = 1
    uselapi       = .TRUE.

    ndays         = 0
    nmax          = 12
    nthreads = 1
    vert_num_threads = 1
    nthreads_accel = -1
    se_ftype = ftype   ! MNL: For non-CAM runs, ftype=0 in control_mod
    phys_tscale=0
    nsplit = 1
    pertlim = 0.0_real_kind

    sub_case      = 1
    numnodes      = -1
    restartfreq   = -100
    restartdir    = "./restart/"
    runtype       = 0
    statefreq     = 1
    remapfreq     = 240
    remap_type    = "parabolic"
    tasknum       =-1
    integration   = "explicit"
    moisture      = "dry"
    nu_top=0
    initial_total_mass=0
    mesh_file='none'
    ne              = 0

    disable_diagnostics = .false.


    ! =======================
    ! Read namelist variables
    ! =======================

!   write(iulog,*) "masterproc=",par%masterproc

    if (par%masterproc) then

       write(iulog,*)"reading ctl namelist..."

#if defined(OSF1) || defined(_BGL) || defined(_NAMELIST_FROM_FILE)
       open(unit=7,file="input.nl",status="OLD")
       read(unit=7,nml=ctl_nl)
#else
      print *, 'Reading namelist ctl_nl from standard input'
      read(*,nml=ctl_nl)
#endif
#ifndef _USEMETIS
      !=================================
      ! override the selected partition
      ! method and set it to SFCURVE
      !=================================
      PARTMETHOD = SFCURVE
#endif
       ! ========================
       ! if this is a restart run
       ! ========================
       if(runtype .eq. 1) then
          write(iulog,*)"readnl: restartfile = ",restartfile
       else if(runtype < 0) then
          write(iulog,*)'readnl: runtype=', runtype,' interpolation mode '
       endif



       ! ================================================
       ! if ndays is defined in the namelist, then
       ! moviefreq and restartfreq are interpreted to be in units of days.
       ! Both must be converted to numbers of steps.
       ! ================================================
       if (tstep <= 0) then
          call abortmp('tstep must be > 0')
       end if
       if (ndays .gt. 0) then
          nmax = ndays * (secpday/tstep)
          restartfreq  = restartfreq*(secpday/tstep)
       end if
       nEndStep = nmax

       if((integration .ne. "explicit").and.(integration .ne. "runge_kutta")) then
          call abortmp('integration must be explicit or runge_kutta')
       end if

       write(iulog,*)"reading filter namelist..."
       ! Set default mu/freq for advection filtering
       filter_mu_advection   = 0.05_real_kind
       filter_freq_advection = 0
       filter_freq=0

#if defined(OSF1) || defined(_BGL) || defined(_NAMELIST_FROM_FILE)
       read(unit=7,nml=filter_nl)
#else
       read(*,nml=filter_nl)
#endif
       !
       ! A modulo(a,p) where p == 0 is undefined
       if(filter_freq == 0) filter_freq = -1

       write(iulog,*)"reading physics namelist..."
       if (test_case(1:5)=="dcmip" .or. &
           test_case(1:10)=="baroclinic" .or. &
           test_case(1:13)=="jw_baroclinic" .or. &
           test_case(1:12)=="dcmip2_schar" .or. &
           test_case(1:4)=="asp_") then
         write(iulog,*) "reading vertical namelist..."
#if defined(OSF1) || defined(_BGL) || defined(_NAMELIST_FROM_FILE)
         read(unit=7,nml=vert_nl)
#else
         read(*,nml=vert_nl)
#endif
         ! Reformat these strings
         vform = trim(adjustl(vform))
         vfile_mid = trim(adjustl(vfile_mid))
         vfile_int = trim(adjustl(vfile_int))

         write(iulog,*) '  vform =',vform
         write(iulog,*) '  vfile_mid=',vfile_mid
         write(iulog,*) '  vfile_int=',vfile_int

         write(iulog,*)"reading aquaplanet namelist..."

       end if

       infilenames(:) = ''
       output_prefix = ""
       output_start_time=0
       output_end_time =0
       output_frequency=0
       output_timeunits=0
! Default is to write all variables at time 0 and time nendstep.
       output_end_time(1) = -1
       output_frequency(1) = nendstep
       output_dir = "./movies/"
       output_varnames1='     '
       output_varnames2='     '
       output_varnames3='     '
       output_varnames4='     '
       output_varnames5='     '
       io_stride=0
       num_io_procs=0
       output_type = 'netcdf' ! Change by MNL
!     output_type = 'pnetcdf'

       write(iulog,*)"reading analysis namelist..."

#if defined(OSF1) || defined(_BGL) || defined(_NAMELIST_FROM_FILE)
       read(unit=7,nml=analysis_nl)
#else
       read(*,nml=analysis_nl)
#endif

      if (io_stride .eq.0 .and. num_io_procs .eq.0) then
         ! user did not set anything
         io_stride=1
         num_io_procs = par%nprocs
      else if (io_stride.eq.0) then  ! user set num_io_procs
         io_stride = max(1,par%nprocs/num_io_procs)
      else if (num_io_procs .eq. 0) then   ! user set io_stride
         num_io_procs=max(1,par%nprocs/io_stride)
      else ! user set both parameters
      endif
      ! sanity check
      if(num_io_procs*io_stride>par%nprocs) then
         if (io_stride > par%nprocs) io_stride=par%nprocs
         num_io_procs=par%nprocs/io_stride
      end if

      if(output_varnames1(1).eq.'     ') then
         if( runtype>=0) then
            call setvarnames(output_varnames1)
         else ! interpolation mode
            output_varnames1='all'
         end if
      end if
       do i=1,max_output_streams
          if(output_frequency(i)>0 .and. output_end_time(i)==0) output_end_time(i)=-1
          if(output_timeunits(i).eq.1) then  ! per_day
             output_frequency(i) = output_frequency(i)*(secpday/tstep)
             output_start_time(i)= output_start_time(i)*(secpday/tstep)
             output_end_time(i)  = output_end_time(i)*(secpday/tstep)
          else if(output_timeunits(i).eq.2) then  ! per_hour
             output_frequency(i) = output_frequency(i)*(secphr/tstep)
             output_start_time(i)= output_start_time(i)*(secphr/tstep)
             output_end_time(i)  = output_end_time(i)*(secphr/tstep)
          end if
          if(output_end_time(i)<0) then
             output_end_time(i)=nEndStep
          endif
       end do

#if defined(OSF1) || defined(_BGL) || defined(_NAMELIST_FROM_FILE)
       close(unit=7)
#endif
    end if

    call MPI_barrier(par%comm,ierr)

    npart  = par%nprocs

    ! =====================================
    !  Spread the namelist stuff around
    ! =====================================

    call MPI_bcast(PARTMETHOD ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(TOPOLOGY     ,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(test_case    ,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(tasknum ,1,MPIinteger_t,par%root,par%comm,ierr)

    call MPI_bcast( ne        ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(qsize     ,1,MPIinteger_t,par%root,par%comm,ierr)

    call MPI_bcast(sub_case ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(remapfreq ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(remap_type, MAX_STRING_LEN, MPIChar_t, par%root, par%comm, ierr)
    call MPI_bcast(statefreq ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(accumfreq ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(accumstart,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(accumstop ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(restartfreq,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(multilevel ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(useframes ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(runtype   ,1,MPIinteger_t,par%root,par%comm,ierr)

    call MPI_bcast(omega     ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(pertlim   ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(tstep     ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(nmax      ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(NTHREADS  ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(vert_num_threads,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(nthreads_accel  ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(ndays     ,1,MPIinteger_t,par%root,par%comm,ierr)

    nEndStep = nmax
call MPI_bcast(smooth    ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(phys_tscale,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(NSPLIT,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(limiter_option  ,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(se_ftype     ,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(energy_fixer,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(vert_remap_q_alg,1,MPIinteger_t   ,par%root,par%comm,ierr)

    call MPI_bcast(fine_ne    ,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(max_hypervis_courant,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(nu         ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(nu_s         ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(nu_q         ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(nu_div       ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(nu_p         ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(nu_top   ,1,MPIreal_t   ,par%root,par%comm,ierr)

    call MPI_bcast(disable_diagnostics,1,MPIlogical_t,par%root,par%comm,ierr)
    call MPI_bcast(psurf_vis,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_order,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_power,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_scaling,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_subcycle,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(hypervis_subcycle_q,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_phis_numcycle,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_sgh_numcycle,1,MPIinteger_t   ,par%root,par%comm,ierr)
    call MPI_bcast(smooth_phis_nudt,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(initial_total_mass ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(u_perturb     ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(rotate_grid   ,1,MPIreal_t   ,par%root,par%comm,ierr)
    call MPI_bcast(integration,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)
    call MPI_bcast(mesh_file,MAX_FILE_LEN,MPIChar_t ,par%root,par%comm,ierr)
    call MPI_bcast(tracer_advection_formulation,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(tstep_type,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(npdg,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(compute_mean_flux,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(cubed_sphere_map,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(qsplit,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(rsplit,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(physics,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(rk_stage_user,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(LFTfreq,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(prescribed_wind,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(moisture,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)

    call MPI_bcast(restartfile,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)
    call MPI_bcast(restartdir,MAX_STRING_LEN,MPIChar_t ,par%root,par%comm,ierr)

    call MPI_bcast(uselapi,1,MPIlogical_t,par%root,par%comm,ierr)

    call MPI_bcast(filter_type   ,8,MPIChar_t    ,par%root,par%comm,ierr)
    call MPI_bcast(transfer_type ,8,MPIChar_t    ,par%root,par%comm,ierr)
    call MPI_bcast(filter_mu     ,1,MPIreal_t    ,par%root,par%comm,ierr)
    call MPI_bcast(filter_freq   ,1,MPIinteger_t ,par%root,par%comm,ierr)
    call MPI_bcast(filter_mu_advection     ,1,MPIreal_t    ,par%root,par%comm,ierr)
    call MPI_bcast(filter_freq_advection   ,1,MPIinteger_t ,par%root,par%comm,ierr)

    if (transfer_type == "bv") then
       call MPI_bcast(p_bv      ,1,MPIreal_t    ,par%root,par%comm,ierr)
       call MPI_bcast(s_bv      ,1,MPIreal_t    ,par%root,par%comm,ierr)
    else if (transfer_type == "fm") then
       call MPI_bcast(kcut_fm   ,1,MPIinteger_t,par%root,par%comm,ierr)
       call MPI_bcast(wght_fm   ,1,MPIreal_t   ,par%root,par%comm,ierr)
    end if

    call MPI_bcast(vform    ,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(vfile_mid,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(vfile_int,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(se_prescribed_wind_2d,1 ,MPIlogical_t  ,par%root,par%comm,ierr)

    call MPI_bcast(output_prefix,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(output_timeunits ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(output_start_time ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(output_end_time ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(output_frequency ,max_output_streams,MPIinteger_t,par%root,par%comm,ierr)

    call MPI_bcast(output_dir   ,MAX_STRING_LEN,MPIChar_t  ,par%root,par%comm,ierr)
    call MPI_bcast(output_varnames1 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(output_varnames2 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(output_varnames3 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(output_varnames4 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(output_varnames5 ,varname_len*max_output_variables,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(io_stride , 1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(num_io_procs , 1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(output_type , 9,MPIChar_t,par%root,par%comm,ierr)
    call MPI_bcast(infilenames ,160*MAX_INFILES ,MPIChar_t,par%root,par%comm,ierr)

! Set and broadcast tracer transport type
    if (trim(tracer_transport_method) == 'se_gll') then
      tracer_transport_type = TRACERTRANSPORT_SE_GLL
      tracer_grid_type = TRACER_GRIDTYPE_GLL
    else
      call abortmp('Unknown tracer transport method: '//trim(tracer_transport_method))
    end if
    call MPI_bcast(tracer_transport_type,1,MPIinteger_t,par%root,par%comm,ierr)
    call MPI_bcast(tracer_grid_type,1,MPIinteger_t,par%root,par%comm,ierr)

#ifdef IS_ACCELERATOR
    if (nthreads_accel > 0) then
        nthreads = nthreads_accel
    end if
#endif

    ! sanity check on thread count
    ! HOMME will run if if nthreads > max, but gptl will print out GB of warnings.
    if (NThreads*vert_num_threads > omp_get_max_threads()) then
       if(par%masterproc) write(iulog,*) "Main:NThreads=",NThreads
       if(par%masterproc) print *,'omp_get_max_threads() = ',OMP_get_max_threads()
       if(par%masterproc) print *,'requested threads exceeds OMP_get_max_threads()'
       call abortmp('stopping')
    endif

    ! if user sets hypervis_subcycle=-1, then use automatic formula
    if (hypervis_subcycle==-1) then
       if (np==4) then
          ! 1.25d23 worked out by testing, for nv==4
          ! a little confusing:
          ! u,v:  nu and hypervis_subcycle
          ! T:    nu_s and hypervis_subcycle
          ! Q:    nu and hypervis_subcycle_q
          dt_max = 1.25d23/(nu*ne**4.0)
          hypervis_subcycle_q = ceiling( tstep/dt_max )
          hypervis_subcycle   = ceiling( tstep/dt_max )
       else
          call abortmp('hypervis_subcycle auto determine only supported for nv==4')
       endif
    endif

    if (mesh_file /= "none" .AND. ne /=0) then
      write (*,*) "namelist_mod: mesh_file:",mesh_file, &
                  " and ne:",ne," are both sepcified in the input file."
      write (*,*) "Specify one or the other, but not both."
      call abortmp("Do not specify ne if using a mesh file input.")
    end if
    if (par%masterproc) write (iulog,*) "Mesh File:", trim(mesh_file)
    if (ne.eq.0) then
       if (par%masterproc) write (iulog,*) "Opening Mesh File:", trim(mesh_file)
      call set_mesh_dimensions()
      call MeshOpen(mesh_file, par)
    end if
    ! set map
    if (cubed_sphere_map<0) then
       cubed_sphere_map=0  ! default is equi-angle gnomonic
       if (ne.eq.0) cubed_sphere_map=2  ! element_local for var-res grids
    endif
    if (par%masterproc) write (iulog,*) "Reference element projection: cubed_sphere_map=",cubed_sphere_map

!logic around different hyperviscosity options
    if (hypervis_power /= 0) then
      if (hypervis_scaling /= 0) then
        print *,'Both hypervis_power and hypervis_scaling are nonzero.'
        print *,'(1) Set hypervis_power=1, hypervis_scaling=0 for HV based on an element area.'
        print *,'(2) Set hypervis_power=0 and hypervis_scaling=1 for HV based on a tensor.'
        print *,'(3) Set hypervis_power=0 and hypervis_scaling=0 for constant HV.'
          call abortmp("Error: hypervis_power>0 and hypervis_scaling>0")
      endif
    endif


    ftype = se_ftype

    rk_stage_user=3  ! 3d PRIM code only supports 3 stage RK tracer advection
    ! CHECK timestepping options
     if (tstep_type == 0) then
        ! pure leapfrog mode, mostly used for debugging
        if (ftype>0) then
           call abortmp('adjustment type forcing (se_ftype>0) not allowed with tstep_type=0')
        endif
        if (qsplit>1) then
          call abortmp('tracer/dynamics subcycling requires tstep_type=1(RK timestepping)')
        endif
        if (rsplit>0) then
          call abortmp('vertically lagrangian code requires tstep_type=1(RK timestepping)')
        endif
    endif


    ! CHECK phys timescale, requires se_ftype=0 (pure tendencies for forcing)
    if (phys_tscale/=0) then
       if (ftype>0) call abortmp('user specified se_phys_tscale requires se_ftype<=0')
    endif
    if (limiter_option==8 .or. limiter_option==84) then
       if (hypervis_subcycle_q/=1) then
          call abortmp('limiter 8,84 requires hypervis_subcycle_q=1')
       endif
    endif

    if((prescribed_wind/=0).and.(prescribed_wind/=1))then
          call abortmp('prescribed_wind should be either 0 or 1')
    endif

    ! some default diffusion coefficiets
    if(nu_s<0) nu_s=nu
    if(nu_q<0) nu_q=nu
    if(nu_div<0) nu_div=nu


    if (multilevel <= 0) then
      nmpi_per_node = 1
    endif

    nnodes = npart/nmpi_per_node

    if(numnodes > 0 .and. multilevel .eq. 1) then
        nnodes = numnodes
        nmpi_per_node = npart/nnodes
    endif

    ! ====================================================================
    !  Do not perform node level partitioning if you are only on one node
    ! ====================================================================
    ! PartitionForNodes=.FALSE.
    if((nnodes .eq. 1) .and. PartitionForNodes) PartitionForNodes=.FALSE.

    if (par%masterproc) then
       write(iulog,*)"done reading namelist..."

       write(iulog,*)"readnl: accum         = ",accumfreq,accumstart,accumstop

       write(iulog,*)"readnl: topology      = ",TRIM( TOPOLOGY )
       write(iulog,*)"readnl: test_case     = ",TRIM(test_case)
       write(iulog,*)"readnl: omega         = ",omega
       write(iulog,*)"readnl: sub_case      = ",sub_case
       write(iulog,*)"readnl: ndays         = ",ndays
       write(iulog,*)"readnl: nmax          = ",nmax
       write(iulog,*)"readnl: pertlim       = ",pertlim

       write(iulog,*)"readnl: qsize,qsize_d = ",qsize,qsize_d
       if (qsize>qsize_d) then
          call abortmp('user specified qsize > qsize_d parameter in dimensions_mod.F90')
       endif
       write(iulog,*)"readnl: NThreads      = ",NTHREADS
       write(iulog,*)"readnl: vert_num_threads = ",vert_num_threads
       write(iulog,*)"readnl: nthreads_accel = ",nthreads_accel

       write(iulog,*)"readnl: ne,np         = ",NE,np
       if (npdg>0) write(iulog,*)"readnl: npdg       = ",npdg
       write(iulog,*)"readnl: partmethod    = ",PARTMETHOD
       write(iulog,*)'readnl: nmpi_per_node = ',nmpi_per_node
       write(iulog,*)'readnl: multilevel    = ',multilevel
       write(iulog,*)'readnl: useframes     = ',useframes
       write(iulog,*)'readnl: nnodes        = ',nnodes
       write(iulog,*)'readnl: npart         = ',npart
       write(iulog,*)'readnl: test_cfldep   = ',test_cfldep

       print *
       write(iulog,*)"readnl: integration   = ",trim(integration)
       if (integration == "explicit" ) then
          write(iulog,*)"readnl: LF-trapazoidal freq= ",LFTfreq
       endif
       if (integration == "runge_kutta" .or. tstep_type>0 ) then
          write(iulog,*)"readnl: rk_stage_user   = ",rk_stage_user
       endif
       write(iulog,*)"readnl: tracer_advection_formulation  = ",tracer_advection_formulation
       write(iulog,*)"readnl: tstep_type    = ",tstep_type
       write(iulog,*)"readnl: vert_remap_q_alg  = ",vert_remap_q_alg

       write(iulog,*)"readnl: tstep          = ",tstep
       write(iulog,*)"readnl: ftype          = ",ftype
       write(iulog,*)"readnl: limiter_option = ",limiter_option
       write(iulog,*)"filter: smooth         = ",smooth

       write(iulog,*)"readnl: qsplit        = ",qsplit
       write(iulog,*)"readnl: vertical remap frequency rsplit (0=disabled): ",rsplit
       write(iulog,*)"readnl: physics       = ",physics

       write(iulog,*)"readnl: energy_fixer  = ",energy_fixer
       write(iulog,*)"readnl: runtype       = ",runtype

       write(iulog,*)"readnl: se_prescribed_wind_2d = ", se_prescribed_wind_2d

       if (hypervis_power /= 0)then
          write(iulog,*)"Variable scalar hyperviscosity: hypervis_power=",hypervis_power
          write(iulog,*)"max_hypervis_courant = ", max_hypervis_courant
          write(iulog,*)"Equivalent ne in fine region = ", fine_ne
       elseif(hypervis_scaling /=0)then
          write(iulog,*)"Tensor hyperviscosity:  hypervis_scaling=",hypervis_scaling
       else
          write(iulog,*)"Constant (hyper)viscosity used."
       endif

       write(iulog,*)"hypervis_subcycle, hypervis_subcycle_q = ",&
            hypervis_subcycle,hypervis_subcycle_q
       !write(iulog,*)"psurf_vis: ",psurf_vis
       write(iulog,'(a,2e9.2)')"viscosity:  nu (vor/div) = ",nu,nu_div
       write(iulog,'(a,2e9.2)')"viscosity:  nu_s      = ",nu_s
       write(iulog,'(a,2e9.2)')"viscosity:  nu_q      = ",nu_q
       write(iulog,'(a,2e9.2)')"viscosity:  nu_p      = ",nu_p
       write(iulog,'(a,2e9.2)')"viscosity:  nu_top      = ",nu_top
       write(iulog,*)"PHIS smoothing:  ",smooth_phis_numcycle,smooth_phis_nudt
       write(iulog,*)"SGH  smoothing:  ",smooth_sgh_numcycle

       if(initial_total_mass>0) then
          write(iulog,*) "initial_total_mass = ",initial_total_mass
       end if

       if (filter_freq>0 .or. filter_freq_advection>0) then
       write(iulog,*)"Filter Method is ",filter_type
       write(iulog,*)"filter: viscosity (mu)  = ",filter_mu
       write(iulog,*)"filter: frequency       = ",filter_freq
       write(iulog,*)"filter_advection: viscosity (mu)  = ",filter_mu_advection
       write(iulog,*)"filter_advection: frequency       = ",filter_freq_advection

       write(iulog,*)"filter: transfer_type   = ",transfer_type
       if (transfer_type == "bv") then
          print *
          write(iulog,*)"with Boyd-Vandeven Transfer Fn Parameters:"
          write(iulog,*)"     filter: order     (p)   = ",p_bv
          write(iulog,*)"     filter: lag coeff (s)   = ",s_bv
       else if (transfer_type == "fm") then
          print *
          write(iulog,*)"with Fischer-Mullen Transfer Fn Parameters:"
          write(iulog,*)"     filter: clipped wave nos kc = ",kcut_fm
          write(iulog,*)"     filter: amount of clipping  = ",wght_fm
       end if
       endif

       write(iulog,*)"  analysis: output_prefix = ",TRIM(output_prefix)
       write(iulog,*)"  analysis: io_stride = ",io_stride
       write(iulog,*)"  analysis: num_io_procs = ",num_io_procs

       do i=1,max_output_streams
          if(output_frequency(i) .gt. 0) then
             write(iulog,*)"  analysis stream     :",i
             write(iulog,*)"  analysis : start_time", output_start_time(i)
             write(iulog,*)"  analysis : end_time  ", output_end_time(i)
             write(iulog,*)"  analysis : frequency ", output_frequency(i)
             select case (i)
             case (1)
                write(*,'(10(a,2x))')"   analysis : variables ",(trim(output_varnames1(j)),j=1,max_output_variables)
             case (2)
                write(*,'(10(a,2x))')"   analysis : variables ",(trim(output_varnames2(j)),j=1,max_output_variables)
             case (3)
                write(*,'(10(a,2x))')"   analysis : variables ",(trim(output_varnames3(j)),j=1,max_output_variables)
             case (4)
                write(*,'(10(a,2x))')"   analysis : variables ",(trim(output_varnames4(j)),j=1,max_output_variables)
             case (5)
                write(*,'(10(a,2x))')"   analysis : variables ",(trim(output_varnames5(j)),j=1,max_output_variables)
             end select
          end if
       end do

       select case (tracer_transport_type)
       case (TRACERTRANSPORT_SE_GLL)
         write(iulog, *) 'Eulerian tracer advection on GLL grid'
        end select

    endif

  end subroutine readnl

end module namelist_mod
