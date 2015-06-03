&ctl_nl
  NThreads      = 1
  partmethod    = 4
  topology      = "cube"                  ! initial mesh type
  test_case     = "dcmip1-2"
  ne            = NE                       ! elements per cube face
  qsize         = 2                       ! num tracer fields
  ndays         = 1                       ! num simulation days, 0=>use nmax
  statefreq     = 10
  accumfreq     = -1
  accumstart    = 200
  accumstop     = 1200
  restartfreq   = -100
  restartfile   = "./R0001"
  runtype       = 0
  tstep         = TIME_STEP
  tstep_type    = 1
  qsplit            = 1
  rsplit            = 0
  integration       = 'explicit'
  smooth            = 0
  nu                = 0
  nu_p              = 0
  nu_s              = 0
  nu_q              = 0
  limiter_option    = 8
  hypervis_order    = 2
  hypervis_subcycle = 1
  prescribed_wind   = 1
  energy_fixer      = -1
/
&filter_nl
  filter_type   = "taylor"
  transfer_type = "bv"
  filter_freq   = 0
  filter_mu     = 0.04
  p_bv          = 12.0D0
  s_bv          = .80
  wght_fm       = 0.10D0
  kcut_fm       = 2
/
&vert_nl
  vform             = "ccm"
  vfile_mid         = "vcoord/12k_top-64m.ascii"
  vfile_int         = "vcoord/12k_top-64i.ascii"
/
&analysis_nl
 interp_gridtype  = 1
 output_dir       = "./movies/"                   ! destination dir for netcdf file
 output_timeunits = 2,                            ! 1=days, 2=hours, 0=timesteps
 output_frequency = 12,                           ! interval between outputs
 output_varnames1 = 'Q2','geo'                    ! variables to write to file
 interp_type      = 0                             ! 0=native grid, 1=bilinear
 output_type      ='netcdf'                       ! netcdf or pnetcdf
 io_stride        = 8
/
&prof_inparm
  profile_outpe_num   = 512
  profile_single_file	= .true.
/