&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "asp_tracer"
ne= 8
qsize         = 2
ndays          = 1
rotate_grid   = 0
statefreq     = 30
accumfreq     = -1
accumstart    = 200
accumstop     = 1200
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = 0
tstep= 1440
tstep_type    = 1
qsplit        = 1
integration   = "explicit"
smooth        = 0.00               ! disabled
nu_q= 2e16
limiter_option = 8
hypervis_order = 2
hypervis_subcycle= 1
prescribed_wind = 1
energy_fixer = -1
/
&solver_nl
precon_method = "identity"
maxits        = 500
tol           = 1.e-9
/
&filter_nl
filter_type   = "taylor"
transfer_type = "bv"
filter_freq   = 0
filter_mu     = 0.04D0
filter_freq_advection   = 0
filter_mu_advection   = 0.00
p_bv          = 12.0D0
s_bv          = .80
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
!vfile_mid     = "vcoord/aspL60_mid.ascii"
!vfile_int     = "vcoord/aspL60_int.ascii"
vfile_mid     = "vcoord/sabm-64.ascii"
vfile_int     = "vcoord/sabi-64.ascii"
/
&analysis_nl
 interp_gridtype=1
! every 6h
 output_timeunits=2
 output_frequency=24
 output_start_time=0
 output_end_time=240
 output_varnames1='Q','Q2'
 io_stride=8
 output_type = 'netcdf'
/


&prof_inparm
profile_outpe_num = 512
profile_single_file		= .true.
/
