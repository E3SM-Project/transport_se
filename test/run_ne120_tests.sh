#!/bin/tcsh
#PBS -q regular
#PBS -A acme
#PBS -l walltime=01:00:00
#PBS -j oe
#PBS -o out_ne120_$PBS_JOBID
#PBS -e err_ne120_$PBS_JOBID

# Be sure to change MPI/threads here and NCPU/NTHREADS below
#PBS -l mppwidth=1920

#_______________________________________________________________________
#
# This file executes a set of DCMIP test to verify tracer perfrormance
#
#  1) Edit the path to the configuration script, if necessary
#  2) Submit this script to the queue or execute it an interactive session
#_______________________________________________________________________

set NE       = 120        # number of elements per cube-edge
set NCPU     = 1920       # number of CPUs to use
set NTHREADS = 1          # number of openMP threads
set TSTEP    = 75         # timestep size
set NU       = 1e13       # hyperviscosity coefficient
set CONFIGURE_DIR = ../../


# diagnostic output every 2day
set statefreq = 48     
@ statefreq *= 3600
@ statefreq /= $TSTEP


# Run length (hours) for DCMIP1-1
# 288h (12 days) needed for verification.  expensive!
# 6h is probably good for performance testing.  
set nhours = 288
#set nhours = 1

# convert to timesteps
set nmax = $nhours
@ nmax *= 3600
@ nmax /= $TSTEP



#_______________________________________________________________________
# get path variables from configure script:

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

cd ${CONFIGURE_DIR}; source configure.sh

set QSIZE1 = $QSIZE
#set QSIZE1 = 50  # override default set in configure.sh
set QSIZE2 = $QSIZE


set TEST1_DIR = $BLD_DIR/test/dcmip1-1  # test case directory
set TEST2_DIR = $BLD_DIR/test/dcmip1-2  # test case directory
set VCOORD    = $REPO/test/vcoord       # location of vertical coordinate files

#_______________________________________________________________________
# create directory for simulation output

mkdir -p $RUN_DIR/movies
cd $RUN_DIR
cp -a $VCOORD vcoord



#_______________________________________________________________________
# create namelist for DCMIP test 1-1
cd $TEST1_DIR
sed s/NE.\*/$NE/ dcmip1-1.nl          |\
sed s/TIME_STEP.\*/$TSTEP/            |\
sed s/statefreq.\*/statefreq=$statefreq/        |\
sed s/ndays.\*/nmax=$nmax/            |\
sed s/qsize.\*/qsize=$QSIZE1/          |\
sed s/NThreads.\*/NThreads=$NTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-1_NE120.nl

#_______________________________________________________________________
# create namelist for DCMIP test 1-2
cd $TEST2_DIR
sed s/NE.\*/$NE/ dcmip1-2.nl          |\
sed s/TIME_STEP.\*/$TSTEP/            |\
sed s/statefreq.\*/statefreq=$statefreq/        |\
sed s/qsize.\*/qsize=$QSIZE2/          |\
sed s/NThreads.\*/NThreads=$NTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-2_NE120.nl

#_______________________________________________________________________
# launch the executable

cd $RUN_DIR
setenv OMP_NUM_THREADS $NTHREADS
date

# run dcmip test 1-1
echo "executing dcmip test 1-1"
echo "${RUN_COMMAND}  $NCPU $EXE < dcmip1-1_NE120.nl"
${RUN_COMMAND}  $NCPU $EXE < dcmip1-1_NE120.nl
if($status) exit
mv HommeTime_stats HommeTime_stats_DCMIP1-1_NE120
date

# run dcmip test 1-2
echo "executing dcmip test 1-2"
echo "${RUN_COMMAND} $NCPU $EXE < dcmip1-2_NE120.nl"
${RUN_COMMAND}  $NCPU $EXE < dcmip1-2_NE120.nl
if($status) exit
mv HommeTime_stats HommeTime_stats_DCMIP1-2_NE120
date


# plot results
echo
echo "Running analysis scripts on data in native-grid format"
echo
cp $TEST1_DIR/dcmip1-1_lat_lon_ng.ncl .
cp $TEST2_DIR/dcmip1-2_lat_height_ng.ncl .
ncl dcmip1-1_lat_lon.ncl
ncl dcmip1-2_lat_height_ng.ncl

# print timing info
cat HommeTime_stats_DCMIP1-1_NE120 | grep walltotal
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE120 | grep prim_run`"
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE120 | grep prim_advance_exp`"
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE120 | grep prim_advec_tracers`"
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE120 | grep vertical_remap`"
echo
echo "DCMIP1-2 `cat HommeTime_stats_DCMIP1-2_NE120 | grep prim_run`"
echo "DCMIP1-2 `cat HommeTime_stats_DCMIP1-2_NE120 | grep prim_advance_exp`"
echo "DCMIP1-2 `cat HommeTime_stats_DCMIP1-2_NE120 | grep prim_advec_tracers`"
echo "DCMIP1-2 `cat HommeTime_stats_DCMIP1-2_NE120 | grep vertical_remap`"
echo



# print error norms
cp $TEST1_DIR/dcmip1-1_error_norm_ng.ncl .
cp $TEST2_DIR/dcmip1-2_error_norm_ng.ncl .
ncl dcmip1-1_error_norm_ng.ncl  | tail -n 1
ncl dcmip1-2_error_norm_ng.ncl  | tail -n 1
echo

exit





