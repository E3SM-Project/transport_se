#!/bin/tcsh
#
#PBS -l walltime=0:05:00
#PBS -l mppwidth=48
#PBS -j oe
#PBS -o out_ne8_$PBS_JOBID
#PBS -e err_ne8_$PBS_JOBID

#_______________________________________________________________________
#
# This file executes a set of DCMIP test to verify tracer perfrormance
#
#  1) Edit the path to the configuration script, if necessary
#  2) Submit this script to the queue or execute it an interactive session
#_______________________________________________________________________

set NE       = 8          # number of elements per cube-edge
set NMPI     = 12         # total number of MPI tasks
set NTHREADS = 4          # number of openMP threads per MPI task
set TSTEP    = 1440       # timestep size
set NU       = 2e16       # hyperviscosity coefficient
set CONFIGURE_DIR = ../../

#_______________________________________________________________________
# get path variables from configure script:

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

cd ${CONFIGURE_DIR}; source configure.sh

set TEST2_DIR = $BLD_DIR/test/dcmip1-2  # test case directory
set VCOORD    = $REPO/test/vcoord       # location of vertical coordinate files

#_______________________________________________________________________
# create directory for simulation output

mkdir -p $RUN_DIR/movies
cd $RUN_DIR
cp -a $VCOORD vcoord

#_______________________________________________________________________
# create namelist for DCMIP test 1-2
cd $TEST2_DIR
sed s/NE.\*/$NE/ dcmip1-2.nl          |\
sed s/TIME_STEP.\*/$TSTEP/            |\
sed s/qsize.\*/qsize=$QSIZE/          |\
sed s/NThreads.\*/NThreads=$NTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-2_NE8.nl

#_______________________________________________________________________
# launch the executable

cd $RUN_DIR
setenv OMP_NUM_THREADS $NTHREADS
date

# run dcmip test 1-2
echo "executing dcmip test 1-2"
echo "${RUN_COMMAND} $NMPI -d $NTHREADS $EXE < dcmip1-2_NE8.nl"
${RUN_COMMAND} $NMPI -d $NTHREADS $EXE < dcmip1-2_NE8.nl
if($status) exit
mv HommeTime_stats HommeTime_stats_DCMIP1-2_NE8
date

# plot results
echo
echo "Running analysis scripts for data in native-grid format"
echo
cp $TEST2_DIR/dcmip1-2_lat_height_ng.ncl .
ncl dcmip1-2_lat_height_ng.ncl NE=$NE

# print timing info
cat HommeTime_stats_DCMIP1-2_NE8 | grep walltotal
echo "DCMIP1-2 `cat HommeTime_stats_DCMIP1-2_NE8 | grep prim_run`"
echo

# print error norms
cp $TEST2_DIR/dcmip1-2_error_norm_ng.ncl .
ncl dcmip1-2_error_norm_ng.ncl NE=$NE | tail -n 1
echo

exit




