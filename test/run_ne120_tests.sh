#!/bin/tcsh
#
#PBS -l walltime=0:10:00
#PBS -l mppwidth=2048
#PBS -j oe
#PBS -o out_ne120_$PBS_JOBID
#PBS -e err_ne120_$PBS_JOBID

#_______________________________________________________________________
#
# This file executes a set of DCMIP test to verify tracer perfrormance
#
#  1) Edit the path to the configuration script, if necessary
#  2) Submit this script to the queue or execute it an interactive session
#_______________________________________________________________________

set NE       = 120        # number of elements per cube-edge
set NCPU     = 2048       # number of CPUs to use
set NTHREADS = 1          # number of openMP threads
set TSTEP    = 75         # timestep size
set NU       = 1e13       # hyperviscosity coefficient
set CONFIGURE_DIR = ../../

#_______________________________________________________________________
# get path variables from configure script:

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

cd ${CONFIGURE_DIR}; source configure.sh

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
sed s/qsize.\*/qsize=$QSIZE/          |\
sed s/NThreads.\*/NThreads=$NTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-1_NE120.nl

#_______________________________________________________________________
# create namelist for DCMIP test 1-2
cd $TEST2_DIR
sed s/NE.\*/$NE/ dcmip1-2.nl          |\
sed s/TIME_STEP.\*/$TSTEP/            |\
sed s/qsize.\*/qsize=$QSIZE/          |\
sed s/NThreads.\*/NThreads=$NTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-2_NE120.nl

#_______________________________________________________________________
# launch the executable

cd $RUN_DIR
setenv OMP_NUM_THREADS $NTHREADS
date

# run the shorter dcmip test 1-2 first
echo "running dcmip test 1-2"
${RUN_COMMAND} $NCPU $EXE < dcmip1-2_NE120.nl
if($status) exit

mv HommeTime_stats HommeTime_stats_DCMIP1-2_NE120
date

# run longer dcmip test 1-1
echo "running dcmip test 1-1"
${RUN_COMMAND}  $NCPU $EXE < dcmip1-1_NE120.nl
if($status) exit

mv HommeTime_stats HommeTime_stats_DCMIP1-1_NE120
date

# post-processes the output
echo "post-processing the data"
date

exit



