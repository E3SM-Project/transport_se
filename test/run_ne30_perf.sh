#!/bin/tcsh
#PBS -q regular
#PBS -l walltime=0:05:00
#PBS -l mppwidth=21600
#PBS -j oe
#PBS -o out_ne30_perf_$PBS_JOBID
#PBS -e err_ne30_perf_$PBS_JOBID

# Runs DCMIP 1-1 with NE=30 for a few model hours and measures run times.
# mppwdith = (Ne*Ne*6)*threads_per_colum = 30*30*6*4 = 21600

#_______________________________________________________________________
# call configure script to ensure executable exists

set CONFIGURE_DIR = ../..     # location of configure script

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif
cd ${CONFIGURE_DIR}; source configure.sh

#_______________________________________________________________________
# set test parameters

set HTHREADS  = 1              # number of horizontal threads
set VTHREADS  = 4              # number of vertical threads (column_omp)

set NE        = 30             # number of elements per cube-edge
set TSTEP     = 300            # time step size, in second
set NU        = 1e15           # hyperviscosity coefficient
set QSIZE     = 35             # number of tracers
set NHOURS    = 144            # total simulation time (12 days)

set TEST_NAME = run_ne30_perf  # name of test for run directory

#_______________________________________________________________________
# compute run parameters from number of procs and number of threads

if ( ${?PBS_NP} == 0) then
  set PBS_NP = 24;                                # set default NP
endif
@ NUM_NODES     = $PBS_NP / 24                    # compute number of nodes from mppwidth
@ NTHREADS      = $HTHREADS * $VTHREADS           # get total number of threads needed
@ NCPU          = $NUM_NODES * 24 / $NTHREADS     # get total number of MPI procs
@ NCPU_PER_NODE = 24 / $NTHREADS                  # get number of MPI procs per node
@ NUM_NUMA      = $NCPU_PER_NODE / 2              # edison has 2 sockets per node
@ statefreq     = 144 * 3600 / $TSTEP             # set diagnostic display frequency
@ nmax          = $NHOURS * 3600 / $TSTEP         # get max number of timesteps

#set RUN_COMMAND = "aprun -n $NCPU -N $NCPU_PER_NODE -d $NTHREADS -S $NUM_NUMA -ss -cc numa_node"
set RUN_COMMAND = "aprun -n $NCPU -N $NCPU_PER_NODE -d $NTHREADS -S $NUM_NUMA -ss -cc cpu"

setenv OMP_NUM_THREADS $NTHREADS

echo "PBS_NP        = $PBS_NP"
echo "NUM_NODES     = $NUM_NODES"
echo "NTHREADS      = $NTHREADS"
echo "NUM_CPU       = $NCPU"
echo "NCPU_PER_NODE = $NCPU_PER_NODE"
echo "NUM_NUMA      = $NUM_NUMA"
echo "statefreq     = $statefreq"

#_______________________________________________________________________
# check for some common errors

set OMP_STATUS        = `cat $BLD_DIR/CMakeCache.txt | grep ENABLE_OPENMP        | grep TRUE`
set COLUMN_OMP_STATUS = `cat $BLD_DIR/CMakeCache.txt | grep ENABLE_COLUMN_OPENMP | grep TRUE`

if( $NTHREADS > 1 && ${%OMP_STATUS} == 0 ) then
  echo "Error: NTHREADS > 1 requires ENABLE_OPENMP=TRUE"; exit
endif

if( $VTHREADS > 1 && ${%COLUMN_OMP_STATUS} == 0 ) then
  echo "Error: VTHREADS > 1 requires ENABLE_COLUMN_OPENMP=TRUE"; exit
endif

#_______________________________________________________________________
# create directories for simulation output

if ( ${?PBS_JOBID} == 0) then
  set PBS_JOBID=`date "+%y%m%d%H%M%S"`
endif

set RUN_DIR = $WORK/${TEST_NAME}_$PBS_JOBID
mkdir -p $RUN_DIR/movies
cd $RUN_DIR
cp -a $REPO/test/vcoord vcoord

echo "RUN_DIR: $RUN_DIR"

#_______________________________________________________________________
# create namelist for DCMIP test 1-1

set TEST1_DIR  = $BLD_DIR/test/dcmip1-1

cd $TEST1_DIR
sed s/NE.\*/$NE/ dcmip1-1.nl                    |\
sed s/TIME_STEP.\*/$TSTEP/                      |\
sed s/statefreq.\*/statefreq=$statefreq/        |\
sed s/ndays.\*/nmax=$nmax/                      |\
sed s/qsize.\*/qsize=$QSIZE/                    |\
sed s/NThreads.\*/NThreads=$HTHREADS/           |\
sed s/vert_num_threads.\*/vert_num_threads=$VTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-1.nl

#_______________________________________________________________________
# launch the executable

cd $RUN_DIR
date

# run dcmip test 1-1
echo "${RUN_COMMAND} $EXE < dcmip1-1.nl"
${RUN_COMMAND} $EXE < dcmip1-1.nl
if($status) exit

date

# print timing info
cat HommeTime_stats | grep walltotal
echo "DCMIP1-1 `cat HommeTime_stats | grep prim_run`"
echo "DCMIP1-1 `cat HommeTime_stats | grep prim_advance_exp`"
echo "DCMIP1-1 `cat HommeTime_stats | grep prim_advec_tracers`"
echo "DCMIP1-1 `cat HommeTime_stats | grep vertical_remap`"
echo

# print error norms
cp $TEST1_DIR/dcmip1-1_error_norm_ng.ncl .
ncl dcmip1-1_error_norm_ng.ncl | tail -n 1
echo

exit




