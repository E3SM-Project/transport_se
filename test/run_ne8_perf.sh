#!/bin/tcsh
#PBS -q debug
#PBS -l walltime=0:05:00
#PBS -l mppwidth=384
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


set HTHREADS = 2          # number of horizontal threads
set VTHREADS = 2          # number of vertical threads

set NE       = 8          # number of elements per cube-edge
set TSTEP    = 1200       # note: remap time should not exceed 1800, so set rsplit=1
set NU       = 6e16       # hyperviscosity coefficient
set QSIZE1   = 50
set CONFIGURE_DIR = ../../

@ NUM_NODES     = $PBS_NP / 24                    # compute number of nodes from mppwidth
@ NUM_THREADS   = $HTHREADS * $VTHREADS           # get total number of threads needed
@ NCPU          = $NUM_NODES * 24 / $NUM_THREADS  # get total number of MPI procs
@ NCPU_PER_NODE = 24 / $NUM_THREADS               # get number of MPI procs per node
@ NUM_NUMA      = $NCPU_PER_NODE / 2              # edison has 2 sockets per node
@ statefreq     = 144 * 3600 / $TSTEP             # set diagnostic display frequency

echo "NUM_NODES     = $NUM_NODES"
echo "NUM_THREADS   = $NUM_THREADS"
echo "NUM_CPU       = $NCPU"
echo "NCPU_PER_NODE = $NCPU_PER_NODE"
echo "NUM_NUMA      = $NUM_NUMA"
echo "statefreq     = $statefreq"

#_______________________________________________________________________
# get path variables from configure script:

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

cd ${CONFIGURE_DIR}; source configure.sh

# override default RUN_COMMAND set in configure.sh
set RUN_COMMAND="aprun -n $NCPU -N $NCPU_PER_NODE -d $NUM_THREADS -S $NUM_NUMA -ss -cc numa_node"

echo "RUN_COMMAND: $RUN_COMMAND"

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
sed s/qsize.\*/qsize=$QSIZE1/          |\
sed s/rsplit.\*/rsplit=1/          |\
sed s/NThreads.\*/NThreads=$HTHREADS/ |\
sed s/vert_num_threads.\*/vert_num_threads=$VTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-1_NE8.nl

#_______________________________________________________________________
# launch the executable

cd $RUN_DIR
setenv OMP_NUM_THREADS $NUM_THREADS
#setenv OMP_NESTED TRUE
date

# run dcmip test 1-1
echo "executing dcmip test 1-1"
echo "${RUN_COMMAND} $EXE < dcmip1-1_NE8.nl"
${RUN_COMMAND}  $EXE < dcmip1-1_NE8.nl
if($status) exit
mv HommeTime_stats HommeTime_stats_DCMIP1-1_NE8
date

# print timing info
cat HommeTime_stats_DCMIP1-1_NE8 | grep walltotal
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE8 | grep prim_run`"
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE8 | grep prim_advance_exp`"
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE8 | grep prim_advec_tracers`"
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE8 | grep vertical_remap`"
echo

# print error norms
cp $TEST1_DIR/dcmip1-1_error_norm_ng.ncl .
ncl dcmip1-1_error_norm_ng.ncl | tail -n 1
echo

exit




