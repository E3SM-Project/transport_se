#!/bin/tcsh
#PBS -q debug
#PBS -l walltime=0:05:00
#PBS -l mppwidth=1536
#PBS -j oe
#PBS -o out_ne8_$PBS_JOBID
#PBS -e err_ne8_$PBS_JOBID

#_______________________________________________________________________
#
# This file executes a set of DCMIP test to verify tracer perfrormance
#
#  1) Edit the path to the configuration script, if necessary
#  2) Submit this script to the queue or execute it an interactive session
#________________
#  _______________________________________________________
set NCPU     = 384     # number of CPUs to use
set NTHREADS = 4          # number of openMP threads per MPI task
set NCPU_PER_NODE = 6     # number of MPI tasks per node

set NE       = 8          # number of elements per cube-edge
set TSTEP    = 1200       # note: remap time should not exceed 1800, so set rsplit=1
set NU       = 6e16       # hyperviscosity coefficient
set CONFIGURE_DIR = ../../

set QSIZE1 = 50

# edison has 2 sockets per node
set NUM_NUMA = $NCPU_PER_NODE
@ NUM_NUMA /= 2



# diagnostic output every so-many hours
set statefreq = 144     
@ statefreq *= 3600
@ statefreq /= $TSTEP

#_______________________________________________________________________
# get path variables from configure script:

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

cd ${CONFIGURE_DIR}; source configure.sh

# override default RUN_COMMAND set in configure.sh
set RUN_COMMAND="aprun -n $NCPU -N $NCPU_PER_NODE -d $NTHREADS -S $NUM_NUMA -ss -cc numa_node "
echo RUN_COMMAND:
echo $RUN_COMMAND


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
sed s/NThreads.\*/NThreads=$NTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-1_NE8.nl

#_______________________________________________________________________
# launch the executable

cd $RUN_DIR
setenv OMP_NUM_THREADS $NTHREADS
date

# run dcmip test 1-1
echo "executing dcmip test 1-1"
echo "${RUN_COMMAND}  $NCPU $EXE < dcmip1-1_NE8.nl"
${RUN_COMMAND}  $EXE < dcmip1-1_NE8.nl
if($status) exit
mv HommeTime_stats HommeTime_stats_DCMIP1-1_NE8
date

# print timing info
cat HommeTime_stats_DCMIP1-2_NE8 | grep walltotal
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




