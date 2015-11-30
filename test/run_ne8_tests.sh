#!/bin/tcsh
#PBS -q debug
#PBS -l mppwidth=48
#PBS -l walltime=00:10:00
#PBS -A acme
#SBATCH -p debug
#SBATCH -t 00:05:00
#SBATCH -N 2

# Runs DCMIP 1-1 and 1-2 and plots results
#_______________________________________________________________________
# call configure script to ensure executable exists
set CONFIGURE_DIR = ../..     # location of configure script

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif
cd ${CONFIGURE_DIR}; source configure.sh

#_______________________________________________________________________
# set test parameters

set TEST_NAME = run_ne8_tests # name of test for run directory
set NE        = 8             # number of elements per cube-edge
set TSTEP     = 400           # time step size, in second.  
                              # Should be TSTEP=1200, but that sets remap=3600 which is too long
set NU        = 6e16          # hyperviscosity coefficient
set QSIZE     = 4             # number of tracers
@ statefreq   = 12 * 3600 / $TSTEP            # set diagnostic display frequency

#_______________________________________________________________________
# compute run parameters from number of procs and number of threads

set HTHREADS  = 1             # number of horizontal threads
set VTHREADS  = 1             # number of vertical threads (column_omp)
@ NTHREADS    = $HTHREADS * $VTHREADS           # get total number of threads needed
setenv OMP_NUM_THREADS $NTHREADS


#set MAX_TASKS_NODE = 24       # Edison.  48 with hyperhtreading
set MAX_TASKS_NODE = 32        # Cori     64 with hyperhtreading
#set MAX_TASKS_NODE = 64       # Mira




set NNODES = 1
if ( ${?SLURM_NNODES} ) then
   set NNODES = $SLURM_NNODES
endif
if ( ${?SLURM_JOB_NUM_NODES} ) then
   set NNODES = $SLURM_JOB_NUM_NODES
endif
if ( ${?PBS_NP} ) then
  # hack for edison because PBS_NUM_NODES is always 1
  @ NNODES = $PBS_NP / 24
endif


@ NMPI = $NNODES * $MAX_TASKS_NODE / $NTHREADS
@ NMPI_PER_NODE = $NMPI / $NNODES              # get number of MPI procs per node
@ NUM_NUMA      = $NMPI_PER_NODE / 2           # edison has 2 sockets per node

#set RUN_COMMAND = "aprun -n $NMPI -N $NMPI_PER_NODE -d $NTHREADS -S $NUM_NUMA -ss -cc numa_node"
set RUN_COMMAND = "srun -n $NMPI -c $NTHREADS"
#set RUN_COMMAND = "mpirun -n $NMPI"

echo "NNODES        = $NNODES"
echo "NMPI          = $NMPI"
echo "NMPI_PER_NODE = $NMPI_PER_NODE"
echo "NTHREADS      = $NTHREADS"
echo "HTHREADS      = $HTHREADS"
echo "VTHREADS      = $VTHREADS"
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

set JDATE=`date "+%y-%m-%d_%H%M%S"`
set RUN_DIR = $WORK/${TEST_NAME}_$JDATE
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
sed s/qsize.\*/qsize=$QSIZE/                    |\
sed s/NThreads.\*/NThreads=$HTHREADS/           |\
sed s/vert_num_threads.\*/vert_num_threads=$VTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-1.nl

#_______________________________________________________________________
# create namelist for DCMIP test 1-2

set TEST2_DIR  = $BLD_DIR/test/dcmip1-2

cd $TEST2_DIR
sed s/NE.\*/$NE/ dcmip1-2.nl                    |\
sed s/TIME_STEP.\*/$TSTEP/                      |\
sed s/statefreq.\*/statefreq=$statefreq/        |\
sed s/qsize.\*/qsize=$QSIZE/                    |\
sed s/NThreads.\*/NThreads=$HTHREADS/           |\
sed s/vert_num_threads.\*/vert_num_threads=$VTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-2.nl

#_______________________________________________________________________
# launch the executable

cd $RUN_DIR
date

# run dcmip test 1-1
echo "executing dcmip test 1-1"
echo "${RUN_COMMAND} $EXE < dcmip1-1.nl"
${RUN_COMMAND} $EXE < dcmip1-1.nl
if($status) exit
mv HommeTime_stats HommeTime_stats_DCMIP1-1
date

# run dcmip test 1-2
echo "executing dcmip test 1-2"
echo "${RUN_COMMAND} $EXE < dcmip1-2.nl"
${RUN_COMMAND} $EXE < dcmip1-2.nl
if($status) exit
mv HommeTime_stats HommeTime_stats_DCMIP1-2
date

# plot results
echo
echo "Running analysis scripts on data in native-grid format"
echo
cp $TEST1_DIR/dcmip1-1_lat_lon_ng.ncl .
cp $TEST2_DIR/dcmip1-2_lat_height_ng.ncl .
ncl dcmip1-1_lat_lon_ng.ncl
ncl dcmip1-2_lat_height_ng.ncl

# print timing info
cat HommeTime_stats_DCMIP1-2 | grep walltotal
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1 | grep prim_run`"
echo "DCMIP1-2 `cat HommeTime_stats_DCMIP1-2 | grep prim_run`"
echo

# print error norms
cp $TEST1_DIR/dcmip1-1_error_norm_ng.ncl .
cp $TEST2_DIR/dcmip1-2_error_norm_ng.ncl .
ncl dcmip1-1_error_norm_ng.ncl | tail -n 1
ncl dcmip1-2_error_norm_ng.ncl | tail -n 1
echo

exit

