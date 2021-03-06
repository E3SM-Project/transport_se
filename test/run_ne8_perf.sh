#!/bin/tcsh
#SBATCH -p debug
#SBATCH -t 00:05:00
#SBATCH -N 12
#PBS -q debug
#PBS -l walltime=0:30:00
#PBS -l mppwidth=384
#PBS -j oe


# Runs DCMIP 1-1 with NE=8 and analyzes run times.
#_______________________________________________________________________
# call configure script to ensure executable exists
if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

set CONFIGURE_DIR = ../..     # location of configure script
cd ${CONFIGURE_DIR}; source configure.sh

#_______________________________________________________________________
# set test parameters

set TEST_NAME = run_ne8_perf  # name of test for run directory
set NE        = 8             # number of elements per cube-edge
set TSTEP     = 1200          # time step size, in second
set NU        = 6e16          # hyperviscosity coefficient
set QSIZE     = 35            # number of tracers
@ statefreq   = 480 * 3600 / $TSTEP             # set diagnostic display frequency


#_______________________________________________________________________
# compute run parameters from number of procs and number of threads
# 384 elements
set HTHREADS  = 4             # number of horizontal threads
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

# Edison
#set RUN_COMMAND = "aprun -n $NMPI -N $NMPI_PER_NODE -d $NTHREADS -S $NUM_NUMA -ss -cc numa_node"
# Cori
set RUN_COMMAND = "srun -n $NMPI -c $NTHREADS"
# Mira
#set RUN_COMMAND="qsub -t 15 -n $NNODES --proccount $NMPI --mode c$NMPI_PER_NODE "

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

# for Intel
if ( ${%OMP_STATUS} != 0  && ${%COLUMN_OMP_STATUS} != 0  ) then
  setenv OMP_NESTED 1 
  setenv KMP_HOT_TEAMS 1 
  setenv KMP_HOT_TEAMS_MAX_LEVEL 2 
  setenv OMP_PROC_BIND spread,close 
  setenv | grep KMP
endif


if( $NTHREADS > 1 && ${%OMP_STATUS} == 0 ) then
  echo "Error: NTHREADS > 1 requires ENABLE_OPENMP=TRUE"; exit
endif

if( $VTHREADS > 1 && ${%COLUMN_OMP_STATUS} == 0 ) then
  echo "Error: VTHREADS > 1 requires ENABLE_COLUMN_OPENMP=TRUE"; exit
endif

#_______________________________________________________________________
# create directories for simulation output
#
set JDATE=`date "+%y-%m-%d_%H%M%S"`
set PCONF="${NNODES}x${NMPI_PER_NODE}x${HTHREADS}x${VTHREADS}"
set RUN_DIR = $WORK/${TEST_NAME}_${PCONF}_$JDATE
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
# launch the executable

cd $RUN_DIR
date

# run dcmip test 1-1
#${RUN_COMMAND} -i dcmip1-1.nl $EXE # on Mira
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


# NCL on cori doesn't like threads
setenv OMP_NUM_THREADS 1

# print error norms
cp $TEST1_DIR/dcmip1-1_error_norm_ng.ncl .
ncl dcmip1-1_error_norm_ng.ncl | tail -n 1
echo

# produce plot of results
cp $TEST1_DIR/dcmip1-1_lat_lon_ng.ncl .
ncl dcmip1-1_lat_lon_ng.ncl

exit

