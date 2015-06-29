#!/bin/tcsh
#PBS -q debug
#PBS -A acme
#PBS -l walltime=00:20:00
#PBS -j oe
#PBS -o out_ne120_$PBS_JOBID
#PBS -e err_ne120_$PBS_JOBID

# Be sure to change MPI/threads here and NCPU/NTHREADS below
#XXX -l mppwidth=960
#PBS -l mppwidth=10800

#_______________________________________________________________________
#
# This file executes a set of DCMIP test to verify tracer perfrormance
#
#  1) Edit the path to the configuration script, if necessary
#  2) Submit this script to the queue or execute it an interactive session
#_______________________________________________________________________
set NCPU     = 10800      # number of CPUs to use
set NTHREADS = 1          # number of openMP threads per MPI task
set NCPU_PER_NODE = 24     # number of MPI tasks per node


set NE       = 120        # number of elements per cube-edge
set TSTEP    = 75         # timestep size
set NU       = 1e13       # hyperviscosity coefficient
set QSIZE1 = 50  
set CONFIGURE_DIR = ../../


# diagnostic output every 2day
set statefreq = 48     
@ statefreq *= 3600
@ statefreq /= $TSTEP


# Run length (hours) for DCMIP1-1
# 288h (12 days) needed for verification.  expensive!
# 6h is probably good for performance testing.  
#set nhours = 288
set nhours = 1

# convert to timesteps
set nmax = $nhours
@ nmax *= 3600
@ nmax /= $TSTEP

# edison has 2 sockets per node
set NUM_NUMA = $NCPU_PER_NODE
@ NUM_NUMA /= 2



#_______________________________________________________________________
# get path variables from configure script:

if ($?PBS_O_WORKDIR) then
  cd $PBS_O_WORKDIR
endif

cd ${CONFIGURE_DIR}; source configure.sh
# override RUN_DIR, so we can que multiple runs at once:
set RUN_DIR = ${RUN_DIR}_$$

# override default RUN_COMMAND set in configure.sh
set RUN_COMMAND="aprun -n $NCPU -N $NCPU_PER_NODE -d $NTHREADS -S $NUM_NUMA -ss -cc numa_node "
echo RUN_COMMAND:
echo $RUN_COMMAND





set TEST1_DIR = $BLD_DIR/test/dcmip1-1  # test case directory
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
sed s/output_frequency.\*/output_frequency=0/          |\
sed s/NThreads.\*/NThreads=$NTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/dcmip1-1_NE120.nl


#_______________________________________________________________________
# launch the executable

cd $RUN_DIR
setenv OMP_NUM_THREADS $NTHREADS
setenv OMP_STACKSIZE 64M   
date

# run dcmip test 1-1
echo "executing dcmip test 1-1"
echo "${RUN_COMMAND}  $EXE < dcmip1-1_NE120.nl"
${RUN_COMMAND}  $EXE < dcmip1-1_NE120.nl
if($status) exit
mv HommeTime_stats HommeTime_stats_DCMIP1-1_NE120
date


# print timing info
cat HommeTime_stats_DCMIP1-1_NE120 | grep walltotal
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE120 | grep prim_run`"
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE120 | grep prim_advance_exp`"
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE120 | grep prim_advec_tracers`"
echo "DCMIP1-1 `cat HommeTime_stats_DCMIP1-1_NE120 | grep vertical_remap`"
echo



exit





