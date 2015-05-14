#!/bin/tcsh
#PBS -l walltime=0:10:00
#PBS -j oe
#PBS -l nodes=35:ppn=4

# This file will configure, build, and run test DCMIP 1-1
#
# 1) Copy this file to your work directory
# 2) Create a cmake machine file for your platform if necessary
# 3) Edit the user-defined configuration parameters below
# 4) Execute this script from the work directory. (It will configure and build the code.)

#_______________________________________________________________________
# USER-DEFINED PARAMETERS

echo "Please edit this file to specify user defined parameters";exit

set REPO=$HOME/CODE/transport_se                  # transport_se code repository
set WORK=$GSCRATCH/transport_se                   # work directory, for building and running
set MACH=$REPO/cmake/machineFiles/edison.cmake    # machine specific CMAKE directory
set QSIZE=4                                       # max number of tracers
set NLEV=64                                       # number of vertical levels

set NE=8                                          # number of elements per cube-edge
set NCPU=2                                        # number of CPUs to use
set NTHREADS=1                                    # number of openMP threads

module load cmake netcdf hdf5

#_______________________________________________________________________
# script derived parameters

set BLD       = $WORK/bld               # build directory
set RUN_DIR   = $WORK/run               # run directory
set EXE       = $BLD/src/preqx/preqx    # location of executable
set TEST_DIR  = $REPO/test/             # test case directory
set VCOORD    = $TEST_DIR/vcoord        # location of vertical coordinate files

#_______________________________________________________________________
# configure the build

if ( -f $EXE  ) then
  echo "Found executable.  Assuming configuration was already successful"

else
  echo "running CMAKE to configure the model"
  mkdir -p $BLD
  cd $BLD
  rm -rf CMakeFiles CMakeCache.txt src

  cmake -C $MACH                \
   -DQSIZE_D=$QSIZE             \
   -DPREQX_PLEV=$NLEV           \
   -DPREQX_NP=4                 \
   -DPREQX_USE_PIO=TRUE         \
   -DBUILD_HOMME_SWDGX=FALSE    \
   -DBUILD_HOMME_SWEQX=FALSE    \
   -DBUILD_HOMME_PRIMDGX=FALSE  \
   $REPO

  make -j4 clean
endif

#_______________________________________________________________________
# compile the code

echo "Compiling the transport_se executable"
make -j4 preqx

#_______________________________________________________________________
# create directory for simulation output

mkdir -p $RUN_DIR/movies
cd $RUN_DIR
cp -a $VCOORD vcoord

#_______________________________________________________________________
# set time-step and hyperviscosity from resolution

set TSTEP
if ( $NE == 8 ) then
# 3.75 degree
set TSTEP  = 1440
set NU     = 2e16
endif

if ( $NE == 30 ) then
# 1 degree
set TSTEP  = 300
set NU     = 1e15
endif

if ( $NE == 120 ) then
# 0.25 degree. ACME target resolution.
set TSTEP  = 75
set NU     = 1e13
endif

echo "Custom run paramteres:"
echo "ne    = $NE"
echo "tstep = $TSTEP"
echo "nu_q  = $NU"
echo

#_______________________________________________________________________
# create namelist from user defined parameters

cd $TEST_DIR
sed s/NE.\*/$NE/ dcmip1-1.nl          |\
sed s/qsize.\*/qsize=$QSIZE/          |\
sed s/NThreads.\*/NThreads=$NTHREADS/ |\
sed s/nu_q.\*/nu_q=$NU/  >  $RUN_DIR/input.nl

#_______________________________________________________________________
# launch the executable

cd $RUN_DIR
setenv OMP_NUM_THREADS $NTHREADS
date

# runjob -np $NCPU $EXE < input.nl
# mpirun -np $NCPU $EXE < input.nl
aprun  -n  $NCPU $EXE < input.nl

date
exit




