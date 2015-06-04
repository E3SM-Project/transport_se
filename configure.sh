#!/bin/csh
#
#  To configure and build the transport_se tracer mini-app:
#
#  1) Copy this script to your build directory
#  2) Create a cmake machine file specific to your platform, if necessary
#  2) Edit the user-defined path variables and parameters below
#  3) Execute this script
#
#  note: this script is also used by test/run_tracer_tests.sh
#_______________________________________________________________________

echo "Please edit this file to specify user defined parameters"; exit # delete this line

# USER-DEFINED PARAMETERS:

module load cmake netcdf hdf5 ncl                 # load required modules

set REPO=$HOME/CODE/transport_se                  # transport_se repository
set WORK=$cwd                                     # work directory for building and running
set MACH=$REPO/cmake/machineFiles/edison.cmake    # machine specific cmake file

set QSIZE=4                                       # max number of tracers
set NLEV=64                                       # number of vertical levels

set RUN_COMMAND="aprun -n"                        # command to launch parallel executable
#set RUN_COMMAND="mpirun -np"

echo "REPO=$REPO"
echo "WORK=$WORK"

#_______________________________________________________________________
# check for some common errors
@ err=0

if ($cwd == $REPO) then
  echo "please copy this script to your build directory"
  @ err=$err + 1
endif

if (! -f $MACH  ) then
  echo "machine file not found: $MACH"
  @ err=$err + 1
endif

if(`where ncl` == "") then
  echo "ncl command not found"
  @ err=$err + 1
endif

if($err > 0) exit

#_______________________________________________________________________
# script derived parameters

set BLD_DIR   = $WORK/bld                         # build directory
set RUN_DIR   = $WORK/run                         # run directory
set EXE       = $BLD_DIR/src/preqx/preqx          # location of executable

#_______________________________________________________________________
# configure
#
#   for a debug configuration add these flags:
#     -DCMAKE_BUILD_TYPE=Debug \
#     -DCMAKE_Fortran_FLAGS_DEBUG="-g -O0 -fbounds-check" \
#
#   for native-grid (non-interpolated) output use:
#     -DPREQX_USE_PIO=TRUE \
#

mkdir -p $BLD_DIR
mkdir -p $RUN_DIR
cd $BLD_DIR

if ( -f $EXE  ) then
  echo "found $EXE executable. assuming configuration was already successful"

else
  echo "running CMAKE to configure the model"
  rm -rf CMakeFiles CMakeCache.txt src

  cmake -C $MACH                \
   -DQSIZE_D=$QSIZE             \
   -DPREQX_PLEV=$NLEV           \
   -DPREQX_NP=4                 \
   -DPREQX_USE_PIO=FALSE        \
   -DBUILD_HOMME_SWEQX=FALSE    \
   -DENABLE_OPENMP=TRUE         \
   $REPO

  make -j8 clean
endif

#_______________________________________________________________________
# build
#

echo "compiling the transport_se executable"
make -j8 preqx

exit

