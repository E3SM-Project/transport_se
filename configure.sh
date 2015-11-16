#!/bin/csh
#
#  To configure and build the transport_se tracer mini-app:
#
#  1) Copy this script to your build directory
#  2) Create a cmake machine file specific to your platform, if necessary
#  2) Edit the user-defined path variables and parameters below
#  3) Execute this script
#
#  note: this script is also called by scripts in the test/ dir
#_______________________________________________________________________

echo "Please edit this file to specify user defined parameters"; exit # delete this line

# set USER-DEFINED PARAMETERS:

# load required modules
module load cmake cray-hdf5-parallel cray-netcdf-hdf5parallel ncl

set REPO=$HOME/CODE/transport_se                  # set transport_se repository
set WORK=$cwd                                     # set work directory for building and running
set MACH=$REPO/cmake/machineFiles/edison.cmake    # set machine specific cmake file

# if you change these, remove the executable to force this script to reconfigure
set QSIZE_D = 35                                  # set max number of tracers (array size)
set NLEV    = 72                                  # set number of vertical levels

echo "REPO=$REPO"
echo "WORK=$WORK"

#_______________________________________________________________________
# set script derived parameters

set BLD_DIR   = $WORK/bld                         # build directory
set EXE       = $BLD_DIR/src/preqx/preqx          # location of executable

#_______________________________________________________________________
# check for common errors

if ($cwd == $REPO) then
  echo "please copy this script to your build directory"; exit
endif

if (! -f $MACH  ) then
  echo "machine file not found: $MACH"; exit
endif

if(`where ncl` == "") then
  echo "ncl command not found"; exit
endif

#_______________________________________________________________________
# configure transport_se

mkdir -p $BLD_DIR
cd $BLD_DIR

if ( -f $EXE  ) then
  echo "found $EXE executable. assuming configuration was already successful"

else
  echo "running CMAKE to configure the model"
  rm -rf CMakeFiles CMakeCache.txt src

  cmake -C $MACH                \
   -DQSIZE_D=$QSIZE_D           \
   -DPREQX_PLEV=$NLEV           \
   -DPREQX_NP=4                 \
   -DENABLE_OPENMP=TRUE         \
   -DENABLE_COLUMN_OPENMP=FALSE \
   $REPO

  make -j8 clean
endif

#_______________________________________________________________________
# build the executable

echo "compiling the transport_se executable"
make -j8 preqx

exit

