#!/bin/bash
module load intel/2021.1.2
module load openmpi/intel-2021.1/4.1.0
module load fftw/intel-2021.1/openmpi-4.1.0/3.3.9
module load hdf5/intel-2021.1/openmpi-4.1.0/1.10.6
module load netcdf/intel-2021.1/hdf5-1.10.6/openmpi-4.1.0/4.7.4
export PSPLINE_HOME=/projects/TRANSP/opt/pspline-2.0.0/intel-2021.1.2
export NTCC_HOME=$HOME/ntcc
export I2MEX_HOME=$HOME/ntcc

export IMAS_PREFIX=$HOME/IMAS_3.34
export LD_LIBRARY_PATH=$IMAS_PREFIX/lib:$LD_LIBRARY_PATH
export IMAS_VERSION=3.34.0
export PYTHONPATH=$IMAS_PREFIX/python/lib

export TRANSP2FAR3D_HOME=$HOME/Transp2Far3d

make all
