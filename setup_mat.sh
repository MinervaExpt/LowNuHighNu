#!/bin/bash

export EXPERIMENT=minerva
export IFDH_DEBUG=0

source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
spack load root@6.28.12
spack load cmake
spack load gcc
spack load fife-utils

source /cvmfs/minerva.opensciencegrid.org/minerva/setup/setup_minerva_products.sh

export TOPDIR=/exp/minerva/app/users/$USER/MATAna/
WORKINGDIR=/exp/minerva/app/users/$USER/MATAna/LowNuHighNu
cd $WORKINGDIR

source ../opt/bin/setup.sh
source ../opt/buildGENIEXSecExtract/setup_GENIEXSecExtract.sh 

export IFDH_DEBUG=0

export PYTHONPATH=$PYTHONPATH:$WORKINGDIR/py_classes
export PLOTROOT=$WORKINGDIR/plotting
export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}
