#!/bin/bash

export EXPERIMENT=minerva
export IFDH_DEBUG=0

source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3
spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v3
spack load ifdhc-config@2.6.20%gcc@11.4.1 arch=linux-almalinux9-x86_64_v3
spack load gcc
#spack load fife-utils@3.7.4
spack load python@3.9.15
spack load py-numpy@1.24.3%gcc@12.2.0

htgettoken -a htvaultprod.fnal.gov -i minerva
export BEARER_TOKEN_FILE=/run/user/`id -u`/bt_u`id -u`

# The checkout this script lives in, and its parent as TOPDIR (loadLibs.C
# compiles $TOPDIR/LowNuHighNu/includes).  For the main checkout these are
# /exp/minerva/app/users/$USER/MATAna/{LowNuHighNu,} as before; a second
# checkout (e.g. MATAna-wk45/LowNuHighNu, with MATAna-wk45/opt -> ../MATAna/opt)
# gets its own TOPDIR, so it never compiles another checkout's includes.
WORKINGDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export TOPDIR="$(dirname "$WORKINGDIR")/"
cd "$WORKINGDIR"

source ../opt/bin/setup.sh
source ../opt/buildGENIEXSecExtract/setup_GENIEXSecExtract.sh 

export IFDH_DEBUG=0

export PYTHONPATH=$PYTHONPATH:$WORKINGDIR/py_classes
export PLOTROOT=$WORKINGDIR/plotting
export LD_LIBRARY_PATH=${ROOTSYS}/lib/root:${LD_LIBRARY_PATH}


# New authentication protocol following the "death of proxy", circa May 2025
htgettoken -a htvaultprod.fnal.gov -i minerva

# this will renew my proxy if needed
jobsub_q -G minerva --user finer
