#!/bin/bash

voms-proxy-destroy
kx509
voms-proxy-init -rfc --voms=fermilab:/fermilab/minerva/Role=Analysis --noregen -valid 24:0

export EXPERIMENT=minerva
export IFDH_DEBUG=0

source /cvmfs/minerva.opensciencegrid.org/minerva/hep_hpc_products/setups
setup root v6_10_04d -q e14:prof
setup cmake v3_7_1

source /cvmfs/minerva.opensciencegrid.org/minerva/setup/setup_minerva_products.sh

export TOPDIR=/minerva/app/users/norrick/MATAna/
CMTUSER=/minerva/app/users/norrick/cmtuser
W=/minerva/app/users/norrick/MATAna/LowNuHighNu
d=/pnfs/minerva/scratch/users/norrick
D=/pnfs/minerva/persistent/users/norrick/LowNuHighNu/2023-06-05
b=/minerva/data/users/norrick/
m=/minerva/app/users/norrick/MATAna/MAT/
mm=/minerva/app/users/norrick/MATAna/MAT-MINERvA/

cd $W

source ../opt/bin/setupROOT6OnGPVMs.sh
source ../opt/bin/setup.sh
source ../opt/buildGENIEXSecExtract/setup_GENIEXSecExtract.sh 

export IFDH_DEBUG=0
