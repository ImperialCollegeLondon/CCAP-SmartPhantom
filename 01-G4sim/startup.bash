#!/bin/bash

[ -d "/cvmfs/sft.cern.ch/lcg/views/" ] &&                        \
    source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_96     \
	   x86_64-centos7-gcc62-opt

source /Users/kennethlong/cernbox/G4/setup.bash

source $g4PATH-install/bin/geant4.sh

export GEANT4_DIR=$g4PATH


