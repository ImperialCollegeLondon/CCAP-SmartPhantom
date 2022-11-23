#!/bin/bash
ulimit -s unlimited
 
#####################
# Run G4 Simulation #
#####################
./SmartPhantom run1.mac

####################
# Compile MatMaker #
####################
cd MatMaker
h5c++ -o SP-MatMakerAS SP-MatMakerAS.cpp ProcessRoot.cpp DataManager.cpp KWaveInput.cpp MatFile/MatFile.cpp Polyfit/PolyFit.cpp -L/vols/ccap/users/htl17/FFTW/lib -L/usr/lib64 -lfftw3 -lm `root-config --cflags --glibs`
cd ../

################
# Run MatMaker #
################
# Arguments are: Input File, Matlab Filename (energy), Matlab Filename(energy density), k-Wave Input Filename, 
#                Voxel-Radius, Voxel-Depth(axis z in Geant4),
#                Analysis-Radius, Analysis-Depth(axis z in Geant4)

./MatMaker/SP-MatMakerAS "waterhits.root" "energy_data2_cpp_AS.mat" "energy_density_data2_cpp_AS.mat" "kwave_input_cpp_AS.h5" 0.5 5 30 300

##############
# Run K-Wave #
##############
./kspaceFirstOrder-OMP -i "kwave_input_cpp_AS.h5" -o "kwave_output_cpp_AS.h5" --p_final
