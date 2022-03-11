
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
h5c++ -o SP-MatMaker SP-MatMaker.cpp ProcessRoot.cpp DataManager.cpp KWaveInput.cpp MatFile/MatFile.cpp -L/vols/ccap/users/htl17/FFTW/lib -L/usr/lib64 -lfftw3 -lm `root-config --cflags --glibs`
cd ../

################
# Run MatMaker #
################
# Arguments are: Input File, Matlab Filename to create, k-Wave Input Filename to create, 
#                   Voxel-Height(axis y in Geant4), Voxel-Width(axis x in Geant4), Voxel-Depth(axis z in Geant4),
#                   Analysis-Height(axis y in Geant4), Analysis-Width(axis x in Geant4), Analysis-Depth(axis z in Geant4)

./MatMaker/SP-MatMaker "waterhits.root" "energy_density_data2_cpp.mat" "kwave_input_cpp.h5" 10 10 5 300 300 300

##############
# Run K-Wave #
##############
./kspaceFirstOrder-OMP -i "kwave_input_cpp.h5" -o "kwave_output_cpp.h5" --p_final
