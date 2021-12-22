Example Geant4 simulation to run the SmartPhantom simulation
============================================================

Environment variables set by running:

source startup.bash

Then, build and run the sumulation:
  * Within the build directory type: "cmake ../G4sim"
  * Then type: "make" 
     - (there will be some warnings but if it goes to 100% it should be fine)

To run simulation: 
  * From build directry run:
       ./SmartPhantom run1.mac
    This produces waterhits.root.

  * To convert results to ascii execute: 
       root -l root2Ascii.cpp
    This produces waterhits.dat.
