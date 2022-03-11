============================================================
Example Geant4 simulation to run the SmartPhantom simulation
============================================================

Environment variables set by running:
    source startup.bash

Then, build and run the sumulation:
  * Within the build directory type: 
        cmake ../G4sim
  * Then type: 
        make 
     * (there will be some warnings but if it goes to 100% it should be fine)

To run simulation: 
  * From build directry run:
        . SP-MatMaker.sh
    (This runs the Geant4 simulation, conversion to Matlab format, and produces 
    a K-Wave input file.)  
    
  * The bash script just contains all the commands in one place, parts 
    of the bash script can be run separately.
    
    * To just run the Geant4 simulation
          ./SmartPhantom (macrofile)
      (This produces waterhits.root. The macrofile by default is run1.mac, but can 
      be modified to adjust some aspects of the simulation. However, this is not 
      fully implemented nor tested) 
      
      Some useful basic commands implemented:
          /SP/IO/inputFile (filename)               -- To specify input coordinate file
          
          /SP/IO/outputFile (filename)              -- To specify output root file
          
          /SP/detector/enablePhantomWall (bool)     -- To add phantom walls with a 10 mm
                                                       thickness made of G4_PLEXIGLASS.
                                                       Also contains an entrance window of 
                                                       thickness consisting of G4_AIR.
                                                       
          /SP/detector/enableScifi (bool)           -- To add scintillating fibre   
                                                       detectors, adds four stations 
                                                       (each consisting of two planes) 
                                                       placed at locations suited
                                                       for 200 MeV protons
                                                       
    * To run the conversion from root to Matlab and K-Wave
          ./MatMaker/SP-MatMaker (rootfile) (matlab filename) (kwave filename) 
                                    (voxelHeight) (voxelWidth) (voxelDepth)
                                    (analysisHeight) (analysisWidth) (analysisDepth)

      This produces reads the root file and outputs the Matlab file and K-Wave input file. 
      * The first argument is the root file name (by default is waterhits.root)
      * The second argument is the name of the produced Matlab file.
      * The third argument is the name of the K-Wave input file (HDF5 format).
      * The fourth argument controls the size of the Voxel Height (in mm).
      * The fifth argument controls the size of the Voxel Width (in mm).
      * The sixth argument controls the size of the Voxel Depth (in mm).
      * The seventh argument controls the size of the analysis region height (i.e. phantom height = 300 mm)
      * The eighth argument controls the size of the analysis region width (i.e. phantom width = 300 mm)
      * The ninth argument controls the size of the analysis region depth (i.e. phantom depth = 300 mm)
      
    * To adjust the K-Wave input variables edit 'SP-MatMaker.cpp'
      * Implemented commands for K-Wave Input are generally similar to Matlab
      * But only basic functions are implemented so far    

    * To run K-Wave from c++
          ./kspaceFirstOrder-OMP -i (kwave filename) -o (kwave outputname) [optional parameters]
