=================================================================================================
=====            Example Geant4 simulation to run the SmartPhantom simulation               ===== 
=================================================================================================

=================================================================================================

Environment variables set by running:
    source startup.bash

Then, build and run the sumulation:
  * Within the build directory type: 
        cmake ../G4sim

  * Then type: 
        make
           
    * (there will be a bunch of warnings but if it goes to 100% it should be fine)
     
If not running on lx machines, you will need the following installed and adjust the compile 
arguments to direct to the correct directories:
        ROOT
        Geant4 
        HDF5 library 
        FFTW library

=================================================================================================

To run simulation: 
  * From build directry run:
        . SP-MatMaker.sh                or                . SP-MatMakerAS.sh      
        
  * (This runs the Geant4 simulation, conversion to Matlab format, and produces 
    a K-Wave input file.)  
    
  * The bash script just contains all the commands in one place, parts 
    of the bash script can be run separately edit the arguments in script as specified below.
    
  * The two bash scripts:
                SP-MatMaker.sh                      -- For 3D cartesian coordinates
                SP-MatMakerAS.sh                    -- For axisymmetric(cylindrical) coordinates
                
=================================================================================================
                
To just run the Geant4 simulation
    ./SmartPhantom (macrofile)

  * (This produces waterhits.root by default. The macrofile by default is run1.mac, but can 
    be modified to adjust some aspects of the simulation. However, this is not fully implemented) 
      
    * Some useful basic commands implemented:
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

=================================================================================================

To run the conversion from ROOT to Matlab (an energy and energy density file) and k-Wave for 3D cartesian
    ./MatMaker/SP-MatMaker (rootfile), (matlab energy filename), (matlab energy density filename), (kwave filename), 
                           (voxelHeight), (voxelWidth), (voxelDepth),
                           (analysisHeight), (analysisWidth), (analysisDepth)

  * This produces reads the root file and outputs the Matlab file and k-Wave input file. 
      * The 1st argument is the root file name (by default is waterhits.root)
      * The 2nd argument is the name of the produced Matlab (energy) file.
      * The 3rd argument is the name of the produced Matlab (energy density) file.
      * The 4th argument is the name of the K-Wave input file (HDF5 format).
      * The 5th argument controls the size of the voxel height (in mm).
      * The 6th argument controls the size of the voxel width (in mm).
      * The 7th argument controls the size of the voxel depth (in mm).
      * The 8th argument controls the size of the analysis region height (i.e. phantom height = 300 mm)
      * The 9th argument controls the size of the analysis region width (i.e. phantom width = 300 mm)
      * The 10th argument controls the size of the analysis region depth (i.e. phantom depth = 300 mm)

Alternatively, to run the conversion from ROOT to Matlab and k-Wave for axisymmetric
    ./MatMaker/SP-MatMakerAS (rootfile), (matlab energy filename), (matlab energy density filename), (kwave filename), 
                             (voxelRadius), (voxelDepth),
                             (analysisRadius), (analysisDepth)

  * This produces reads the root file and outputs the Matlab file and k-Wave input file. 
      * The 1st argument is the root file name (by default is waterhits.root)
      * The 2nd argument is the name of the produced Matlab (energy) file.
      * The 3rd argument is the name of the produced Matlab (energy density) file.
      * The 4th argument is the name of the K-Wave input file (HDF5 format).
      * The 5th argument controls the size of the voxel radius (in mm).
      * The 6th argument controls the size of the voxel depth (in mm).
      * The 7th argument controls the size of the analysis radius (in mm).
      * The 8th argument controls the size of the analysis depth (i.e. phantom depth = 300 mm)

=================================================================================================
      
To adjust the K-Wave input variables edit 'SP-MatMaker.cpp'
    * Implemented commands for K-Wave Input are generally similar to Matlab
    * But only basic functions are implemented so far    
    
    * To run K-Wave from C++
          ./kspaceFirstOrder-OMP -i (kwave filename) -o (kwave outputname) [optional parameters]

=================================================================================================
