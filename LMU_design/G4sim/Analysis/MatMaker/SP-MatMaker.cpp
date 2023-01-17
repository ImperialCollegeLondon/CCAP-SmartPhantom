#include "ProcessRoot.h"
#include "KWaveInput.h"
#include <stdlib.h>

int main(int argc, char* argv[])
{
    std::cout << "\n|-- SP-MatMaker --|\n" << std::endl;

    const std::string rootName = argv[1];             // Input Root File Name
    const std::string energyMatName = argv[2];        // Output Matlab File Name (Energy)
    const std::string energyDensityMatName = argv[3]; // Output Matlab File Name (Energy Density)
    const std::string kwaveName = argv[4];            // K-Wave File Name
    float voxelHeight = strtof(argv[5],NULL);         // Voxel Height
    float voxelWidth = strtof(argv[6],NULL);          // Voxel Width
    float voxelDepth = strtof(argv[7],NULL);          // Voxel Depth
    double analysisHeight = strtod(argv[8],NULL);     // Analysis Height
    double analysisWidth = strtod(argv[9],NULL);      // Analysis Width
    double analysisDepth = strtod(argv[10],NULL);     // Analysis Depth    
    
    if( analysisHeight < voxelHeight || analysisWidth < voxelWidth || analysisDepth < voxelDepth )
    {
        std::cout << "ERROR: Analysis volume must be larger than or equal to a single voxel volume." << std::endl;
        return 1;
    }
    
    std::vector<float> voxelSize{voxelWidth, voxelHeight, voxelDepth};
    std::vector<double> analysisVolume{analysisWidth, analysisHeight, analysisDepth};

    ///////////////////////////////////////////
    // Process ROOT file + Write Matlab File //
    ///////////////////////////////////////////
    ProcessRoot rootFile(rootName);
    rootFile.Initialise(voxelSize,analysisVolume);
    rootFile.ReadFile();
    rootFile.WriteMatFile(energyMatName, energyDensityMatName);
    rootFile.CloseFile();

    ////////////////////////
    // Write K-Wave Input //
    ////////////////////////
    KWaveInput kwaveInput(kwaveName);  

    // Define bin dimensions
    // **************************************************************
    int binDim[3];
    binDim[0] = rootFile.GetBinHeight();
    binDim[1] = rootFile.GetBinWidth();
    binDim[2] = rootFile.GetBinDepth();
    
    size_t binDimSize = sizeof(binDim)/sizeof(*binDim);
    
    kwaveInput.SetBinDim(binDim,binDimSize);
    
    double* energyDensityData = rootFile.GetEnergyDensity();   // Energy Density [J/m^3]
    kwaveInput.SetEnergyDensity(energyDensityData);

    // Apply smoothening by multiplying by window filter
    // **************************************************************
    kwaveInput.SetGrueneisen(0.11);                 // Grueneisen parameter for water
    kwaveInput.SetSmooth(true);                     // Whether or not to smoothen data
    kwaveInput.SetRestoreMagnitude(true);           // Whether to restore magnitude after smooth
    kwaveInput.SetRotWindow(true);                  // Type of window function 
                                                    // (K-Wave uses rotational window by default)
    
    // Grid Size (should be consistent with Matlab convention)
    // **************************************************************
    unsigned long Nx = binDim[0];                   // Height
    unsigned long Ny = binDim[1];                   // Width
    unsigned long Nz = binDim[2];                   // Depth
    float dx = voxelHeight*1.0e-3;                  // grid point spacing in the x direction [m]
    float dy = voxelWidth*1.0e-3;                   // grid point spacing in the y direction [m]
    float dz = voxelDepth*1.0e-3;                   // grid point spacing in the z direction [m]
    
    kwaveInput.SetNx(Nx);
    kwaveInput.SetNy(Ny);
    kwaveInput.SetNz(Nz);
    kwaveInput.SetDx(dx);
    kwaveInput.SetDy(dy);
    kwaveInput.SetDz(dz);
    
    // Thickness of absorbing boundary (PML) around the simulation
    // **************************************************************
    kwaveInput.SetPMLSizeX(10);                 // Not sure why despite pml_size = 1 in matlab 
                                                // it still writes default (default = 10)
    kwaveInput.SetPMLSizeY(10);
    kwaveInput.SetPMLSizeZ(10);
    kwaveInput.SetPMLAlphaX(2);                 // Absorption within the perfectly matched 
                                                // layer in Nepers per grid point (default = 2)
    kwaveInput.SetPMLAlphaY(2);
    kwaveInput.SetPMLAlphaZ(2);
    
    // Define the properties of the propagation medium
    // **************************************************************
    float mediumSoundSpeed = 1500;                      // [m/s]
    kwaveInput.SetMediumSoundSpeed(mediumSoundSpeed);
    kwaveInput.SetCRef(mediumSoundSpeed);               // reference sound speed used within the
                                                        // k-space operator (phase correction term)
    
    // Set the time limit and time steps
    // **************************************************************
    float tEnd = (float)Nx*dx/mediumSoundSpeed;  // [s]
    float cfl = 0.3;
    kwaveInput.SetTEnd(tEnd);
    kwaveInput.SetCFL(cfl); 
    
    // Sensor Mask
    // **************************************************************
    // Following should give the same indices as with: sensor.mask(Nx/2, Ny/4:4:3*Ny/4, 3*Nz/4) = 1
    float senMinX = (float)Nx/2.;       // Minimum for x
    float senMaxX = (float)Nx/2.;       // Maximum for x
    int   senDx   = 0;                  // Stepsize for x
    
    float senMinY = (float)Ny/4.;       // Minimum for y
    float senMaxY = (float)3*Ny/4.;     // Maximum for y
    int   senDy   = 4;                  // Stepsize for y
    
    float senMinZ = (float)3*Nz/4.;     // Minimum for z
    float senMaxZ = (float)3*Nz/4.;     // Maximum for z
    int   senDz   = 0;                  // Stepsize for z
    
    std::vector<unsigned long> sensorMask = kwaveInput.SensorMaskIndices(Nx, senMinX, senMaxX, senDx,
                                                                         Ny, senMinY, senMaxY, senDy,
                                                                         Nz, senMinZ, senMaxZ, senDz);
    kwaveInput.SetSensorMask(sensorMask);
    
    // Write HDF5 File
    // **************************************************************    
    kwaveInput.WriteFile();
    
    std::cout << "Finished Main" << std::endl;
}
