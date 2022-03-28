#include "ProcessRoot.h"
#include "KWaveInput.h"
#include <stdlib.h>

int main(int argc, char* argv[])
{
    std::cout << "\n|-- SP-MatMakerAS --|\n" << std::endl;

    const std::string inputName = argv[1];            // Input Root File Name
    const std::string energyMatName = argv[2];        // Output Matlab File Name (Energy)
    const std::string energyDensityMatName = argv[3]; // Output Matlab File Name (Energy Density)
    const std::string kwaveName = argv[4];            // K-Wave File Name
    float voxelRadius = strtof(argv[5],NULL);         // Voxel Radius
    float voxelDepth = strtof(argv[6],NULL);          // Voxel Depth
    double analysisRadius = strtod(argv[7],NULL);     // Analysis Radius
    double analysisDepth = strtod(argv[8],NULL);      // Analysis Depth
    
    if( analysisRadius < voxelRadius || analysisDepth < voxelDepth )
    {
        std::cout << "ERROR: Analysis volume must be larger than or equal to a single voxel volume." << std::endl;
        return 1;
    }
    
    std::vector<float> voxelSize{voxelRadius, voxelDepth};
    std::vector<double> analysisVolume{analysisRadius, analysisDepth};

    ///////////////////////////////////////////
    // Process ROOT file + Write Matlab File //
    ///////////////////////////////////////////
    ProcessRoot rootFile(inputName);
    rootFile.InitialiseAS(voxelSize,analysisVolume);
    rootFile.ReadFileAS();
    rootFile.WriteMatFileAS(energyMatName, energyDensityMatName);
    rootFile.CloseFile();
    
    ////////////////////////
    // Write K-Wave Input //
    ////////////////////////
    KWaveInput kwaveInput(kwaveName);  

    // Define bin dimensions
    // **************************************************************
    int binDim[2];
    binDim[0] = rootFile.GetBinDepth();
    binDim[1] = rootFile.GetBinRadius();

    size_t binDimSize = sizeof(binDim)/sizeof(*binDim);

    kwaveInput.SetBinDim(binDim,binDimSize);
    
    double* energyDensityData = rootFile.GetEnergyDensity();   // Energy Density [J/m^2]
    kwaveInput.SetEnergyDensity(energyDensityData);

    // Apply smoothening by multiplying by window filter
    // **************************************************************
    kwaveInput.SetGrueneisen(0.11);                 // Grueneisen parameter for water
    kwaveInput.SetSmooth(true);                     // Whether or not to smoothen data
    kwaveInput.SetRestoreMagnitude(true);           // Whether to restore magnitude after smooth
    kwaveInput.SetRotWindow(true);                  // Type of window function (some discrepancy if false) 
    
    // Grid Size (should be consistent with Matlab/k-Wave convention)
    // **************************************************************
    unsigned long Nx = binDim[0];                   // Depth
    unsigned long Ny = binDim[1];                   // Radius
    unsigned long Nz = 1;               
    float dx = voxelDepth*1.0e-3;                   // grid point spacing in the x direction [m]
    float dy = voxelRadius*1.0e-3;                  // grid point spacing in the y direction [m]
    
    kwaveInput.SetNx(Nx);
    kwaveInput.SetNy(Ny);
    kwaveInput.SetNz(Nz);
    kwaveInput.SetDx(dx);
    kwaveInput.SetDy(dy);
    
    // Thickness of absorbing boundary (PML) around the simulation
    // **************************************************************
    kwaveInput.SetPMLSizeX(20);
    kwaveInput.SetPMLSizeY(20);
    kwaveInput.SetPMLSizeZ(0);
    kwaveInput.SetPMLAlphaX(2);
    kwaveInput.SetPMLAlphaY(2);
    
    // Define the properties of the propagation medium
    // **************************************************************
    float mediumSoundSpeed = 1500.;                     // [m/s]
    kwaveInput.SetMediumSoundSpeed(mediumSoundSpeed);
    kwaveInput.SetCRef(mediumSoundSpeed);               // reference sound speed used within the
                                                        // k-space operator (phase correction term)
    
    // Set the time limit and time steps
    // **************************************************************
    float tEnd = (float)Ny*dy/mediumSoundSpeed;  // [s]
    float cfl = 0.3;
    kwaveInput.SetTEnd(tEnd);
    kwaveInput.SetCFL(cfl); 
    
    // Sensor Mask
    // **************************************************************
    // Following should give the same indices as with: sensor.mask(Nx, Ny/4:4:3*Ny/4) = 1
    float senMinX = (float)Nx;          // Minimum for x
    float senMaxX = (float)Nx;          // Maximum for x
    int   senDx   = 0;                  // Stepsize for x
    
    float senMinY = (float)Ny/4.;       // Minimum for y
    float senMaxY = (float)3*Ny/4.;     // Maximum for y
    int   senDy   = 4;                  // Stepsize for y
    
    std::vector<unsigned long> sensorMask = kwaveInput.SensorMaskIndices(Nx, senMinX, senMaxX, senDx,
                                                                         Ny, senMinY, senMaxY, senDy);
    
    kwaveInput.SetSensorMask(sensorMask);
    
    // Write HDF5 File
    // **************************************************************
    kwaveInput.WriteFileAS();
    
    std::cout << "Finished Main" << std::endl;
}
