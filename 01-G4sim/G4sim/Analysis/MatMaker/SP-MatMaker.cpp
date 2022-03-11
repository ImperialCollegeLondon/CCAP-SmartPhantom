#include "ProcessRoot.h"
#include "KWaveInput.h"
#include <stdlib.h>

int main(int argc, char* argv[])
{
    std::cout << "\n|-- SP-MatMaker --|\n" << std::endl;

    const std::string inputName = argv[1];          // Input Root File Name
    const std::string outputName = argv[2];         // Output Matlab File Name
    const std::string kwaveName = argv[3];          // K-Wave File Name
    float voxelHeight = strtof(argv[4],NULL);       // Voxel Height
    float voxelWidth = strtof(argv[5],NULL);        // Voxel Width
    float voxelDepth = strtof(argv[6],NULL);        // Voxel Depth
    double analysisHeight = strtod(argv[7],NULL);   // Analysis Height
    double analysisWidth = strtod(argv[8],NULL);    // Analysis Width
    double analysisDepth = strtod(argv[9],NULL);    // Analysis Depth    
    
    if( analysisHeight < voxelHeight || analysisWidth < voxelWidth || analysisDepth < voxelDepth )
    {
        std::cout << "ERROR: Analysis volume must be larger than or equal to a single voxel volume." << std::endl;
        return 1;
    }
    
    // Careful of order here (width, height, depth), optimise later to be more consistent with ordering of axes
    std::vector<float> voxelSize{voxelWidth, voxelHeight, voxelDepth};
    std::vector<double> analysisVolume{analysisWidth, analysisHeight, analysisDepth};

    ///////////////////////////////////////////
    // Process ROOT file + Write Matlab File //
    ///////////////////////////////////////////
    ProcessRoot rootFile(inputName, outputName, voxelSize, analysisVolume);

    ////////////////////////
    // Write K-Wave Input //
    ////////////////////////
    KWaveInput kwaveInput(kwaveName);  

    // Define bin dimensions
    int binDim[3];
    binDim[0] = rootFile.GetBinHeight();
    binDim[1] = rootFile.GetBinWidth();
    binDim[2] = rootFile.GetBinDepth();

    kwaveInput.SetBinDim(binDim);
    
    double* energyData = rootFile.GetEnergySum();
    kwaveInput.SetEnergyData(energyData);

    // Apply smoothening by multiplying by window filter
    kwaveInput.SetGrueneisen(0.11);             // Grueneisen parameter for water
    kwaveInput.SetSmooth(true);                 // Whether or not to smoothen data
    kwaveInput.SetRestoreMagnitude(true);       // Whether to restore magnitude after smooth
    kwaveInput.SetRotWindow(true);              // Type of window function 
                                                // (K-Wave uses rotational window by default)
    
    // Grid Size (should be consistent with Matlab convention)
    unsigned long Nx = binDim[0];               // Height
    unsigned long Ny = binDim[1];               // Width
    unsigned long Nz = binDim[2];               // Depth
    float dx = 1e-3;                            // grid point spacing in the x direction [m]
    float dy = 1e-3;                            // grid point spacing in the y direction [m]
    float dz = 1e-3;                            // grid point spacing in the z direction [m]
    
    kwaveInput.SetNx(Nx);
    kwaveInput.SetNy(Ny);
    kwaveInput.SetNz(Nz);
    kwaveInput.SetDx(dx);
    kwaveInput.SetDy(dy);
    kwaveInput.SetDz(dz);
    
    // Thickness of absorbing boundary (PML) around the simulation
    kwaveInput.SetPMLSizeX(10);                 // Not sure why despite pml_size = 1 in matlab 
                                                // it still writes default (default = 10)
    kwaveInput.SetPMLSizeY(10);
    kwaveInput.SetPMLSizeZ(10);
    kwaveInput.SetPMLAlphaX(2);                 // Absorption within the perfectly matched 
                                                // layer in Nepers per grid point (default = 2)
    kwaveInput.SetPMLAlphaY(2);
    kwaveInput.SetPMLAlphaZ(2);
    
    // Define the properties of the propagation medium
    float mediumSoundSpeed = 1500;              // [m/s]
    kwaveInput.SetMediumSoundSpeed(mediumSoundSpeed);
    kwaveInput.SetCRef(1500);                   // reference sound speed used within the
                                                // k-space operator (phase correction term)
    
    // Set the time limit and time steps
    float tEnd = (float)Nx*dx/mediumSoundSpeed;  // [s]
    kwaveInput.SetTEnd(tEnd);
    kwaveInput.SetCFL(0.3); 
    
    // Sensor Mask
    // Following should give the same indices as with: sensor.mask(Nx/2, Ny/4:4:3*Ny/4, 3*Nz/4) = 1
    std::vector<unsigned long> sensorMask = kwaveInput.SensorMaskIndices(Nx, Nx/2, Nx/2, 0,
                                                                         Ny, Ny/4, 3*Ny/4, 4,
                                                                         Nz, 3*Nz/4, 3*Nz/4, 0);
    kwaveInput.SetSensorMask(sensorMask);
    
    kwaveInput.WriteFile();
    
    std::cout << "Finished Main" << std::endl;
}
