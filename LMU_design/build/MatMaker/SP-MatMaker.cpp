#include "ProcessRoot.h"
#include <stdlib.h>

int main(int argc, char* argv[])
{
    std::cout << "\n|-- SP-MatMaker --|\n" << std::endl;

    const std::string rootName = argv[1];             // Input Root File Name
    const std::string energyMatName = argv[2];        // Output Matlab File Name (Energy)
    float voxelHeight = strtof(argv[3],NULL);         // Voxel Height
    float voxelWidth = strtof(argv[4],NULL);          // Voxel Width
    float voxelDepth = strtof(argv[5],NULL);          // Voxel Depth
    double analysisHeight = strtod(argv[6],NULL);     // Analysis Height
    double analysisWidth = strtod(argv[7],NULL);      // Analysis Width
    double analysisDepth = strtod(argv[8],NULL);      // Analysis Depth    
    
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
    rootFile.WriteMatFile(energyMatName);
    rootFile.CloseFile();
    
    std::cout << "Finished Main" << std::endl;
}
