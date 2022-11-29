#ifndef DATAMANAGER_H
#define DATAMANAGER_H 

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>

#include "TString.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TCanvas.h"

// Class to handle and manipulate data  
class DataManager
{
public:
    DataManager();
    ~DataManager();

    std::vector<float> DefinePhantom(double width, double height, double depth);
    std::vector<float> DefineVoxel(float width, float height, float depth);
    float              GetVoxVol(std::vector<float>& voxBox);

    std::vector<int> GetNumberOfVoxels(std::vector<float>& edges, std::vector<float>& vox);

    int GetNumberOfTimeSteps(float tofMax=100., float tStep=2.); 
    int GetBin(int nBins, float lowEdge, float highEdge, float val);

    std::vector<std::vector<double>> SumEBin(double* binData, int axis, std::vector<int>& numVoxels,   
                                            std::vector<float>& waterBox, std::vector<float>& voxelBox,
                                            float offset);
    
    void PlotGraph(int bins, TString pltTitle, double* hAxis, TString hAxisName, double* vAxis, TString vAxisName);
    void Plot3Coord(double* binData, std::vector<int>& numVoxels, std::vector<float>& waterBox, std::vector<float>& voxelBox);
  
};
#endif // DATAMANAGER_H
