#include "DataManager.h"

DataManager::DataManager()
{
    
}

DataManager::~DataManager()
{
    
}

std::vector<float> DataManager::DefinePhantom(double width, double height, double depth)
{
    //  Input: Width, Height, Depth of box in mm.
    //
    //   Return: Cartesian coordinates of extremes of box as distance
    //             from the centre:
    //             edges[2]: edges[0] -- distance of x edges from centre
    //                       edges[1] -- distance of y edges from centre
    //                       edges[2] -- distance of z edges from centre
    
    float halfWidth = width/2.;
    float halfHeight = height/2.;
    float halfDepth = depth/2.;
    
     std::vector<float> edges{halfWidth, halfHeight, halfDepth}; 
     
     std::cout << " DefineWaterBox: edges of water box at +/- x, y, z (mm): " << edges[0] << ", " << edges[1] << ", " << edges[2] << std::endl;
     return edges;
}

std::vector<float> DataManager::DefineVoxel(float width, float height, float depth)
{
    // Input: Width, Height, Depth of voxel in mm.
    //
    // Return: Numpy array with voxel size in mm
    //             Vox[2]: Vox[0] --  x length of voxel
    //                     Vox[1] --  y length of voxel
    //                     Vox[2] --  z length of voxel

    std::vector<float> vox{width, height, depth};
    
    std::cout << " DefineVoxel: voxel size (mm): " << vox[0] << ", " << vox[1] << ", " << vox[2] << std::endl;
    
     return vox;
}

std::vector<int> DataManager::GetNumberOfVoxels(std::vector<float>& edges, std::vector<float>& vox)
{
    // Input:  edges (3 element float vector) of box (as defined by 
    //            DefineWaterBox) and voxel size (3 element float vector,
    //            as defined by DefineVox).
    //    Return: Number of bins in x, y, and z directions (integers)
    
    int NBinX = 2.*edges[0]/vox[0];
    int NBinY = 2.*edges[1]/vox[1];
    int NBinZ = 2.*edges[2]/vox[2];

    std::vector<int> voxelNum;
    voxelNum.push_back(NBinX);
    voxelNum.push_back(NBinY);
    voxelNum.push_back(NBinZ);
    
    std::cout << " GetNumberOfVoxels: x, y, and z: " << voxelNum[0] << ", " << voxelNum[1] << ", " << voxelNum[2] << std::endl;

    return voxelNum;
}

int DataManager::GetNumberOfTimeSteps(float tofMax, float tStep)
{
    // Input:  
    // - tofMax: "Typical" time period (float) for travel of proton in
    //              water phantom until it comes to rest (ns).
    // - tStep  : Time step (time bin) (ns).
    //            as defined by defineVox).
    //    Return: Number of bins in t (integer)

    int nBinT = tofMax/tStep;
    std::cout << " GetNumberOfTimeSteps: " << nBinT << std::endl;

    return nBinT;
}

int DataManager::GetBin(int nBins, float lowEdge, float highEdge, float val)
{
    // Return bin number.  

    float binWidth = (highEdge - lowEdge) / static_cast<float>(nBins);
    int bin = static_cast<int>( (val - lowEdge) / binWidth);
    return bin;
}

std::vector<std::vector<double>> DataManager::SumEBin(double* binData, int axis, std::vector<int>& numVoxels,
                                                      std::vector<float>& waterBox, std::vector<float>& voxelBox)
{
    // Input:  
    // - binData  : Pointer to binned data (energySum)
    // - axis     : The axis being plotted against energy (0, 1, or 2)
    // 			corresponding to (x, y, or z)
    // - numVoxels: Vector containing number of voxels
    // - waterBox : Vector of the dimensions of phantom
    // - voxelBox : Vector of the dimensions of a voxel
    //   Return   : 2D vector of summed energy bins along axis 
    
    int nBinsX = numVoxels[0];
    int nBinsY = numVoxels[1];
    int nBinsZ = numVoxels[2];

    //std::vector<std::vector<double>> sumBin(numVoxels[axis], std::vector<double>(numVoxels[axis],0));
    std::vector<std::vector<double>> sumBin(2, std::vector<double>(numVoxels[axis],0));
    if(axis==0)
    {
        for(int width=0; width<nBinsX; width++)
        {
            for(int depth=0; depth<nBinsZ; depth++)
            {
                for(int height=0; height<nBinsY; height++)
                {
                    //int index = depth + nBinsZ*(height + nBinsY * width);
                    int index = depth + nBinsZ*(width + nBinsX * height);
                    sumBin[1][width] += binData[index];

                }
            }
            sumBin[0][width] = -waterBox[0] + voxelBox[0]*(0.5 + width); // center of bin
        }
    }
    else if(axis==1)
    {
        for(int height=0; height<nBinsY; height++)
        {
            for(int width=0; width<nBinsX; width++)
            {
                for(int depth=0; depth<nBinsZ; depth++)
                {
                    int index = depth + nBinsZ*(width + nBinsX * height);
                    sumBin[1][height] += binData[index];
                }
            }
            sumBin[0][height] = -waterBox[1] + voxelBox[1]*(0.5 + height); // center of bin
        }
    }
    else if(axis==2)
    {
        for(int depth=0; depth<nBinsZ; depth++)
        {
            for(int height=0; height<nBinsY; height++)
            {
                for(int width=0; width<nBinsX; width++)
                {
                    int index = depth + nBinsZ*(width + nBinsX * height);
                    sumBin[1][depth] += (binData[index]);
                }
            }
            sumBin[0][depth] = -waterBox[2] + voxelBox[2]*(0.5 + depth); // center of bin
        }
    }

    return sumBin;
}

void DataManager::PlotGraph(int bins, TString pltTitle, double* hAxis, TString hAxisName, double* vAxis, TString vAxisName)
{
    // Plots graph of given horizontal and vertical axis
    TGraph* gr = new TGraph(bins, hAxis, vAxis);
    gr->SetMarkerColor(4);
    gr->SetMarkerSize(0.4);
    gr->SetMarkerStyle(21);

    gr->SetTitle(pltTitle);
    gr->GetXaxis()->SetTitle(hAxisName);
    gr->GetXaxis()->SetTitleSize(.04);
    gr->GetYaxis()->SetTitle(vAxisName);
    gr->GetYaxis()->SetTitleSize(.04);

    gr->Draw("AP");
}

void DataManager::Plot3Coord(double* binData, std::vector<int>& numVoxels, std::vector<float>& waterBox, std::vector<float>& voxelBox)
{
    // Create canvas
    TCanvas *c = new TCanvas("","",1000,1000,1100,950);
    c->Divide(1,3);
    
    // Plotting z
    c->cd(1);
    std::vector<std::vector<double>> sumBin = SumEBin(binData, 2, numVoxels, waterBox, voxelBox);
    double* horBin = &sumBin[0][0];
    double* verBin = &sumBin[1][0];
    
    // Plotting z
    TString pltTitle = "Energy distribution vs binned z";
    TString horAxisName = "z [mm]";
    TString verAxisName = "Integrated Energy [MeV]";
    PlotGraph(numVoxels[2],pltTitle,horBin,horAxisName,verBin,verAxisName);
    
    // Plotting x
    c->cd(2);
    sumBin = SumEBin(binData, 0, numVoxels, waterBox, voxelBox);
    horBin = &sumBin[0][0];
    verBin = &sumBin[1][0];
    pltTitle = "Energy distribution vs binned x";
    horAxisName = "x [mm]";
    verAxisName = "Integrated Energy [MeV]";
    PlotGraph(numVoxels[0],pltTitle,horBin,horAxisName,verBin,verAxisName);
    
    // Plotting y
    c->cd(3);
    sumBin = SumEBin(binData, 1, numVoxels, waterBox, voxelBox);
    horBin = &sumBin[0][0];
    verBin = &sumBin[1][0];
    pltTitle = "Energy distribution vs binned y";
    horAxisName = "y [mm]";
    verAxisName = "Integrated Energy [MeV]";
    PlotGraph(numVoxels[1],pltTitle,horBin,horAxisName,verBin,verAxisName);

    c->Print("energyplot-protons.png");
}
