#define _USE_MATH_DEFINES
#include "DataManager.h"

DataManager::DataManager()
{}

DataManager::~DataManager()
{}

std::vector<float> DataManager::DefinePhantom(double width, double height, double depth)
{
    /*
        Define phantom dimensions

        Input:      width  -- Width [mm]
                    height -- Height [mm]
                    depth  -- Depth [mm]
                    
        Return:     Vector of cartesian coordinates of extremes of box as distances
                    from centre:
                    
                    edges[2]: edges[0] -- distance of x edges from centre
                              edges[1] -- distance of y edges from centre
                              edges[2] -- distance of z edges from centre
    */

    float halfWidth = width/2.;
    float halfHeight = height/2.;
    float halfDepth = depth/2.;
    
     std::vector<float> edges{halfWidth, halfHeight, halfDepth}; 
     
     std::cout << " DefineWaterBox: edges of water box at +/- x, y, z (mm): " << edges[0] << ", " << edges[1] << ", " << edges[2] << std::endl;
     return edges;
}

std::vector<float> DataManager::DefineCylinder(float radius, float length)
{
    /*
        Define cylinder dimensions

        Input:      radius  -- Radius [mm]
                    length  -- Length [mm]
                    
        Return:     Edges of cylinder from centre of water phantom:
                    
                    DimCyl[1]: DimCyl[0] -- distance to extreme Radius
                               DimCyl[1] -- half the maximum phi
                               DimCyl[2] -- distance of z edges from centre
    */
    
    std::vector<float> dimCyl{radius, 360./2, length/2};
    std::cout << " DefineCylinder: extremes of cylinder at R=" << dimCyl[0] << " mm, z = +/-" << dimCyl[1] << " mm" << std::endl;
    return dimCyl;
}

std::vector<float> DataManager::DefineVoxel(float width, float height, float depth)
{
    /*
        Define voxel dimensions

        Input:      width  -- Width [mm]
                    height -- Height [mm]
                    depth  -- Depth [mm]
                    
        Return:     Vector of voxel size in mm:
                    
                    Vox[2]: Vox[0] --  x length of voxel
                            Vox[1] --  y length of voxel
                            Vox[2] --  z length of voxel
    */
    
    std::vector<float> vox{width, height, depth};
    
    std::cout << " DefineVoxel: voxel size (mm): " << vox[0] << ", " << vox[1] << ", " << vox[2] << std::endl;
    
     return vox;
}

float DataManager::GetVoxVol(std::vector<float>& voxBox)
{
    /*
        Calculate voxel volume 
        
        Input:      voxBox -- Voxel Box coming from DefineVoxel
        
        Output:     Volume of voxel [mm^3]
    */
    
    float voxVol = voxBox[0]*voxBox[1]*voxBox[2];
    return voxVol;
}

std::vector<float> DataManager::DefineCylVoxel(float dR, float dL, float radius)
{
    /*
        Define cylindrical voxel dimensions

        Input:      dR      -- Radius Voxel [mm]
                    dL      -- Length Voxel [mm]
                    radius  -- Maximum Radius [mm]
                    
        Return:     Vector of cylindrical voxel:
                    
                    Vox[2]: Vox[0] --  r length of voxel
                            Vox[1] --  phi length of voxel (degrees)
                            Vox[2] --  z length of voxel
    */

    float dPhi0 = dR/radius;
    int   nBin0 = static_cast<int>(2.*M_PI / dPhi0);
    float dPhi1 = 2.*M_PI / nBin0;
    
    float     c = 2.*M_PI * radius;
    float   dR1 = c / static_cast<float>(nBin0);
    
          dPhi1 = dR1/radius;
    
    float phiVox = static_cast<float>(dPhi1*180./M_PI);
          
    std::vector<float> vox{dR, phiVox, dL};
    std::cout << " DefineCylVoxel: " << vox[0] << ", " << vox[1] << ", " << vox[2] << std::endl;

    return vox;
}

float DataManager::GetCylArea(std::vector<float>& cylVox)
{
    /*
        Define cylindrical voxel area (Double check this!!)
        
        Input:      cylVox      -- Cylinder Voxel from DefineCylVoxel
                    
        Return:     Area element of cylinder voxel
    */
    
    float dA = cylVox[0] * cylVox[2]; // (Is this right? dA ?= dr*dz)
    return dA;
}

float DataManager::GetCylVol(int rBin, std::vector<float>& cylVox)
{
    /*
        Define cylindrical voxel volume (Double check this!!)
        
        Input:      rBin        -- Radius Bin
                    cylVox      -- Cylinder Voxel from DefineCylVoxel
                    
        Return:     Volume element of cylinder voxel
    */
    
    float r = -cylVox[0]/2. + (rBin+1) * cylVox[0];
    float dV = r * cylVox[0] * cylVox[1]*M_PI/180. * cylVox[2]; // (Is this right? dV ?= r*dr*dphi*dz)
    return dV;
}

std::vector<int> DataManager::GetNumberOfVoxels(std::vector<float>& edges, std::vector<float>& vox)
{
    /*
        Calculate the number of voxels for each edge 

        Input:      edges      -- Edges of box (as defined by DefineWaterBox)
                    vox        -- Voxel size (as defined by DefineVoxel)
                    
        Return:     Number of bins in x, y, and z directions (integers)                    
    */
    
    int nBinX = 2.*edges[0]/vox[0];
    int nBinY = 2.*edges[1]/vox[1];
    int nBinZ = 2.*edges[2]/vox[2];

    std::vector<int> voxelNum{nBinX, nBinY, nBinZ};
    
    std::cout << " GetNumberOfVoxels: x, y, and z: " << voxelNum[0] << ", " << voxelNum[1] << ", " << voxelNum[2] << std::endl;

    return voxelNum;
}

std::vector<int> DataManager::GetNumberOfVoxelsCyl(std::vector<float> &dimCyl, std::vector<float> &cylVox)
{
    /*
        Calculate the number of voxels for cylinder in each dimension

        Input:      dimCyl      -- Cylinder dimensions (as defined by DefineCylinder)
                    cylVox      -- Voxel size of cylinder (as defined by DefineCylVoxel)
                    
        Return:     Number of bins in R, Phi, and Z directions(integers)                    
    */

    int   nBinR = std::round(dimCyl[0] / cylVox[0]);
    int nBinPhi = std::round(360. / cylVox[1]);
    int   nBinZ = std::round(2.*dimCyl[2]/cylVox[2]);
    
    std::vector<int> vox{nBinR, nBinPhi, nBinZ};
    std::cout << " GetNumberOfVoxelsCyl: r, phi, and z: " << vox[0] << ", " << vox[1] << ", " << vox[2] << std::endl;
    
    return vox;
}

int DataManager::GetNumberOfTimeSteps(float tofMax, float tStep)
{
    /*
        Calculate the time step

        Input:      tofMax      -- "Typical" time period (float) for travel of proton in
                                   water phantom until it comes to rest (ns).
                    tStep       -- Time step (time bin) (ns).
                    
        Return:     Number of bins in t (integer)                   
    */

    int nBinT = tofMax/tStep;
    std::cout << " GetNumberOfTimeSteps: " << nBinT << std::endl;

    return nBinT;
}

int DataManager::GetBin(int nBins, float lowEdge, float highEdge, float val)
{
    /*
        Calculate Bin Number

        Input:      nBins    -- Number of bins
                    lowEdge  -- Minimum
                    highEdge -- Maximum
                    val      -- Value
                    
        Return:     Bin number                 
    */

    float binWidth = (highEdge - lowEdge) / static_cast<float>(nBins);
    int bin = static_cast<int>( (val - lowEdge) / binWidth);
    return bin;
}

std::vector<std::vector<double>> DataManager::SumEBin(double* binData, int axis, std::vector<int>& numVoxels,
                                                      std::vector<float>& waterBox, std::vector<float>& voxelBox, float offset)
{
    /*
        Bin energy density along axes

        Input:      binData    -- Pointer to binned data
                    axis       -- Axis being plotted against energy
                    numVoxels  -- Vector containing number of voxels
                    waterBox   -- Vector of the dimensions of the phantom
                    voxelBox   -- Vector of the dimensions of a voxel
                    offset     -- Offset for plotting
                    
        Return:     2D vector of summed energy density bins along axes
    */    
    
    int nBinsX = numVoxels[0];
    int nBinsY = numVoxels[1];
    int nBinsZ = numVoxels[2];

    std::vector<std::vector<double>> sumBin(2, std::vector<double>(numVoxels[axis],0));

    if(axis==0)
    {
        for(int width=0; width<nBinsX; width++)
        {
            for(int depth=0; depth<nBinsZ; depth++)
            {
                for(int height=0; height<nBinsY; height++)
                {
                    int index = depth + nBinsZ*(width + nBinsX * height);
                    sumBin[1][width] += binData[index];

                }
            }
            sumBin[0][width] = -waterBox[0] + voxelBox[0]*(0.5 + width) + offset; // center of bin
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
            sumBin[0][height] = -waterBox[1] + voxelBox[1]*(0.5 + height) + offset; // center of bin
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
            sumBin[0][depth] = -waterBox[2] + voxelBox[2]*(0.5 + depth) + offset; // center of bin
        }
    }
    else
        std::cout << "Unknown axis provided, index = " << axis << std::endl;

    return sumBin;
}

void DataManager::PlotGraph(int bins, TString pltTitle, double* hAxis, TString hAxisName, double* vAxis, TString vAxisName)
{
    /*
        Plots graph of given horizontal and vertical axis

        Input:      bins          -- Total number of bins
                    pltTitle      -- Title of plot
                    hAxis         -- Horizontal data
                    hAxisName     -- Title of horizontal axis
                    vAxis         -- Vertical data
                    vAxisName     -- Title of vertical axis
    */    

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
    /*
        3 plots against the 3 cartesian axes

        Input:      binData      -- Energy density data
                    numVoxels    -- Number of voxels
                    waterBox     -- Dimensions of phantom
                    voxelBox     -- Dimensions of voxel
    */    

    // Create canvas
    // **************************************************************
    TCanvas *c = new TCanvas("","",1000,1000,1100,950);
    c->Divide(1,3);
    
    std::vector<float> offset{0., 0., 0.};
        
    // Plotting x
    // **************************************************************
    c->cd(1);
    std::vector<std::vector<double>> sumBin = SumEBin(binData, 0, numVoxels, waterBox, voxelBox, offset[0]);
    double* horBin = &sumBin[0][0];
    double* verBin = &sumBin[1][0];
    TString pltTitle = "Energy density distribution vs binned x";
    TString horAxisName = "x [mm]";
    TString verAxisName = "Energy Density [MeV/mm^{3}]";
    PlotGraph(numVoxels[0],pltTitle,horBin,horAxisName,verBin,verAxisName);
    
    // Plotting y
    // **************************************************************
    c->cd(2);
    sumBin = SumEBin(binData, 1, numVoxels, waterBox, voxelBox, offset[1]);
    horBin = &sumBin[0][0];
    verBin = &sumBin[1][0];
    pltTitle = "Energy density distribution vs binned y";
    horAxisName = "y [mm]";
    verAxisName = "Energy Density [MeV/mm^{3}]";
    PlotGraph(numVoxels[1],pltTitle,horBin,horAxisName,verBin,verAxisName);

    // Plotting z
    // **************************************************************
    c->cd(3);
    sumBin = SumEBin(binData, 2, numVoxels, waterBox, voxelBox, offset[2]);
    horBin = &sumBin[0][0];
    verBin = &sumBin[1][0];
    pltTitle = "Energy density distribution vs binned z";
    horAxisName = "z [mm]";
    verAxisName = "Energy Density [MeV/mm^{3}]";
    PlotGraph(numVoxels[2],pltTitle,horBin,horAxisName,verBin,verAxisName);
    
    // Create file
    // **************************************************************
    c->Print("energydensity-protons.png");
}

void DataManager::Plot3CoordAS(double* binData, std::vector<int>& numVoxels, std::vector<float>& waterBox, std::vector<float>& voxelBox)
{
    /*
        3 plots against the 3 cylindrical axes

        Input:      binData      -- Energy density data
                    numVoxels    -- Number of voxels
                    waterBox     -- Dimensions of cylinder
                    voxelBox     -- Dimensions of cylindrical voxel
    */    

    // Create canvas
    // **************************************************************
    TCanvas *c = new TCanvas("","",1000,1000,1100,950);
    c->Divide(1,3);
        
    std::vector<float> offset{waterBox[0], waterBox[1], 0.};

    // Plotting radius
    // **************************************************************
    c->cd(1);
    std::vector<std::vector<double>> sumBin = SumEBin(binData, 0, numVoxels, waterBox, voxelBox, offset[0]);
    double* horBin = &sumBin[0][0];
    double* verBin = &sumBin[1][0];
    TString pltTitle = "Energy density distribution vs binned radius";
    TString horAxisName = "radius [mm]";
    TString verAxisName = "Energy Density [MeV/mm^{3}]";
    PlotGraph(numVoxels[0],pltTitle,horBin,horAxisName,verBin,verAxisName);

    // Plotting phi
    // **************************************************************
    c->cd(2);
    sumBin = SumEBin(binData, 1, numVoxels, waterBox, voxelBox, offset[1]);
    horBin = &sumBin[0][0];
    verBin = &sumBin[1][0];
    pltTitle = "Energy density distribution vs binned phi";
    horAxisName = "phi [degree]";
    verAxisName = "Energy Density [MeV/mm^{3}]";
    PlotGraph(numVoxels[1],pltTitle,horBin,horAxisName,verBin,verAxisName);

    // Plotting depth
    // **************************************************************
    c->cd(3);
    sumBin = SumEBin(binData, 2, numVoxels, waterBox, voxelBox, offset[2]);
    horBin = &sumBin[0][0];
    verBin = &sumBin[1][0];
    pltTitle = "Energy density distribution vs binned z";
    horAxisName = "z [mm]";
    verAxisName = "Energy Density [MeV/mm^{3}]";
    PlotGraph(numVoxels[2],pltTitle,horBin,horAxisName,verBin,verAxisName);

    // Create file
    // **************************************************************
    c->Print("energydensity-protonsAS.png");
}
