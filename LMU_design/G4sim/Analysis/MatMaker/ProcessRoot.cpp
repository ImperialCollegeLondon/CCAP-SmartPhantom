#define _USE_MATH_DEFINES

#include "MatFile/MatFile.h"
#include "ProcessRoot.h"
#include "KWaveInput.h"
#include "fftw3.h"

#include "Linterp/linterp.h"

//ProcessRoot::ProcessRoot(const std::string& fName, const std::string& outName, std::vector<float>& voxelSize, std::vector<double>& analysisSize)
ProcessRoot::ProcessRoot(const std::string& fName)
{
    // Specify Name of Trees containing data (i.e. detectors)
    rTree.push_back("Phantom");
    rTree.push_back("Station1_Epoxy");
    rTree.push_back("Station1_PlaneH");
    rTree.push_back("Station1_PlaneV");
    rTree.push_back("Station2_Epoxy");
    rTree.push_back("Station2_PlaneH");
    rTree.push_back("Station2_PlaneV");
    rTree.push_back("Station3_Epoxy");
    rTree.push_back("Station3_PlaneH");
    rTree.push_back("Station3_PlaneV");
    rTree.push_back("Station4_Epoxy");
    rTree.push_back("Station4_PlaneH");
    rTree.push_back("Station4_PlaneV");
    
    debug = true; // Debug output
    
    rootFile = new TFile(fName.c_str(),"READ"); // Open root file
    
    if(debug)
        std::cout << "Opened " << fName << '\n' << std::endl;

    fReader = nullptr;
    fEdep = nullptr;
    fHitTime = nullptr;
    fEvent = nullptr;
    fPosX = nullptr;
    fPosY = nullptr;
    fPosZ = nullptr;
    fDeltaT = nullptr;
    fPName = nullptr;
}

ProcessRoot::~ProcessRoot()
{
    if(debug)
        std::cout << "Beginning destructor for ProcessRoot." << std::endl;
    
    delete rootFile;
    delete fReader;
    delete fEdep;
    delete fHitTime;
    delete fEvent;
    delete fPosX;
    delete fPosY;
    delete fPosZ;
    delete fDeltaT;
    delete fPName;
    delete dataMan;
    delete energyDensitySum;
   
    std::cout << "End of destructor for ProcessRoot.\n" << std::endl;
}

void ProcessRoot::CloseFile()
{
    /*
        Closes ROOT file.
    */
    rootFile->Close();
    std::cout << "Closed ROOT File.\n" << std::endl;
}

void ProcessRoot::Initialise(std::vector<float>& voxelSize, std::vector<double>& analysisSize)
{
    /*
        Initialises some of the variables for 3D
        
        Input:      voxelSize       -- Voxel dimensions
                    analysisSize    -- Analytical volume dimensions
    */

    if(debug)
        std::cout << "Defining Phantom, Voxel Size, and Timesteps." << std::endl;
    
    dataMan = new DataManager();
    
    tofMax = 100.;
    float tStep = 1.;
    nBinsT = dataMan->GetNumberOfTimeSteps(tofMax, tStep); // Doesn't do anything at the moment

    // Specify phantom and voxels
    wBox = dataMan->DefinePhantom(analysisSize[0], analysisSize[1], analysisSize[2]);
    voxBox = dataMan->DefineVoxel(voxelSize[0], voxelSize[1], voxelSize[2]);
    numVox = dataMan->GetNumberOfVoxels(wBox, voxBox);

    nBinsWidth = numVox[0];
    nBinsHeight = numVox[1];
    nBinsDepth = numVox[2];

    // Initialise array size
    energyDensitySum    = new double[nBinsHeight*nBinsWidth*nBinsDepth]{0};
    energySumMat        = new double[nBinsDepth*nBinsHeight*nBinsWidth]{0};
    energyDensitySumMat = new double[nBinsDepth*nBinsHeight*nBinsWidth]{0};
}

void ProcessRoot::ReaderForTree(TString treeName)
{
    /*
        Prepares TTree readers
        
        Input:      treeName       -- Name of tree
    */

    delete fReader;
    delete fEdep;
    delete fHitTime;
    delete fEvent;
    delete fPosX;
    delete fPosY;
    delete fPosZ;
    delete fDeltaT;
    delete fPName;
    
    if(debug)
        std::cout << "Assigning TTreeReader for tree name: " << treeName << std::endl;
    
    TString branch = "Data"; // Branch Name

    // Reassign Reader variables for given tree
    fReader = new TTreeReader(treeName, rootFile);
    fEdep = new TTreeReaderValue<double>(*fReader, branch+".Edep");
    fHitTime = new TTreeReaderValue<double>(*fReader, branch+".HitTime");
    fEvent = new TTreeReaderValue<double>(*fReader, branch+".EventN");
    fPosX = new TTreeReaderValue<double>(*fReader, branch+".PosX");
    fPosY = new TTreeReaderValue<double>(*fReader, branch+".PosY");
    fPosZ = new TTreeReaderValue<double>(*fReader, branch+".PosZ");
    fDeltaT = new TTreeReaderValue<double>(*fReader, branch+".DeltaT");
    fPName = new TTreeReaderValue<std::string>(*fReader, "Name");    
}

void ProcessRoot::BinTreeData()
{
    /*
        Bins data for TTree   
    */
    while(fReader->Next())
    {
        int letters = (*fPName).Get()->size(); // Number of letters in particle name
        
        // Assemble particle name
        std::string pname = "";
        for(int i=0; i<letters; i++)
        {
            pname += (*fPName)->at(i);
        }        
    
        // Bin Data
        //if(pname == "proton") // If you only want interactions involving protons as particle making deposition
                                // i.e. ignoring secondary interactions of other particles like e-
        //{
            if( std::abs(**fPosX) <= wBox[0] && std::abs(**fPosY) <= wBox[1])
            {
                int widthBin = dataMan->GetBin(nBinsWidth, -wBox[0], wBox[0], (**fPosX));
                int heightBin = dataMan->GetBin(nBinsHeight, -wBox[1], wBox[1], (**fPosY));
                int depthBin = dataMan->GetBin(nBinsDepth, -wBox[2], wBox[2], (**fPosZ));
                
                float dV = dataMan->GetVoxVol(voxBox);          // Volume element for 3D
                
                if(heightBin >= 0 && widthBin >= 0 && depthBin >= 0 &&
                    heightBin < nBinsHeight && widthBin < nBinsWidth && depthBin < nBinsDepth) 
                {
                    // Dividing by voxel volume to get energy density
                    energyDensitySum[depthBin + nBinsDepth * (widthBin + nBinsWidth * heightBin)] += (**fEdep)/dV;
                    energySumMat[heightBin + nBinsHeight * (widthBin + nBinsWidth * depthBin)] += (**fEdep);
                    energyDensitySumMat[heightBin + nBinsHeight * (widthBin + nBinsWidth * depthBin)] += (**fEdep)/dV; 
                }
                else
                {
                    std::cout << "Bin Under/Overflow! Ignoring Event." << std::endl;
                }
            }
        //}
    }    
}

void ProcessRoot::ReadFile()
{
    /*
        Read ROOT File for 3D
    */
    if(debug)
        std::cout << "\nBeginning to process root file." << std::endl;
    
    // Bin Data
    for(int i=0; i<rTree.size(); i++)
    {
        ReaderForTree(rTree[i]);
        BinTreeData();
    }

    // Plot energy vs depth
     if(debug)
        std::cout << "\nPlotting energy density against binned axes." << std::endl;
    dataMan->Plot3Coord(energyDensitySum, numVox, wBox, voxBox);

    // Convert energy to Joules
    if(debug)
        std::cout << "\nConverting energy density from MeV/mm^3 to Joules/m^3." << std::endl;
    
    double scalingFactor = 1.60218e-13;           // Conversion factor to Joules
    double siFactor = 1000.*1000.*1000.;          // Conversion factor from mm^3 to m^3
    for(int i=0; i<(nBinsHeight*nBinsWidth*nBinsDepth); i++)
    {
        double scaledEDen = energyDensitySum[i]*scalingFactor*siFactor;
        double scaledEMat = energySumMat[i]*scalingFactor;
        double scaledEDenMat = energyDensitySumMat[i]*scalingFactor*siFactor;
        
        // Arbitrary cut to energy for values < 1e-15 
        if(scaledEDen < 1e-15)
            scaledEDen = 0;
        if(scaledEMat < 1e-15)
            scaledEMat = 0;
        if(scaledEDenMat < 1e-15)
            scaledEDenMat = 0;
        
	    energyDensitySum[i] = scaledEDen;
        energySumMat[i] = scaledEMat;
        energyDensitySumMat[i] = scaledEDenMat;
    }

    if(debug)
        std::cout << "Finished processing root file.\n" << std::endl;
}

void ProcessRoot::WriteMatFile(const std::string& matEnergyName)
{
    /*
        Outputs the Matlab file containing binned energy for 3D
        
        Input:      matEnergyName        -- Name of Matlab file (energy)       
    */

    // Writing Energy Matlab File    
    if(debug) 
        std::cout << "Beginning WriteMatFile to: " << matEnergyName << std::endl;
    
    size_t extIndex = matEnergyName.find_last_of(".");
    const std::string datasetName = matEnergyName.substr(0,extIndex);                           // Name of matlab dataset
    MatFile matFile(matEnergyName);                                                             // Create specified .mat file
    matFile.Write3Dim(datasetName, nBinsDepth, nBinsWidth, nBinsHeight, energySumMat);          // Write data to .mat
    matFile.Close();
    delete energySumMat;

    if(debug)
        std::cout << "Finished WriteMatFile.\n" << std::endl;
}

std::vector<double> ProcessRoot::CentreBeam(double* zData, double* xData, double* yData, int numElements, size_t polyOrder)
{
    /*
        Centre beam for axisymmetric coordinate transformation
    
        Input:      zData       -- Binned z data (depth)
                    xData       -- Binned x data (width)
                    yData       -- Binned y data (height)
                    numElements -- Number of bins
                    polyOrder   -- Polynomial order
        
        Return:     Coefficients to centre beam
    */
   
    PolyFit pFit;
    
    // Input values 
    size_t k = polyOrder;                            // Polynomial order
    bool fixedinter = false;                         // Fixed the intercept (coefficient A0)
    int wtype = 0;                                   // Weight: 0 = none (default), 1 = sigma, 2 = 1/sigma^2
    double fixedinterval = 0.;                       // The fixed intercept value (if applicable)
    size_t n = numElements;
    
    double erry[n] = {0.};

    // Definition of other variables 
    double coefbetaX[k+1];                           // Coefficients of the polynomial (x)
    double coefbetaY[k+1];                           // Coefficients of the polynomial (y)
    double **XTWXInv;                                // Matrix XTWX Inverse [k+1,k+1]
    double **Weights;                                // Matrix Weights [n,n]

    XTWXInv = pFit.Make2DArray(k+1,k+1);
    Weights = pFit.Make2DArray(n, n);

    // Calculate weight (I don't sepcify any)  
    pFit.CalculateWeights(erry, Weights, n, wtype);

    // Polynomial fit to centre beam for transverse planes  
    pFit.Fit(zData,xData,n,k,fixedinter,fixedinterval,coefbetaX,Weights,XTWXInv);
    pFit.Fit(zData,yData,n,k,fixedinter,fixedinterval,coefbetaY,Weights,XTWXInv);
    
    std::vector<double> coeff{coefbetaX[0],coefbetaX[1],coefbetaY[0],coefbetaY[1]};
    std::cout << "Coefficients X: Intercept = " << coefbetaX[0] << '\t' << "Slope = " << coefbetaX[1] << std::endl;
    std::cout << "Coefficients Y: Intercept = " << coefbetaY[0] << '\t' << "Slope = " << coefbetaY[1] << std::endl;
    std::cout << std::endl;
        
    return coeff;
}
