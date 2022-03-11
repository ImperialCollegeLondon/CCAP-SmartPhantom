#include "MatFile/MatFile.h"
#include "ProcessRoot.h"
#include "KWaveInput.h"
#include "fftw3.h"

#include "Linterp/linterp.h"

ProcessRoot::ProcessRoot(const std::string& fName, const std::string& outName, std::vector<float>& voxelSize, std::vector<double>& analysisSize)
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
    const std::string datasetName = "energy_density_data2_cpp"; // Name of matlab dataset
    
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
    
    Initialise(voxelSize, analysisSize);
    ReadFile();
    WriteMatFile(outName, datasetName);
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
    delete energySum;
   
    std::cout << "End of destructor for ProcessRoot.\n" << std::endl;
}

void ProcessRoot::Initialise(std::vector<float>& voxelSize, std::vector<double>& analysisSize)
{
/*
    if(debug)
        std::cout << "Assigning TTreeReader and TTreeReaderValue for Geometry.\n" << std::endl;
    
    // Reader for Geometry tree (phantom geometry)
    TTreeReader fReaderGeom("Geometry", rootFile);
    
    TTreeReaderValue<double> fPhantomHeight(fReaderGeom, "Size.Height");
    TTreeReaderValue<double> fPhantomWidth(fReaderGeom, "Size.Width");
    TTreeReaderValue<double> fPhantomDepth(fReaderGeom, "Size.Depth");

    double phantomHeight, phantomWidth, phantomDepth;
    
    while(fReaderGeom.Next())
    {
        phantomHeight = *fPhantomHeight;
        phantomWidth = *fPhantomWidth;
        phantomDepth = *fPhantomDepth;
    }
*/
    
    if(debug)
        std::cout << "Defining Phantom, Voxel Size, and Timesteps." << std::endl;
    
    dataMan = new DataManager();
    
    tofMax = 100.;
    float tStep = 1.;
    
    wBox = dataMan->DefinePhantom(analysisSize[0], analysisSize[1], analysisSize[2]);
    voxBox = dataMan->DefineVoxel(voxelSize[0], voxelSize[1], voxelSize[2]);
    numVox = dataMan->GetNumberOfVoxels(wBox, voxBox);
    
    nBinsWidth = numVox[0];
    nBinsHeight = numVox[1];
    nBinsDepth = numVox[2];
    nBinsT = dataMan->GetNumberOfTimeSteps(tofMax, tStep); // Doesn't do anything at the moment

    energySum = new double[nBinsHeight*nBinsWidth*nBinsDepth]{0};
    energySumMat = new double[nBinsDepth*nBinsHeight*nBinsWidth]{0};
}

void ProcessRoot::ReaderForTree(TString treeName)
{
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
    while(fReader->Next())
    {
        int letters = (*fPName).Get()->size(); // Number of letters in particle name
        
        // Assemble particle name
        std::string pname = "";
        for(int i=0; i<letters; i++)
        {
            pname += (*fPName)->at(i);
        }        
    
        //if(pname == "proton") // If you only want interactions involving protons as particle making deposition
                                // i.e. ignoring secondary interactions of other particles like e-
        //{
            if( std::abs(**fPosX) <= wBox[0] && std::abs(**fPosY) <= wBox[1]) // Not too sure why..
            {
                int widthBin = dataMan->GetBin(nBinsWidth, -wBox[0], wBox[0], (**fPosX));
                int heightBin = dataMan->GetBin(nBinsHeight, -wBox[1], wBox[1], (**fPosY));
                int depthBin = dataMan->GetBin(nBinsDepth, -wBox[2], wBox[2], (**fPosZ));
                if(heightBin >= 0 && widthBin >= 0 && depthBin >= 0) 
                {
                    energySum[depthBin + nBinsDepth * (widthBin + nBinsWidth * heightBin)] += (**fEdep);
                    //energySumMat[widthBin + nBinsWidth * (heightBin + nBinsHeight * depthBin)] += (**fEdep);
                    energySumMat[heightBin + nBinsHeight * (widthBin + nBinsWidth * depthBin)] += (**fEdep); 
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
    if(debug)
        std::cout << "\nBeginning to process root file." << std::endl;
    
    for(int i=0; i<rTree.size(); i++)
    {
        ReaderForTree(rTree[i]);
        BinTreeData();
    }

    // Bin and Plot energy vs depth
     if(debug)
        std::cout << " Plotting energy against binned axes." << std::endl;
    dataMan->Plot3Coord(energySum, numVox, wBox, voxBox);

    if(debug)
        std::cout << " Converting energy density from MeV to Joules." << std::endl;
    double scalingFactor = 1.60218e-13;
    for(int i=0; i<(nBinsHeight*nBinsWidth*nBinsDepth); i++)
    {
        double scaledE = energySum[i]*scalingFactor;
        double scaledEMat = energySumMat[i]*scalingFactor;
        
        if(scaledE < 1e-15)
        {
            scaledE = 0;
        }
        if(scaledEMat < 1e-15)
        {
            scaledEMat = 0;
        }
        
	    energySum[i] = scaledE;
        energySumMat[i] = scaledEMat;
    }
    
    if(debug)
        std::cout << "Finished processing root file.\n" << std::endl;
}

void ProcessRoot::WriteMatFile(const std::string& matFileName, const std::string& datasetName)
{
    if(debug) 
        std::cout << "Beginning WriteMatFile to: " << matFileName << std::endl;
    
    MatFile matFile(matFileName); // Create specified .mat file
    matFile.Write3Dim(datasetName, nBinsDepth, nBinsWidth, nBinsHeight, energySumMat); // Write data to .mat
    matFile.Close();
    delete energySumMat;

    if(debug)
        std::cout << "Finished WriteMatFile.\n" << std::endl;
} 
