#ifndef PROCESSROOT_H
#define PROCESSROOT_H

#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <cmath>

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TPaveStats.h"
#include "DataManager.h"
#include "H5Cpp.h"
#include "Polyfit/PolyFit.h"

// Class to read and analyse ROOT files from SmartPhantom (Geant4) 
class ProcessRoot
{
public:
    //ProcessRoot(const std::string& fName, const std::string& outName, std::vector<float>& voxelSize, std::vector<double>& analysisSize);
    ProcessRoot(const std::string& fName);
    ~ProcessRoot();
    
    void Initialise(std::vector<float>& voxelSize, std::vector<double>& analysisSize);

    void ReadFile();
    void CloseFile();

    void WriteMatFile(const std::string& matEnergyName);
    
    int GetBinHeight() { return nBinsHeight; };
    int GetBinWidth() { return nBinsWidth; };
    int GetBinDepth() { return nBinsDepth; };
    int GetBinPhi() { return nBinsPhi; };
    int GetBinRadius() { return nBinsR; };
    double* GetEnergyDensity() { return energyDensitySum; };
    
    void ReaderForTree(TString treeName);
    void BinTreeData();
    std::vector<double> CentreBeam(double* zData, double* xData, double* yData, int numElement, size_t polyOrder);
    
private:
    std::string fileName;
    TFile *rootFile;
    bool debug;
    
    TTreeReader* fReader;
    std::vector<TString> rTree;

    float tofMax;
    int nBinsHeight;
    int nBinsWidth;
    int nBinsDepth; 
    int nBinsT;
    int nBinsR;
    int nBinsPhi;
    
    TTreeReaderValue<double>* fEdep;
    TTreeReaderValue<double>* fHitTime;
    TTreeReaderValue<double>* fEvent;
    TTreeReaderValue<double>* fPosX;
    TTreeReaderValue<double>* fPosY;
    TTreeReaderValue<double>* fPosZ;
    TTreeReaderValue<double>* fDeltaT;
    TTreeReaderValue<std::string>* fPName;
    
    DataManager* dataMan;
    double* energyDensitySum;
    double* energySumMat;
    double* energyDensitySumMat;

    std::vector<float> wBox;
    std::vector<float> voxBox;
    std::vector<int> numVox;

};
#endif // PROCESSROOT_H
