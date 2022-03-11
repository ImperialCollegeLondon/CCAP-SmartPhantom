#ifndef PROCESSROOT_H
#define PROCESSROOT_H

#include <fstream>
#include <vector>
#include <string>

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

// Class to read from .root files and 
class ProcessRoot
{
public:
    ProcessRoot(const std::string& fName, const std::string& outName, std::vector<float>& voxelSize, std::vector<double>& analysisSize);
    ~ProcessRoot();
    
    void Initialise(std::vector<float>& voxelSize, std::vector<double>& analysisSize);
    void ReadFile();
    void WriteMatFile(const std::string& matFileName, const std::string& datsetName); // Write matlab file
    int GetBinHeight() { return nBinsHeight; };
    int GetBinWidth() { return nBinsWidth; };
    int GetBinDepth() { return nBinsDepth; };
    double* GetEnergySum() { return energySum; };
    void ReaderForTree(TString treeName);
    void BinTreeData();
    
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
    
    TTreeReaderValue<double>* fEdep;
    TTreeReaderValue<double>* fHitTime;
    TTreeReaderValue<double>* fEvent;
    TTreeReaderValue<double>* fPosX;
    TTreeReaderValue<double>* fPosY;
    TTreeReaderValue<double>* fPosZ;
    TTreeReaderValue<double>* fDeltaT;
    TTreeReaderValue<std::string>* fPName;
    
    DataManager* dataMan;
    double* energySum;
    double* energySumMat;

    std::vector<float> wBox;
    std::vector<float> voxBox;
    std::vector<int> numVox;
};
#endif // PROCESSROOT_H
