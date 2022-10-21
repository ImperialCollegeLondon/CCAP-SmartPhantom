#ifndef INCLUDE_ROOTIO_HH 
#define INCLUDE_ROOTIO_HH 1

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "SciFiHit.hh"
#include "G4String.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"


class RootIO 
{
public: 
    virtual ~RootIO();
  
    static RootIO* GetInstance(G4String name = "");
    TTree* CreateTree(TString treeName);
    void WriteToRoot(SciFiHitsCollection* hsf, double evtID, TTree* &tree);
    void Close();
    static G4String fOutputFile;
    std::vector<G4String> GetFileExt(const G4String& str);
    G4String GetRootName();
    G4String GetCheckedName() { return checkedName; };
    
protected:
    RootIO(); 
  
private:
    
    TFile* fFile;
    int fNevents;
    double data[8];
    std::string name;
    G4String checkedName;
    
};
#endif // INCLUDE_ROOTIO_HH
