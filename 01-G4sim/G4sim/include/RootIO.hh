#ifndef INCLUDE_ROOTIO_HH 
#define INCLUDE_ROOTIO_HH 1

// Include files
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "SciFiHit.hh"
#include "G4String.hh"

class RootIO 
{
public: 
  virtual ~RootIO();
  
  static RootIO* GetInstance();
  void Write_WBox(SciFiHitsCollection* hsf, double evtID);
  void Close();

protected:
  RootIO(); 
  
private:
  TTree* waterBox;

  TFile* fFile;
  int fNevents;
  double data[8];
  std::string name;
};
#endif // INCLUDE_ROOTIO_HH
