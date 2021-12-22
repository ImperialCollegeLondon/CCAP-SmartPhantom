//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "RootIO.hh"
//
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static RootIO* instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootIO::RootIO():fNevents(0)
{
  // initialize ROOT
  TSystem ts;

  if (instance == 0 )
  {
    fFile = new TFile("waterhits.root","RECREATE");
    
    waterBox = new TTree("Event","Event");
    waterBox->Branch("WaterBox", &data, "Edep/D:HitTime/D:EventN/D:PosX/D:PosY/D:PosZ/D:Steplength/D:DeltaT/D"); // Branch to store numerical data
    waterBox->Branch("Name", &name); // Branch to store particle names
    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootIO::~RootIO()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootIO* RootIO::GetInstance()
{
  if (instance == 0 )
  {
    instance = new RootIO();
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootIO::Write_WBox(SciFiHitsCollection* hsf, double evtID)
{
    int n_hit = hsf->entries();
    for(int i=0;i<n_hit;i++){
        SciFiHit* hit = (*hsf)[i];
        double eDep = hit->GetEdep();
        if(eDep>0){
            
            data[0] = eDep;
            data[1] = hit->GetTime();
            data[2] = evtID;
            data[3] = hit->GetPos().x();
            data[4] = hit->GetPos().y();
            data[5] = hit->GetPos().z();
            data[6] = hit->GetStepLength();
            data[7] = hit->GetDeltaT();
            name = hit->GetName();
            
            waterBox->Fill();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootIO::Close()
{ 
  fFile->Write();
  fFile->Close();
  delete fFile;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
