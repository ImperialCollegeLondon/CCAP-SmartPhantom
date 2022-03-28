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

#include "SciFiSD.hh"
#include "SciFiHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "Randomize.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SciFiSD::SciFiSD(G4String name)
: G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
    /*
        Give a name for the collection for all SciFiSD detectors
    */
    
    G4String HCname = "SciFiColl";
    collectionName.insert(HCname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SciFiSD::~SciFiSD()
{
    collectionName.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SciFiSD::Initialize(G4HCofThisEvent* hce)
{
    /*
        Prepare the hits collection for detector
    */
    
    fHitsCollection = new SciFiHitsCollection(SensitiveDetectorName,collectionName[0]);
    if (fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
    hce->AddHitsCollection(fHCID,fHitsCollection);
    
    // Fill hits with zero energy GetTotalEnergyDeposit
    // **************************************************************
    SciFiHit* hit = new SciFiHit(0);
    fHitsCollection->insert(hit);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SciFiSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    /*
        When proccessing a hit specify what data to collect
    */
    G4double edep = step->GetTotalEnergyDeposit();                                  // Collect energy deposited
    if (edep==0.) return true;                                                      // If no energy deposited skip

    // Data related to volume
    // **************************************************************
    G4TouchableHistory* touchable = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* physical = touchable->GetVolume();
    G4int copyNo = physical->GetCopyNo();
            
    G4double steplength = step->GetStepLength();                                    // Steplength
    
    // The point before and after a step
    // **************************************************************
    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4StepPoint* postStepPoint = step->GetPostStepPoint();
    
    // Randomize energy deposition location to between preStep and postStep
    // **************************************************************
    G4ThreeVector prePoint = preStepPoint->GetPosition();
    G4ThreeVector postPoint = postStepPoint->GetPosition();
    G4ThreeVector worldPos = prePoint + G4UniformRand()*(postPoint - prePoint);
    
    G4double postTime = postStepPoint->GetGlobalTime();                             // Get global time at post step
    G4double deltaT = step->GetDeltaTime();                                         // Get the time of event 
                                                                                    // (i.e. between preStep and postStep)
    
    G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();   // Particle Name causing hit
    
    // Write information of hit to SciFiHit for SciFiSD
    // **************************************************************
    SciFiHit* hit = new SciFiHit(copyNo,postTime,edep,steplength,worldPos,deltaT,particleName);
    fHitsCollection->insert(hit);
    
    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
