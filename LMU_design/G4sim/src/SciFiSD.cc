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

#include "G4OpticalPhoton.hh"
#include "G4Scintillation.hh"
#include "G4ParticleTable.hh"
#include "EventAction.hh"

#include <iostream>
#include <fstream>
using namespace std;

SciFiSD::SciFiSD(G4String name)
: G4VSensitiveDetector(name), fHitsCollection(0), fHCID(-1)
{
    // Give a name for the collection for all SciFiSD detectors
    G4String HCname = "SciFiColl";
    collectionName.insert(HCname);
}


SciFiSD::~SciFiSD()
{
    collectionName.clear();
}

void SciFiSD::Initialize(G4HCofThisEvent* hce)
{
    // Prepare the hits collection for detector
    fHitsCollection = new SciFiHitsCollection(SensitiveDetectorName,collectionName[0]);
    if (fHCID<0)
    { fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); }
    hce->AddHitsCollection(fHCID,fHitsCollection);
    
    // Fill hits with zero energy GetTotalEnergyDeposit
    SciFiHit* hit = new SciFiHit(0);
    fHitsCollection->insert(hit);    
}

G4bool SciFiSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
    // When proccessing a hit specify what data to collect
    G4double edep = step->GetTotalEnergyDeposit();      // Collect energy deposited
    if (edep==0.) return true;                          // If no energy deposited skip

    // Data related to volume
    G4TouchableHistory* touchable = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* physical = touchable->GetVolume();
    G4int copyNo = physical->GetCopyNo();
            
    G4double steplength = step->GetStepLength();        // Steplength
    
    // The point before and after a step
    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4StepPoint* postStepPoint = step->GetPostStepPoint();
    
    // Randomize energy deposition location to between preStep and postStep
    G4ThreeVector prePoint = preStepPoint->GetPosition();
    G4ThreeVector postPoint = postStepPoint->GetPosition();
    G4ThreeVector worldPos = prePoint + G4UniformRand()*(postPoint - prePoint);
    
    G4double postTime = postStepPoint->GetGlobalTime();                             // Get global time at post step
    G4double deltaT = step->GetDeltaTime();                                         // Get the time of event 
                                                                                    // (i.e. between preStep and postStep)
    
    G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();   // Particle Name causing hit
    //G4int fibreN = GetFibreN(worldPos,step);
    
    // Write information of hit to SciFiHit for SciFiSD
    SciFiHit* hit = new SciFiHit(copyNo,postTime,edep,steplength,worldPos,deltaT,particleName);
    fHitsCollection->insert(hit);

    G4double xPos = worldPos[0];
    G4double yPos = worldPos[1];
    G4double zPos = worldPos[2];

    G4int fibreN = 0;
    G4double fibreSep = 0.300*mm;       // [mm]
    G4int scifiN = 33;  
    G4double fibreRadius = 0.250*mm/2;  // [mm]
    G4double fibreLength = 10*mm;       // [mm] length of fibres in y direction
    G4double centerYPos = 0*mm;         // [mm]
    G4double maxYPos = centerYPos + (fibreLength/2);
    G4double minYPos = centerYPos - (fibreLength/2);
    G4double centerZPos = (-50+13)*mm;  // [mm] z-position of frame

    for (int i = 1; i <= scifiN; ++i) {
        G4double centerXPos = -(fibreSep * floor(scifiN/2)) + (fibreSep*(i-1)); // calculates centre x-position of fibre i
        G4double rho = pow(xPos - centerXPos, 2) + pow(zPos - centerZPos, 2);   // radial position of energy dep relative to axis of fibre i
        if (rho <= fibreRadius && minYPos <= yPos && yPos <= maxYPos) {
            fibreN = i;
        }
    }

    G4double sciEff = 8000;    // photons per meV
    G4double transEff = 0.03;  // percentage efficiency
        
    G4double lightYield = edep * sciEff * transEff;
    
    ofstream myFile; 
    myFile.open("SciFiHitsFibres.dat", ofstream::app);
    //myFile << "Writing this to a file.\n";
    myFile << " " << copyNo << " " << postTime << " " << edep << " " << steplength << 
        " " << worldPos << " " << deltaT << " " << particleName << " " << fibreN << " " << lightYield << " " << endl;
    myFile.close();

    return true;
}
