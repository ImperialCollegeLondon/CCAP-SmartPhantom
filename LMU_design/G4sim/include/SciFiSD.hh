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

#ifndef SciFiSD_h
#define SciFiSD_h 1

#include "G4VSensitiveDetector.hh"
#include "DetectorConstruction.hh"
#include "SciFiHit.hh"

#include <iostream>

using namespace std;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

// SciFi Sensitive Detector
// **************************************************************
class SciFiSD : public G4VSensitiveDetector
{
public:
    SciFiSD(G4String name);
    virtual ~SciFiSD();
    
    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);

    // define test function for GetFibreN which outputs integer 1
    // G4int GetFibreN(G4ThreeVector worldPos, G4Step* step) {
        
    //     G4VPhysicalVolume* ThisVol = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    //     G4VPhysicalVolume* NextVol = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
        
    //     G4int scifiN = 33;
    //     G4double fibreSep = 0.300*mm; // define fibre separation (centre to centre distance)
    //     G4double fibreRadius = 0.25*mm/2; // define fibre radius
    //     G4double planeZPos = 3*mm; // define z-position of fibre plane
    //     //G4double fibreLength = 10*mm; //fibre length in the y-direction


    //     G4double xPos = worldPos[0];
    //     G4double yPos = worldPos[1];
    //     G4double zPos = worldPos[2];

    //     G4double maxYPos = 5*mm; // top of fibre
    //     G4double minYPos = -5*mm; // bottom of fibre

    //     G4double maxZPos = planeZPos + fibreRadius;
    //     G4double minZPos = planeZPos - fibreRadius;

    //     if (NextVol && ThisVol->GetName()=="scifiStation1Physical" && NextVol->GetName()=="scifiStation1Physical") {
    //         //return 1;
    //         int i = 1;
    //         G4double centerPos = -(fibreSep * floor(scifiN/2)) + (fibreSep*(i-1));
    //         G4double maxXPos = centerPos + fibreRadius;
    //         G4double minXPos = centerPos - fibreRadius;
    //         // if (xPos >= minXPos && xPos <= maxXPos) {
    //         //     return i;
    //         // }
    //         // else {return 0;}
            
    //         for (int i = 1; i <= 33; ++i) {
    //             G4double centerPos = -(fibreSep * floor(scifiN/2)) + (fibreSep*(i-1));
    //             G4double maxXPos = centerPos + fibreRadius;
    //             G4double minXPos = centerPos - fibreRadius; //change to double not Geant

    //             //myFile << "Writing this to a file.\n";
    //             // cout << " " << centerPos << " " << maxXPos << " " << minXPos << " " << i <<  " " << endl;
    //             cout << " " << centerPos << " " << maxXPos << " " << minXPos << " " << i <<  " " << endl;

    //             if (xPos >= minXPos && xPos <= maxXPos) {
    //                 cout << " " << centerPos << " " << maxXPos << " " << minXPos << " " << i <<  " " << endl;
    //                 return i;
    //                 // put another print statement here
    //                 break;  // return fibre number & exit loop
    //             }
    //             else {continue;
    //             // put another print statement here
    //             }
    //             // continue go back through loop     
                
    //     }
    //     else {return 0;}
    //     }}
        
    


    
private:
    SciFiHitsCollection* fHitsCollection;
    G4int fHCID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
