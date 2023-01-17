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

#ifndef SciFiHit_h
#define SciFiHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

class SciFiHit : public G4VHit
{
public:
    SciFiHit();
    SciFiHit(G4int z);
    SciFiHit(G4int i,G4double t,G4double e,G4double steplength, G4ThreeVector pos, G4double deltaT, G4String name, G4int fibreN, G4double lightYield);
    virtual ~SciFiHit();
        
    // Volume Copy Number
    G4int GetID() const {return fID; }
    void SetTime(G4double val) { fTime = val; }
    G4double GetTime() const { return fTime; }
    
    // Edep
    void SetEdep(G4double de) { fEdep = de; }
    void AddEdep(G4double de) { fEdep += de; }
    G4double GetEdep() const { return fEdep; }
    
    // StepLength
    G4double GetStepLength() const { return fSteplength; }
    
    // Position
    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }
    
    // Delta time
    void SetDeltaT(G4double val) { fDeltaT = val; }
    G4double GetDeltaT() const { return fDeltaT; }
    
    // Name
    void SetName(G4String val) { fName = val; }
    G4String GetName() const { return fName; }

    // Fibre number
    void SetFibreN(G4int num) { fFibreN = num; }
    G4int GetFibreN() const { return fFibreN; }

    // Light Yield
    void SetLightYield( G4double yield ) {fLightYield = yield; }
    G4double GetLightYield() const { return fLightYield; }

                
private:
    G4int fID;
    G4double fTime;
    G4double fEdep;
    G4double fSteplength;
    G4ThreeVector fPos;
    G4double fDeltaT;
    G4String fName;
    G4int fFibreN;
    G4int fLightYield;
};

typedef G4THitsCollection<SciFiHit> SciFiHitsCollection;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
