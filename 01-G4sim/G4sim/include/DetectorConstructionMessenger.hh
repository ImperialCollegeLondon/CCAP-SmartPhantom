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

#ifndef DetectorConstructionMessenger_h
#define DetectorConstructionMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4Types.hh"
#include "DetectorConstruction.hh"

class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstructionMessenger : public G4UImessenger
{
public:
    DetectorConstructionMessenger(DetectorConstruction* detector);
    
    ~DetectorConstructionMessenger();
    
    void SetNewValue(G4UIcommand* command, G4String newValue) override; 
    
private:
    DetectorConstruction* m_detector;
    G4UIdirectory* m_detDir;
    G4UIcmdWith3VectorAndUnit* m_setWorldVolumeCmd;
    G4UIcmdWith3VectorAndUnit* m_setPhantomVolumeCmd;
    G4UIcmdWith3VectorAndUnit* m_setWaterVolumeCmd;
    G4UIcmdWithABool* m_enablePhantomCmd;
    G4UIcmdWithADoubleAndUnit* m_setPhantomThicknessCmd;
    G4UIcmdWithAString* m_setWorldMaterialCmd;
    G4UIcmdWithAString* m_setPhantomMaterialCmd;
    G4UIcmdWithAString* m_setWaterMaterialCmd;
    G4UIcmdWithABool* m_enableScifiCmd;
    G4UIcmdWith3VectorAndUnit* m_setStation1PosCmd;
    G4UIcmdWith3VectorAndUnit* m_setStation2PosCmd;
    G4UIcmdWith3VectorAndUnit* m_setStation3PosCmd;
    G4UIcmdWith3VectorAndUnit* m_setStation4PosCmd;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
