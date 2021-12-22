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

#include "DetectorConstruction.hh"
#include "SciFiSD.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"
#include "G4Trd.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(), 
  fVisAttributes(), waterboxLogical(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    for (G4int i=0; i<G4int(fVisAttributes.size()); ++i) 
    {
      delete fVisAttributes[i];
    }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    /////////////////////////
    // Construct Materials //
    /////////////////////////
    ConstructMaterials();
    G4Material* water = G4Material::GetMaterial("G4_WATER");
    G4Material* air = G4Material::GetMaterial("G4_AIR");
        
    // Option to switch on/off checking of volumes overlaps
    G4bool checkOverlaps = false;
    
    ///////////////////////////////
    // Z-Positions of Components //
    ///////////////////////////////
    G4double world_z = 300*mm/2; // half-z of world
    G4double water_z = 300*mm/2; // half-z of phantom
        
    //////////////////////
    // Defining Volumes //
    //////////////////////
    
    // World
    G4VSolid* worldSolid = new G4Box("worldBox",15*cm,15*cm,world_z); // Format: G4Box(name, half-x, half-y, half-z)
    G4LogicalVolume* worldLogical = new G4LogicalVolume(worldSolid,air,"worldLogical");
    G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,false,0,checkOverlaps);

    G4VisAttributes* worldAttributes = new G4VisAttributes(G4Colour(0,0,0));
    worldAttributes->SetVisibility(false);
    worldLogical->SetVisAttributes(worldAttributes);
    
    // Water Box
    G4VSolid* waterboxSolid = new G4Box("waterBox",15*cm,15*cm,water_z); 
    waterboxLogical = new G4LogicalVolume(waterboxSolid,water,"waterboxLogical");
    G4VPhysicalVolume* waterboxPhysical = new G4PVPlacement(0,G4ThreeVector(),waterboxLogical,"waterboxPhysical",worldLogical,false,0,checkOverlaps);

    /////////////////////
    // G4VisAttributes //
    /////////////////////
    G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0,0,1,0.15));
    waterboxLogical->SetVisAttributes(visAttributes);
    fVisAttributes.push_back(visAttributes);
    
    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{

    // sensitive detectors -----------------------------------------------------
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    
    G4VSensitiveDetector* waterSD = new SciFiSD(SDname="/waterBoxSD");
    SDman->AddNewDetector(waterSD);
    waterboxLogical->SetSensitiveDetector(waterSD);

}    


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructMaterials()
{
    G4NistManager* nistManager = G4NistManager::Instance();

    // Air 
    nistManager->FindOrBuildMaterial("G4_AIR");
    
    // Water
    nistManager->FindOrBuildMaterial("G4_WATER");   
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::FindMaterial(G4String name) {
    G4Material* material = G4Material::GetMaterial(name,true);
    return material;
}
