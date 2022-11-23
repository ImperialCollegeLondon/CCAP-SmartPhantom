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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include <vector>


#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4GDMLParser.hh"

#include "DetectorConstructionMessenger.hh"

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class DetectorConstructionMessenger;

/// Detector construction

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();
    
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
          
    G4Material* GetNISTMaterial(G4String name, G4String fallbackName = "G4_AIR");

    G4double GetWorldZ() { return worldZ; };
    G4double GetPhantomX();
    G4double GetPhantomY();
    G4double GetPhantomZ();
    
    void CreateSciFiStation(G4LogicalVolume** stationLogical, G4LogicalVolume** stationLogicalH, G4LogicalVolume** stationLogicalV, G4String stationName, std::vector<G4double> &pos);
    
    G4bool GetTogglePhantom() { return togglePhantom; }; 
    G4bool GetToggleScifi() { return toggleScifi; };
    
    // UI Commands
    // **************************************************************
    void SetWorldVolume(G4ThreeVector vec);
    void SetPhantomVolume(G4ThreeVector vec);
    void SetWaterVolume(G4ThreeVector vec);
    void EnablePhantomWall(G4bool toggle);
    void SetPhantomThickness(G4double val);
    void SetWorldMaterial(G4String name);
    void SetPhantomMaterial(G4String name);
    void SetWaterMaterial(G4String name);    
    void EnableScifi(G4bool toggle);
    void SetStation1Pos(G4ThreeVector vector);
    void SetStation2Pos(G4ThreeVector vector);
    void SetStation3Pos(G4ThreeVector vector);
    void SetStation4Pos(G4ThreeVector vector);    
    
private:        
    G4double worldZ, worldX, worldY;
    G4double phantomZ, phantomX, phantomY;
    G4double eWindowZ, eWindowX, eWindowY;
    G4double waterZ, waterX, waterY;
    
    G4double phantomThickness;
    G4bool togglePhantom;
    G4bool toggleScifi;
    
    G4int scifiN;
    G4double scifiFibreRadius;
    G4double scifiPitch;
    G4double scifiLength;
    G4double scifiStationDepth;
    G4double scifiStationSide;
    G4double scifiPlaneSide;
    G4double scifiEpoxySide;

    G4double* scifiDepth;
    
    G4double cwRot;
    G4double ccRot;
    G4RotationMatrix* clockRot;
    G4RotationMatrix* aclockRot;
    G4RotationMatrix* stationRot;
        
    G4VSolid* worldSolid; 
    G4VSolid* eWindowSolid;
    G4VSolid* phantomSolid; 
    G4VSolid* waterSolid;
    
    G4LogicalVolume* worldLogical;
    G4LogicalVolume* phantomLogical; 
    G4LogicalVolume* eWindowLogical;
    G4LogicalVolume* waterLogical;
    
    G4LogicalVolume* scifiStation1Logical;
    G4LogicalVolume* scifiStation1LogicalHor;
    G4LogicalVolume* scifiStation1LogicalVer;

    G4LogicalVolume* scifiStation2Logical;    
    G4LogicalVolume* scifiStation2LogicalHor;
    G4LogicalVolume* scifiStation2LogicalVer;

    G4LogicalVolume* scifiStation3Logical;    
    G4LogicalVolume* scifiStation3LogicalHor;
    G4LogicalVolume* scifiStation3LogicalVer;

    G4LogicalVolume* scifiStation4Logical;
    G4LogicalVolume* scifiStation4LogicalHor;
    G4LogicalVolume* scifiStation4LogicalVer;
    
    std::vector<double> station1Pos;
    std::vector<double> station2Pos;
    std::vector<double> station3Pos;
    std::vector<double> station4Pos;
    
    DetectorConstructionMessenger* m_detectorConstructionMessenger;
    
    G4Material* worldMaterial;
    G4Material* phantomMaterial;
    G4Material* waterMaterial;
    G4Material* fibreMaterial;
    G4Material* epoxyResinMaterial;
    
    G4VisAttributes* worldAttr;
    G4VisAttributes* phantomAttr; 
    G4VisAttributes* eWindowAttr;
    G4VisAttributes* waterAttr;
    
    G4bool checkOverlaps;
    
    G4VSensitiveDetector* phantomSD;
    G4VSensitiveDetector* scifi1EpoxySD;
    G4VSensitiveDetector* scifi1HSD;
    G4VSensitiveDetector* scifi1VSD;
    G4VSensitiveDetector* scifi2EpoxySD;
    G4VSensitiveDetector* scifi2HSD;
    G4VSensitiveDetector* scifi2VSD;
    G4VSensitiveDetector* scifi3EpoxySD;
    G4VSensitiveDetector* scifi3HSD;
    G4VSensitiveDetector* scifi3VSD;
    G4VSensitiveDetector* scifi4EpoxySD;
    G4VSensitiveDetector* scifi4HSD;
    G4VSensitiveDetector* scifi4VSD;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
