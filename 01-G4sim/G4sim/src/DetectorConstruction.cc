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
#include "G4String.hh"

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
#include "G4RunManagerKernel.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(), 
  worldAttr(nullptr), phantomAttr(nullptr), eWindowAttr(nullptr), waterAttr(nullptr),
  worldSolid(nullptr), phantomSolid(nullptr), eWindowSolid(nullptr), waterSolid(nullptr), worldLogical(nullptr),
  phantomLogical(nullptr), eWindowLogical(nullptr), waterLogical(nullptr), togglePhantom(false), toggleScifi(false),
  worldMaterial(nullptr), phantomMaterial(nullptr), waterMaterial(nullptr), checkOverlaps(false),
  phantomSD(nullptr), scifi1EpoxySD(nullptr), scifi1HSD(nullptr), scifi1VSD(nullptr), scifi2EpoxySD(nullptr),
  scifi2HSD(nullptr), scifi2VSD(nullptr), scifi3EpoxySD(nullptr), scifi3HSD(nullptr), scifi3VSD(nullptr),
  scifi4EpoxySD(nullptr), scifi4HSD(nullptr), scifi4VSD(nullptr)
{
    /*
        Initialises a bunch of parameters for the simulation
    */
    
    // Defining Dimensions
    // **************************************************************
    worldZ = 300*mm/2;                                              // half-z of world
    worldX = worldZ;                                                // half-x of world
    worldY = worldZ;                                                // half-y of world
    
    waterZ = 300*mm/2;                                              // half-z of water
    waterX = waterZ;                                                // half-x of water
    waterY = waterZ;                                                // half-y of water

    phantomZ = waterZ;                                              // half-z of phantom
    phantomX = phantomZ;                                            // half-x of phantom
    phantomY = phantomZ;                                            // half-y of phantom
    phantomThickness = 10*mm;                                       // Full thickness, overrides water dimensions 
                                                                    // if togglePhantom is true

    eWindowZ = 3.05*mm/2;                                           // half-z of entrance window
    eWindowX = 170*mm/2;                                            // half-x of entrance window
    eWindowY = 170*mm/2;                                            // half-y of entrance window
    
    // SmartPhantom Planes
    // **************************************************************
    scifiN = 492;                                                   // Number of fibres
    scifiFibreRadius = 0.25*mm/2;                                   // Half-radius of fibre
    scifiPitch = 0.305*mm;                                          // Fibre pitch
    scifiLength = (scifiPitch*(scifiN/2)+scifiFibreRadius);         // half-length of plane
    scifiStationDepth = 0.8964*mm/2;                                // half-depth of a Station (2 Planes -> 1 Station)
    scifiStationSide = scifiLength;                                 // Transverse face edge length of station
    scifiEpoxySide = scifiPlaneSide;                                // Match epoxy resin layer to plane edge length
    
    // Rotating Stations
    // **************************************************************
    G4double theta = 90*deg;
    cwRot = 0*deg;
    ccRot = 90*deg;
    clockRot = new G4RotationMatrix();                              // Fibre Rotation (clockwise)
    clockRot->rotateX(theta);
    clockRot->rotateY(cwRot);
    aclockRot = new G4RotationMatrix();                             // Fibre Rotation (counter clockwise)
    aclockRot->rotateX(theta);
    aclockRot->rotateY(ccRot);
    stationRot = new G4RotationMatrix();                            // Station Rotation
    stationRot->rotateZ(45*deg);
    
    // Suitable for ~200 MeV Protons
    // **************************************************************
    std::vector<G4double> vec1 = {0*mm, 0*mm, -waterZ + 12.5*mm};
    std::vector<G4double> vec2 = {0*mm, 0*mm, -waterZ + 235*mm};
    std::vector<G4double> vec3 = {0*mm, 0*mm, -waterZ + 240*mm};
    std::vector<G4double> vec4 = {0*mm, 0*mm, -waterZ + 245*mm};

    station1Pos = vec1;                                             // Station 1 depth (at station centre)
    station2Pos = vec2;                                             // Station 2 depth (at station centre)
    station3Pos = vec3;                                             // Station 3 depth (at station centre)
    station4Pos = vec4;                                             // Station 4 depth (at station centre)

    // Construct Materials
    // **************************************************************
    worldMaterial = GetNISTMaterial("G4_AIR");                      // Default material for the world volume is G4_AIR
    phantomMaterial = GetNISTMaterial("G4_PLEXIGLASS");             // Default material for phantom walls is G4_PLEXIGLASS
    waterMaterial = GetNISTMaterial("G4_WATER");                    // Default material for water volume is G4_WATER
    fibreMaterial = GetNISTMaterial("G4_POLYSTYRENE");              // Default material for scintillating fibres is G4_POLYSTYRENE

    G4Element* H = new G4Element("Hydrogen","H",1.,1.01*g/mole);
    G4Element* C = new G4Element("Carbon","C",6.,12.01*g/mole);
    G4Element* O = new G4Element("Oxygen","O",8.,16.00*g/mole);
    G4Element* Cl = new G4Element("Chlorine","Cl",17.,35.453*g/mole);

    // Create Epoxy Resin from Components (based on Araldite I think?)
    // Used for the glue layer in the station
    // **************************************************************
    G4double epoxyResinDensity = 1.2*g/cm3;
    epoxyResinMaterial = new G4Material("EpoxyResin",epoxyResinDensity, 4);
    epoxyResinMaterial->AddElement(C, 21);
    epoxyResinMaterial->AddElement(H, 25);
    epoxyResinMaterial->AddElement(Cl, 1);
    epoxyResinMaterial->AddElement(O, 5);
        
    // Initialise Messenger
    // **************************************************************
    m_detectorConstructionMessenger = new DetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
    delete clockRot;
    delete aclockRot;
    delete stationRot;
    
    delete epoxyResinMaterial;
    
    delete worldSolid;
    delete waterSolid;
    delete eWindowSolid;
    delete phantomSolid;
    
    delete worldLogical;
    delete waterLogical;
    delete eWindowLogical;
    delete phantomLogical;

    delete worldAttr;
    delete phantomAttr;
    delete eWindowAttr;
    delete waterAttr;

    delete m_detectorConstructionMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    /*
        Constructs the various volumes used in the simulation
    */
    
    delete worldSolid;
    delete waterSolid;
    delete eWindowSolid;
    delete phantomSolid;
    
    delete worldLogical;
    delete waterLogical;
    delete eWindowLogical;
    delete phantomLogical;
        
    delete worldAttr;
    delete phantomAttr;
    delete eWindowAttr;
    delete waterAttr;
          
    //////////////////////
    // Defining Volumes //
    //////////////////////
    
    // World
    // **************************************************************
    worldSolid = new G4Box("worldBox", worldX, worldY, worldZ); // Format: G4Box(name, half-x, half-y, half-z)
    worldLogical = new G4LogicalVolume(worldSolid,worldMaterial,"worldLogical");
    G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,false,0,checkOverlaps);
    worldAttr = new G4VisAttributes(G4Colour(0,0,0));
    worldAttr->SetVisibility(false);
    worldLogical->SetVisAttributes(worldAttr);

    if( togglePhantom ) // If phantom walls are enabled
    {
        // Phantom
        // **************************************************************
        phantomSolid = new G4Box("phantomBox", phantomX, phantomY, phantomZ);
        phantomLogical = new G4LogicalVolume(phantomSolid, phantomMaterial, "phantomLogical");
        new G4PVPlacement(0,G4ThreeVector(),phantomLogical,"phantomPhysical",worldLogical,false,0,checkOverlaps);
        phantomAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5,0.15));
        phantomLogical->SetVisAttributes(phantomAttr);
        
        // Entrance Window (Recessed part of wall)
        // **************************************************************
        eWindowSolid = new G4Box("eWindowBox", eWindowX, eWindowY, eWindowZ);
        eWindowLogical = new G4LogicalVolume(eWindowSolid, worldMaterial, "eWindowLogical");
        new G4PVPlacement(0,G4ThreeVector(0,0,-phantomZ+eWindowZ),eWindowLogical,"eWindowPhysical",phantomLogical,false,0,checkOverlaps);
        eWindowAttr = new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.15));
        eWindowLogical->SetVisAttributes(eWindowAttr);
        
        // Water Box
        // **************************************************************
        waterX = phantomX - phantomThickness;
        waterY = phantomY - phantomThickness;
        waterZ = phantomZ - phantomThickness;
        
        waterSolid = new G4Box("waterBox", waterX, waterY, waterZ); 
        waterLogical = new G4LogicalVolume(waterSolid,waterMaterial,"waterboxLogical");
        new G4PVPlacement(0,G4ThreeVector(),waterLogical,"waterboxPhysical",phantomLogical,false,0,checkOverlaps);
        waterAttr = new G4VisAttributes(G4Colour(0,0,1,0.15));
        waterLogical->SetVisAttributes(waterAttr);
    }
    else // If phantom walls not enabled
    {
        // Water Box
        // **************************************************************
        waterSolid = new G4Box("waterBox", waterX, waterY, waterZ); 
        waterLogical = new G4LogicalVolume(waterSolid,waterMaterial,"waterboxLogical");
        new G4PVPlacement(0,G4ThreeVector(),waterLogical,"waterboxPhysical",worldLogical,false,0,checkOverlaps);
        waterAttr = new G4VisAttributes(G4Colour(0,0,1,0.15));
        waterLogical->SetVisAttributes(waterAttr);
    }

    if( toggleScifi ) // If scifi detectors enabled
    {
        CreateSciFiStation(&scifiStation1Logical,&scifiStation1LogicalHor,&scifiStation1LogicalVer,"scifiStation1",station1Pos);
        CreateSciFiStation(&scifiStation2Logical,&scifiStation2LogicalHor,&scifiStation2LogicalVer,"scifiStation2",station2Pos);
        CreateSciFiStation(&scifiStation3Logical,&scifiStation3LogicalHor,&scifiStation3LogicalVer,"scifiStation3",station3Pos);
        CreateSciFiStation(&scifiStation4Logical,&scifiStation4LogicalHor,&scifiStation4LogicalVer,"scifiStation4",station4Pos);
    }
        
    return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::CreateSciFiStation(G4LogicalVolume** scifiStationLogical, G4LogicalVolume** scifiStationLogicalH, G4LogicalVolume** scifiStationLogicalV, G4String stationName, std::vector<G4double> &pos)
{    
    /*
        Creates and assigns volumes for SciFi stations
        
        Input:      scifiStationLogical      -- Logical Volume for Station
                    scifiStationLogicalH     -- Logical volume for horizontal plane of station
                    scifiStationLogicalV     -- Logical volume for vertical plane of station
                    stationName              -- Name of the station
                    pos                      -- Centre position of the station
    */
    
    G4VisAttributes* scifiStationAttributes = new G4VisAttributes(G4Colour(1,0.85,0,0.95));
    G4VisAttributes* scifiFibreAttributes = new G4VisAttributes(G4Colour(0,1,0,0.3));
    
    // SciFi Station
    // **************************************************************
    G4VSolid* scifiStationSolid = new G4Box(stationName,scifiStationSide,scifiStationSide,scifiStationDepth);
    
    G4String stationLogicalName = stationName + "Logical";
    *scifiStationLogical = new G4LogicalVolume(scifiStationSolid,epoxyResinMaterial,stationLogicalName);
    G4String stationPhysicalName = stationName + "Physical";
    new G4PVPlacement(stationRot,G4ThreeVector(pos[0],pos[1],pos[2]),*scifiStationLogical,stationPhysicalName,waterLogical,false,0,checkOverlaps);
    
    // Assign visual attributes
    // **************************************************************
    (*scifiStationLogical)->SetVisAttributes(scifiStationAttributes);

    // Horizontal and Vertical Fibre Planes
    // **************************************************************
    G4String stationFibreSolidName = stationName + "Fibre";
    G4String stationFibreLogicalHName = stationName + "LogicalHor";
    G4String stationFibreLogicalVName = stationName + "LogicalVer";
    G4String stationFibrePhysicalHName = stationName + "PhysicalHor";
    G4String stationFibrePhysicalVName = stationName + "PhysicalVer";

    G4VSolid* scifiStationFibreSolid = new G4Tubs(stationFibreSolidName,0,scifiFibreRadius,scifiLength,0*deg,360.0*deg);
    *scifiStationLogicalH = new G4LogicalVolume(scifiStationFibreSolid,fibreMaterial,stationFibreLogicalHName);
    *scifiStationLogicalV = new G4LogicalVolume(scifiStationFibreSolid,fibreMaterial,stationFibreLogicalVName);

    // Loop to create individual fibres
    // **************************************************************
    for(G4int i=0; i<scifiN; i++)
    {
        G4double x = -(scifiPitch * scifiN/2) + scifiPitch*i; // Centre position of fibre
        G4double xHor = x*cos(cwRot); // Rotate fibre
        G4double xVer = x*cos(ccRot);
        
        G4double yHor = x*sin(cwRot);
        G4double yVer = x*sin(ccRot);
        
        G4double zHor = -scifiFibreRadius;
        G4double zVer = scifiFibreRadius;
        
        new G4PVPlacement(clockRot,G4ThreeVector(xHor,yHor,zHor),*scifiStationLogicalH,stationFibrePhysicalHName,*scifiStationLogical,false,i,checkOverlaps);
        new G4PVPlacement(aclockRot,G4ThreeVector(xVer,yVer,zVer),*scifiStationLogicalV,stationFibrePhysicalVName,*scifiStationLogical,false,i,checkOverlaps);            
    }
    
    (*scifiStationLogicalH)->SetVisAttributes(scifiFibreAttributes);
    (*scifiStationLogicalV)->SetVisAttributes(scifiFibreAttributes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
    /*
        Create/Assign Sensitive Detectors (SD) i.e. volumes where data recorded
    */
    
    delete phantomSD;
    delete scifi1EpoxySD;
    delete scifi1HSD;
    delete scifi1VSD;
    delete scifi2EpoxySD;
    delete scifi2HSD;
    delete scifi2VSD;
    delete scifi3EpoxySD;
    delete scifi3HSD;
    delete scifi3VSD;
    delete scifi4EpoxySD;
    delete scifi4HSD;
    delete scifi4VSD;
 
    // Sensitive detectors
    // **************************************************************
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    
    phantomSD = new SciFiSD(SDname="/phantomSD");
    SDman->AddNewDetector(phantomSD);

    waterLogical->SetSensitiveDetector(phantomSD);
    
    if(togglePhantom) // If phantom walls enabled add to phantomSD
    {
        phantomLogical->SetSensitiveDetector(phantomSD);
        eWindowLogical->SetSensitiveDetector(phantomSD);
    }

    if( toggleScifi ) // If SciFi included add create these sensitive detectors
    {
        // Create sensitive detectors for stations
        // **************************************************************
        scifi1EpoxySD = new SciFiSD(SDname="/scifi1EpoxySD");
        scifi1HSD = new SciFiSD(SDname="/scifi1HSD");
        scifi1VSD = new SciFiSD(SDname="/scifi1VSD");
    
        scifi2EpoxySD = new SciFiSD(SDname="/scifi2EpoxySD");
        scifi2HSD = new SciFiSD(SDname="/scifi2HSD");
        scifi2VSD = new SciFiSD(SDname="/scifi2VSD");
        
        scifi3EpoxySD = new SciFiSD(SDname="/scifi3EpoxySD");
        scifi3HSD = new SciFiSD(SDname="/scifi3HSD");
        scifi3VSD = new SciFiSD(SDname="/scifi3VSD");
        
        scifi4EpoxySD = new SciFiSD(SDname="/scifi4EpoxySD");
        scifi4HSD = new SciFiSD(SDname="/scifi4HSD");
        scifi4VSD = new SciFiSD(SDname="/scifi4VSD");

        // Add to the G4SDManager
        // **************************************************************
        SDman->AddNewDetector(scifi1EpoxySD);
        SDman->AddNewDetector(scifi1HSD);
        SDman->AddNewDetector(scifi1VSD);
        SDman->AddNewDetector(scifi2EpoxySD);
        SDman->AddNewDetector(scifi2HSD);
        SDman->AddNewDetector(scifi2VSD);
        SDman->AddNewDetector(scifi3EpoxySD);
        SDman->AddNewDetector(scifi3HSD);
        SDman->AddNewDetector(scifi3VSD);
        SDman->AddNewDetector(scifi4EpoxySD);
        SDman->AddNewDetector(scifi4HSD);
        SDman->AddNewDetector(scifi4VSD);
        
        // Assign a logical volume to the sensitive detectors
        // **************************************************************
        scifiStation1Logical->SetSensitiveDetector(scifi1EpoxySD);
        scifiStation1LogicalHor->SetSensitiveDetector(scifi1HSD);        
        scifiStation1LogicalVer->SetSensitiveDetector(scifi1VSD);        

        scifiStation2Logical->SetSensitiveDetector(scifi2EpoxySD);
        scifiStation2LogicalHor->SetSensitiveDetector(scifi2HSD);        
        scifiStation2LogicalVer->SetSensitiveDetector(scifi2VSD);        

        scifiStation3Logical->SetSensitiveDetector(scifi3EpoxySD);
        scifiStation3LogicalHor->SetSensitiveDetector(scifi3HSD);        
        scifiStation3LogicalVer->SetSensitiveDetector(scifi3VSD);        

        scifiStation4Logical->SetSensitiveDetector(scifi4EpoxySD);
        scifiStation4LogicalHor->SetSensitiveDetector(scifi4HSD);        
        scifiStation4LogicalVer->SetSensitiveDetector(scifi4VSD);  
    }
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetNISTMaterial(G4String name, G4String fallbackName)
{
    /*
        Get material properties from the NIST database
    */
    
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* material;
    
    try                         // Check if material in database
    {
        G4Material* checkMat = nistManager->FindOrBuildMaterial(name);
        if( checkMat == 0 )
            throw(checkMat);
        material = G4Material::GetMaterial(name);
    }
    catch(G4Material* ex)       // If material name not recgonised revert to default material
    {
        G4ExceptionDescription msg;
        msg << "Could not find material: " << name << ". Reverting material to " << fallbackName << ". Note the name is case sensitive.";
        G4Exception("DetectorConstruction::GetNISTMaterial()","Code001",JustWarning,msg);
        material = G4Material::GetMaterial(fallbackName);
    }
    
    return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldVolume(G4ThreeVector vector)
{
    /*
        Function for UI command to change the world volume        
    */
    
    G4double newX = vector[0];
    G4double newY = vector[1];
    G4double newZ = vector[2];
        
    // Check new volume is bigger than water volume
    // **************************************************************
    if( newX < waterX || newY < waterY || newZ < waterZ )
    {
        G4ExceptionDescription msg;
        msg << "Dimensions specified are smaller than the water volume!\n World volume is left unchanged.\n";
        G4Exception("DetectorConstruction::SetWorldVolume()","Code001",JustWarning,msg);
        return;
    }
    
    // Set new values to world
    // **************************************************************
    worldX = newX;
    worldY = newY;
    worldZ = newZ;
    
    // Modify geometry
    // **************************************************************
    G4RunManagerKernel::GetRunManagerKernel()->DefineWorldVolume(Construct(),true); // To visualise volume change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetPhantomVolume(G4ThreeVector vector)
{
    /*
        Function for UI command to change phantom volume 
    */
    
    G4double newX = vector[0];
    G4double newY = vector[1];
    G4double newZ = vector[2];
        
    // Check new volume is smaller than the world volume
    // **************************************************************
    if( newX > worldX || newY > worldY || newZ > worldZ )
    {
        G4ExceptionDescription msg;
        msg << "Dimensions specified are larger than the world volume!\n Phantom volume is left unchanged.\n";
        G4Exception("DetectorConstruction::SetPhantomVolume()","Code002",JustWarning,msg);
        return;
    }

    // Check whether phantom walls are enabled else it is applied to water volume
    // **************************************************************
    if(!togglePhantom)
    {
        G4ExceptionDescription msg;
        msg << "Note that phantom walls are not enabled! So volume changes are applied to the water volume instead.\n";
        G4Exception("DetectorConstruction::SetPhantomVolume()","Code001",JustWarning,msg);
        
        waterX = newX;
        waterY = newY;
        waterZ = newZ;
    }
    
    // Set new values
    // **************************************************************
    phantomX = newX;
    phantomY = newY;
    phantomZ = newZ;
            
    // Modify geometry
    // **************************************************************
    G4RunManagerKernel::GetRunManagerKernel()->DefineWorldVolume(Construct(),true); // To visualise volume change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetWaterVolume(G4ThreeVector vector)
{
    /*
        Function for UI command to change water volume
    */
    
    G4double newX = vector[0];
    G4double newY = vector[1];
    G4double newZ = vector[2];
    
    // Check new volume is smaller than world volume
    // **************************************************************
    if( newX > worldX || newY > worldY || newZ > worldZ )
    {
        G4ExceptionDescription msg;
        msg << "Dimensions specified are larger than the world volume!\n Water volume is left unchanged.\n";
        G4Exception("DetectorConstruction::SetWaterVolume()","Code001",JustWarning,msg);
        return;
    }
    // Check new volume is smaller than the phantom volume
    // **************************************************************
    else if( (newX > phantomX || newY > phantomY || newZ > phantomZ) && togglePhantom == true)
    {
        G4ExceptionDescription msg;
        msg << "Dimensions specified are larger than the phantom volume!\n Water volume is left unchanged.\n";
        G4Exception("DetectorConstruction::SetWaterVolume()","Code002",JustWarning,msg);
        return;        
    }
    
    // Set new values
    // **************************************************************
    waterX = newX;
    waterY = newY;
    waterZ = newZ;
            
    // Modify geometry
    // **************************************************************
    G4RunManagerKernel::GetRunManagerKernel()->DefineWorldVolume(Construct(),true); // To visualise volume change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetPhantomThickness(G4double val)
{
    /*
        Function for UI command to modify the thickness of phantom walls
    */
    
    G4double newThickness = val;
    
    // Check new thickness is smaller than the phantom dimensions
    // **************************************************************
    if( newThickness > phantomX || newThickness > phantomY || newThickness > phantomZ )
    {
        G4ExceptionDescription msg;
        msg << "Note the thickness specified will cause overlaps!\n";
        G4Exception("DetectorConstruction::SetPhantomThickness()","Code001",JustWarning,msg);
    }
    
    // Check if specified thickness is larger than the entrance window
    // **************************************************************
    if( newThickness <= (2*eWindowZ) )
    {
        G4ExceptionDescription msg;
        msg << "Note the thickness specified is less than or equal to the air gap. \n Due to this, air gap has been resized to specified thickness. \n (i.e. incident beam only encounters G4_AIR and not the phantom material.";
        G4Exception("DetectorConstruction::SetPhantomThickness()","Code002",JustWarning,msg);
        eWindowZ = val/2;
    }
    
    // Set new values
    // **************************************************************
    phantomThickness = val;

    // Modify geometry
    // **************************************************************
    G4RunManagerKernel::GetRunManagerKernel()->DefineWorldVolume(Construct(),true); // To visualise volume change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::EnablePhantomWall(G4bool toggle)
{        
    /*
        Function for UI command to toggle phantom wall
    */
    
    togglePhantom = toggle;
    
    // If true then add walls, else water
    // **************************************************************
    if(togglePhantom == true)
    {
        phantomX = waterX;
        phantomY = waterY;
        phantomZ = waterZ;
    }
    else // If phantom wall is disabled
    {
        waterX = phantomX;
        waterY = phantomY;
        waterZ = phantomZ;
    }
    
    // Modify geometry
    // **************************************************************
    G4RunManagerKernel::GetRunManagerKernel()->DefineWorldVolume(Construct(),true); // To visualise volume change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetWorldMaterial(G4String name)
{
    /*
        Function for UI command to change world material
    */
    worldMaterial = GetNISTMaterial(name, worldMaterial->GetName());
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    G4cout << "Changed world volume material to: " << worldMaterial->GetName() << G4endl;
}

void DetectorConstruction::SetPhantomMaterial(G4String name)
{
    /*
        Function for UI command to change phantom material
    */
    phantomMaterial = GetNISTMaterial(name, phantomMaterial->GetName());
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    G4cout << "Changed phantom volume material to: " << phantomMaterial->GetName() << G4endl;
}

void DetectorConstruction::SetWaterMaterial(G4String name)
{
    /*
        Function for UI command to change water material
    */
    waterMaterial = GetNISTMaterial(name, waterMaterial->GetName());
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
    G4cout << "Changed water volume material to: " << waterMaterial->GetName() << G4endl;
}

void DetectorConstruction::EnableScifi(G4bool toggle)
{        
    /*
        Function for UI command to enable scifi detectors
    */
    toggleScifi = toggle;
        
    G4RunManagerKernel::GetRunManagerKernel()->DefineWorldVolume(Construct(),true); // To visualise volume change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetStation1Pos(G4ThreeVector vector)
{
    /*
        Function for UI command to modify station 1 position
    */
    station1Pos[0] = vector[0];
    station1Pos[1] = vector[1];
    station1Pos[2] = vector[2];
    
    G4RunManagerKernel::GetRunManagerKernel()->DefineWorldVolume(Construct(),true); // To visualise volume change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetStation2Pos(G4ThreeVector vector)
{
    /*
        Function for UI command to modify station 2 position
    */
    station2Pos[0] = vector[0];
    station2Pos[1] = vector[1];
    station2Pos[2] = vector[2];
    
    G4RunManagerKernel::GetRunManagerKernel()->DefineWorldVolume(Construct(),true); // To visualise volume change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetStation3Pos(G4ThreeVector vector)
{
    /*
        Function for UI command to modify station 3 position
    */
    station3Pos[0] = vector[0];
    station3Pos[1] = vector[1];
    station3Pos[2] = vector[2];
    
    G4RunManagerKernel::GetRunManagerKernel()->DefineWorldVolume(Construct(),true); // To visualise volume change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetStation4Pos(G4ThreeVector vector)
{
    /*
        Function for UI command to modify station 4 position
    */
    station4Pos[0] = vector[0];
    station4Pos[1] = vector[1];
    station4Pos[2] = vector[2];
    
    G4RunManagerKernel::GetRunManagerKernel()->DefineWorldVolume(Construct(),true); // To visualise volume change
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Functions to return simulation geometry
// **************************************************************
G4double DetectorConstruction::GetPhantomX()
{
    if(togglePhantom)
        return phantomX;
    else
        return waterX;
}

G4double DetectorConstruction::GetPhantomY()
{
    if(togglePhantom)
        return phantomY;
    else
        return waterY;
}

G4double DetectorConstruction::GetPhantomZ()
{
    if(togglePhantom)
        return phantomZ;
    else
        return waterZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
