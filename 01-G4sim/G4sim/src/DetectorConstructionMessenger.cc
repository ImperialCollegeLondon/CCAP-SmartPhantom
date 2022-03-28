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

#include "DetectorConstructionMessenger.hh"

#include "G4UIcommand.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction* detector)
{
    /*
        Define new UI commands relating to DetectorConstruction
    */
    
    m_detector = detector;
    m_detDir = new G4UIdirectory("/SP/detector/");
    m_detDir->SetGuidance("Detector control");
    
    // Modify World Volume 
    // **************************************************************
    m_setWorldVolumeCmd = new G4UIcmdWith3VectorAndUnit("/SP/detector/setWorldVolume",this);
    m_setWorldVolumeCmd->SetGuidance("Modify dimensions for world volume.");
    m_setWorldVolumeCmd->SetParameterName("Half_X","Half_Y","Half_Z",false);
    m_setWorldVolumeCmd->SetRange("Half_X>0. && Half_Y>0 && Half_Z >0");
    m_setWorldVolumeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // World Material
    // **************************************************************
    m_setWorldMaterialCmd = new G4UIcmdWithAString("/SP/detector/setWorldMaterial",this);
    m_setWorldMaterialCmd->SetGuidance("Change material of world volume to a material from NIST database. Note the name is case sensitive.");
    m_setWorldMaterialCmd->SetParameterName("Material",false);
    m_setWorldMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // Modify Phantom Volume
    // **************************************************************
    m_setPhantomVolumeCmd = new G4UIcmdWith3VectorAndUnit("/SP/detector/setPhantomVolume",this);
    m_setPhantomVolumeCmd->SetGuidance("Modify dimensions for overall phantom volume and resize water volume according to specified phantom wall thickness (Must enable phantom walls).");
    m_setPhantomVolumeCmd->SetParameterName("Half_X","Half_Y","Half_Z",false);
    m_setPhantomVolumeCmd->SetRange("Half_X>0. && Half_Y>0 && Half_Z >0");
    m_setPhantomVolumeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    // Set Phantom Wall
    // **************************************************************
    m_enablePhantomCmd = new G4UIcmdWithABool("/SP/detector/enablePhantomWall",this);
    m_enablePhantomCmd->SetGuidance("Enable phantom walls, unless specified it will have a 10 mm thickness. Adopts the dimensions of the water volume, and resize water volume. If set to false, it will scale the dimensions of the water volume to the current phantom dimensions. \n Modify the dimensions with '/SP/detector/setPhantomVolume' command");
    m_enablePhantomCmd->SetParameterName("Bool",false);
    m_enablePhantomCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // Modify Phantom Wall Thickness
    // **************************************************************
    m_setPhantomThicknessCmd = new G4UIcmdWithADoubleAndUnit("/SP/detector/setPhantomThickness",this);
    m_setPhantomThicknessCmd->SetGuidance("Modify thickness of phantom wall, but retain the same phantom dimensions. Water volume will be resized according to the specified thickness. (Must enable phantom walls).");
    m_setPhantomThicknessCmd->SetParameterName("Thickness",false);
    m_setPhantomThicknessCmd->SetRange("Thickness>0.");
    m_setPhantomThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // Phantom Material
    // **************************************************************
    m_setPhantomMaterialCmd = new G4UIcmdWithAString("/SP/detector/setPhantomMaterial",this);
    m_setPhantomMaterialCmd->SetGuidance("Change material of phantom volume to a material from NIST database. Note the name is case sensitive.");
    m_setPhantomMaterialCmd->SetParameterName("Material",false);
    m_setPhantomMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // Modify Water Volume
    // **************************************************************
    m_setWaterVolumeCmd = new G4UIcmdWith3VectorAndUnit("/SP/detector/setWaterVolume",this);
    m_setWaterVolumeCmd->SetGuidance("Modify dimensions for water volume.");
    m_setWaterVolumeCmd->SetParameterName("Half_X","Half_Y","Half_Z",false);    
    m_setWaterVolumeCmd->SetRange("Half_X>0. && Half_Y>0 && Half_Z >0");
    m_setWaterVolumeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
    // Water Material
    // **************************************************************
    m_setWaterMaterialCmd = new G4UIcmdWithAString("/SP/detector/setWaterMaterial",this);
    m_setWaterMaterialCmd->SetGuidance("Change material of water volume to a material from NIST database. Note the name is case sensitive.");
    m_setWaterMaterialCmd->SetParameterName("Material",false);
    m_setWaterMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);    
    
    // Set SciFi Detectors
    // **************************************************************
    m_enableScifiCmd = new G4UIcmdWithABool("/SP/detector/enableScifi",this);
    m_enableScifiCmd->SetGuidance("Enable scintillating fibre detectors.");
    m_enableScifiCmd->SetParameterName("Bool",false);
    m_enableScifiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // Set Station 1 Position
    // **************************************************************
    m_setStation1PosCmd = new G4UIcmdWith3VectorAndUnit("/SP/detector/setStation1Pos",this);
    m_setStation1PosCmd->SetGuidance("Modify position for station 1. (In local coordinates WRT water volume)");
    m_setStation1PosCmd->SetParameterName("Pos_X","Pos_Y","Pos_Z",false);    
    m_setStation1PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // Set Station 2 Position
    // **************************************************************
    m_setStation2PosCmd = new G4UIcmdWith3VectorAndUnit("/SP/detector/setStation2Pos",this);
    m_setStation2PosCmd->SetGuidance("Modify position for station 2. (In local coordinates WRT water volume)");
    m_setStation2PosCmd->SetParameterName("Pos_X","Pos_Y","Pos_Z",false);    
    m_setStation2PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // Set Station 3 Position
    // **************************************************************
    m_setStation3PosCmd = new G4UIcmdWith3VectorAndUnit("/SP/detector/setStation3Pos",this);
    m_setStation3PosCmd->SetGuidance("Modify position for station 3. (In local coordinates WRT water volume)");
    m_setStation3PosCmd->SetParameterName("Pos_X","Pos_Y","Pos_Z",false);    
    m_setStation3PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

    // Set Station 4 Position
    // **************************************************************
    m_setStation4PosCmd = new G4UIcmdWith3VectorAndUnit("/SP/detector/setStation4Pos",this);
    m_setStation4PosCmd->SetGuidance("Modify position for station 4. (In local coordinates WRT water volume)");
    m_setStation4PosCmd->SetParameterName("Pos_X","Pos_Y","Pos_Z",false);    
    m_setStation4PosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstructionMessenger::~DetectorConstructionMessenger()
{
    delete m_detDir;
    delete m_setWorldVolumeCmd;
    delete m_setPhantomVolumeCmd;
    delete m_setWaterVolumeCmd;
    delete m_enablePhantomCmd;
    delete m_setPhantomThicknessCmd;
    delete m_setWorldMaterialCmd;
    delete m_setPhantomMaterialCmd;
    delete m_setWaterMaterialCmd; 
    delete m_enableScifiCmd;
    delete m_setStation1PosCmd;
    delete m_setStation2PosCmd;
    delete m_setStation3PosCmd;
    delete m_setStation4PosCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
    /*
        Pass values supplied to Messenger (from G4UI) to DetectorConstruction functions 
    */
    
    if( command == m_setWorldVolumeCmd )
        m_detector->SetWorldVolume( m_setWorldVolumeCmd->GetNew3VectorValue(newValue) );

    if( command == m_setPhantomVolumeCmd )
        m_detector->SetPhantomVolume( m_setPhantomVolumeCmd->GetNew3VectorValue(newValue) );

    if( command == m_setWaterVolumeCmd )
        m_detector->SetWaterVolume( m_setWaterVolumeCmd->GetNew3VectorValue(newValue) );

    if( command == m_enablePhantomCmd )
        m_detector->EnablePhantomWall( m_enablePhantomCmd->ConvertToBool(newValue) );

    if( command == m_setPhantomThicknessCmd )
        m_detector->SetPhantomThickness( m_setPhantomThicknessCmd->GetNewDoubleValue(newValue) );

    if( command == m_setWorldMaterialCmd )
        m_detector->SetWorldMaterial( newValue );

    if( command == m_setPhantomMaterialCmd )
        m_detector->SetPhantomMaterial( newValue );

    if( command == m_setWaterMaterialCmd )
        m_detector->SetWaterMaterial( newValue );

    if( command == m_enableScifiCmd )
        m_detector->EnableScifi( m_enableScifiCmd->ConvertToBool(newValue) );
    
    if( command == m_setStation1PosCmd )
        m_detector->SetStation1Pos( m_setStation1PosCmd->GetNew3VectorValue(newValue) );

    if( command == m_setStation2PosCmd )
        m_detector->SetStation2Pos( m_setStation2PosCmd->GetNew3VectorValue(newValue) );

    if( command == m_setStation3PosCmd )
        m_detector->SetStation3Pos( m_setStation3PosCmd->GetNew3VectorValue(newValue) );

    if( command == m_setStation4PosCmd )
        m_detector->SetStation4Pos( m_setStation4PosCmd->GetNew3VectorValue(newValue) );    
}
