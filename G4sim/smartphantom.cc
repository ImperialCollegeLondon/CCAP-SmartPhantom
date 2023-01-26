/*          --------  SmartPhantom G4 simulation main  --------          */

//  G4 header files:
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

//  SmartPhantom header files:
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"


#include "G4UImanager.hh"
#include "G4ScoringManager.hh"
#include "QGSP_BIC.hh"
#include "G4PhysListFactory.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  //Detect interactive mode (if no argument) and define UI session
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) { //No commands line argument
    ui = new G4UIExecutive(argc,argv);
  }

    // Construct the default run manager
    G4RunManager* runManager = new G4RunManager;
    
    // Activate command-based scorer
    G4ScoringManager::GetScoringManager();
    
    // Mandatory user initialization classes
    
    //====================
    //The Geometry
    runManager->SetUserInitialization(new DetectorConstruction);
    
    //====================
    //The Physics
    G4PhysListFactory *physListFactory = new G4PhysListFactory();
    G4VUserPhysicsList *physicsList = physListFactory->GetReferencePhysList("QGSP_BERT");
    runManager->SetUserInitialization(physicsList);
    
    //====================
    // User action initialization
    runManager->SetUserInitialization(new ActionInitialization());
    
    // Visualization manager construction
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    
    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (argc>1) {
        // execute an argument macro file if exist
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }
    else {
      //We have visualization, initialize defaults: look in init_vis.mac macro
      UImanager->ApplyCommand("/control/execute init_vis.mac");
      ui->SessionStart();
      delete ui;
    }
    // Job termination
  
    delete visManager;
    delete runManager;
    
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
