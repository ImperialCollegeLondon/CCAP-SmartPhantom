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

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Analysis.hh"
#include "Run.hh"
#include "RootIO.hh"
#include <sstream>

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    
  G4cout << "Begin of RunAction" << G4endl;  
  
  PrimaryGeneratorAction* fPGA = (PrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
  G4String outName = fPGA->GetOutputName();
  RootIO::GetInstance(outName);
  
  G4cout << "    ----> RootIO::GetInstance() OK" << G4endl;  
  
  G4cout << "End of RunAction" << G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{    

    const Run* myrun = dynamic_cast<const Run*>(run);
    if ( myrun )
    {
        G4int nEvents = myrun->GetNumberOfEvent();
        if ( nEvents < 1 )
        {
            G4ExceptionDescription msg;
            msg << "Run consists of 0 events";
            G4Exception("RunAction::EndOfRunAction()",
                        "Code001", JustWarning, msg);
        }
        
        G4cout << "Data written to: " << RootIO::GetInstance()->GetCheckedName() << G4endl; 

        PrimaryGeneratorAction* fPGA = (PrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
        fPGA->ResetStartLine();
        
        EventAction* evtAct = (EventAction*)G4RunManager::GetRunManager()->GetUserEventAction();
        evtAct->ResetTreeWriteStatus();
        
        RootIO::GetInstance()->Close();
        
        G4cout << "RootIO File Closed." << '\n' << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun() 
{
    return new Run;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
