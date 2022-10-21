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

#include "EventAction.hh"
#include "Analysis.hh"
#include "RootIO.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
  
#include <sstream>
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
: G4UserEventAction()
{
    fManager = G4RunManager::GetRunManager();
    detector = (DetectorConstruction*)fManager->GetUserDetectorConstruction();
    ResetTreeWriteStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::ResetTreeWriteStatus()
{
    /*
        Reset booleans controlling on whether to write phantom/scifi TTree in file
    */
    wrotePhantomTree = false;
    wroteSciFiTree = false;
}

void EventAction::CreateSciFiTree()
{
    /*
        Create the TTree related to SciFi detectors to file and assign to EventData struct
    */
    
    station1Epoxy.tree = RootIO::GetInstance()->CreateTree("Station1_Epoxy");
    station1PlaneH.tree = RootIO::GetInstance()->CreateTree("Station1_PlaneH");
    station1PlaneV.tree = RootIO::GetInstance()->CreateTree("Station1_PlaneV");

    station2Epoxy.tree = RootIO::GetInstance()->CreateTree("Station2_Epoxy");
    station2PlaneH.tree = RootIO::GetInstance()->CreateTree("Station2_PlaneH");
    station2PlaneV.tree = RootIO::GetInstance()->CreateTree("Station2_PlaneV");

    station3Epoxy.tree = RootIO::GetInstance()->CreateTree("Station3_Epoxy");
    station3PlaneH.tree = RootIO::GetInstance()->CreateTree("Station3_PlaneH");
    station3PlaneV.tree = RootIO::GetInstance()->CreateTree("Station3_PlaneV");

    station4Epoxy.tree = RootIO::GetInstance()->CreateTree("Station4_Epoxy");
    station4PlaneH.tree = RootIO::GetInstance()->CreateTree("Station4_PlaneH");
    station4PlaneV.tree = RootIO::GetInstance()->CreateTree("Station4_PlaneV");
}

void EventAction::DeleteSciFiTree()
{
    /*
        Delete TTree from file
    */
    
    gDirectory->Delete("Station1_Epoxy");
    gDirectory->Delete("Station1_PlaneH");
    gDirectory->Delete("Station1_PlaneV");
    gDirectory->Delete("Station2_Epoxy");
    gDirectory->Delete("Station2_PlaneH");
    gDirectory->Delete("Station2_PlaneV");
    gDirectory->Delete("Station3_Epoxy");
    gDirectory->Delete("Station3_PlaneH");
    gDirectory->Delete("Station3_PlaneV");
    gDirectory->Delete("Station4_Epoxy");
    gDirectory->Delete("Station4_PlaneH");
    gDirectory->Delete("Station4_PlaneV");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
    /*
        At the start of event collect the CollectionID for each sensitive detector
    */
    
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
 
    if( !wrotePhantomTree )
    {
        phantom.tree = RootIO::GetInstance()->CreateTree("Phantom"); // Assign to Phantom EventData struct
        wrotePhantomTree = true;
    }
    
    phantom.id = sdManager->GetCollectionID("phantomSD/SciFiColl");
    
    if( detector->GetToggleScifi() )
    {
        if( !wroteSciFiTree )   // Write SciFi Related Trees to Root File
        {
            CreateSciFiTree();
            wroteSciFiTree = true;
        }
        
        station1Epoxy.id = sdManager->GetCollectionID("scifi1EpoxySD/SciFiColl");
        station1PlaneH.id = sdManager->GetCollectionID("scifi1HSD/SciFiColl");
        station1PlaneV.id = sdManager->GetCollectionID("scifi1VSD/SciFiColl");

        station2Epoxy.id = sdManager->GetCollectionID("scifi2EpoxySD/SciFiColl");
        station2PlaneH.id = sdManager->GetCollectionID("scifi2HSD/SciFiColl");
        station2PlaneV.id = sdManager->GetCollectionID("scifi2VSD/SciFiColl");

        station3Epoxy.id = sdManager->GetCollectionID("scifi3EpoxySD/SciFiColl");
        station3PlaneH.id = sdManager->GetCollectionID("scifi3HSD/SciFiColl");
        station3PlaneV.id = sdManager->GetCollectionID("scifi3VSD/SciFiColl");

        station4Epoxy.id = sdManager->GetCollectionID("scifi4EpoxySD/SciFiColl");
        station4PlaneH.id = sdManager->GetCollectionID("scifi4HSD/SciFiColl");
        station4PlaneV.id = sdManager->GetCollectionID("scifi4VSD/SciFiColl");        
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{    
    /*
        At end of event get the collection of hits then write to ROOT file
    */
    
    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce)
    {
        G4ExceptionDescription msg;
        msg << "No hits collection of this event found.\n";
        G4Exception("EventAction::EndOfEventAction()","Code001",JustWarning,msg);
        return;
    }
    
    // Get hits collections 
    // **************************************************************
    phantom.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(phantom.id));

    if( detector->GetToggleScifi() )
    {
        station1Epoxy.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station1Epoxy.id));
        station1PlaneH.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station1PlaneH.id));
        station1PlaneV.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station1PlaneV.id));

        station2Epoxy.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station2Epoxy.id));
        station2PlaneH.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station2PlaneH.id));
        station2PlaneV.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station2PlaneV.id));

        station3Epoxy.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station3Epoxy.id));
        station3PlaneH.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station3PlaneH.id));
        station3PlaneV.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station3PlaneV.id));

        station4Epoxy.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station4Epoxy.id));
        station4PlaneH.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station4PlaneH.id));
        station4PlaneV.hColl = static_cast<SciFiHitsCollection*>(hce->GetHC(station4PlaneV.id));
    }

    // Write to Root file for each TTree
    // **************************************************************
    eventID = event->GetEventID();
    RootIO* instance = RootIO::GetInstance();
    instance->WriteToRoot(phantom.hColl,eventID,phantom.tree);

    if( detector->GetToggleScifi() )
    {
        instance->WriteToRoot(station1Epoxy.hColl,eventID,station1Epoxy.tree);
        instance->WriteToRoot(station1PlaneH.hColl,eventID,station1PlaneH.tree);
        instance->WriteToRoot(station1PlaneV.hColl,eventID,station1PlaneV.tree);
        
        instance->WriteToRoot(station2Epoxy.hColl,eventID,station2Epoxy.tree);
        instance->WriteToRoot(station2PlaneH.hColl,eventID,station2PlaneH.tree);
        instance->WriteToRoot(station2PlaneV.hColl,eventID,station2PlaneV.tree);

        instance->WriteToRoot(station3Epoxy.hColl,eventID,station3Epoxy.tree);
        instance->WriteToRoot(station3PlaneH.hColl,eventID,station3PlaneH.tree);
        instance->WriteToRoot(station3PlaneV.hColl,eventID,station3PlaneV.tree);

        instance->WriteToRoot(station4Epoxy.hColl,eventID,station4Epoxy.tree);
        instance->WriteToRoot(station4PlaneH.hColl,eventID,station4PlaneH.tree);
        instance->WriteToRoot(station4PlaneV.hColl,eventID,station4PlaneV.tree);        
    }    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
