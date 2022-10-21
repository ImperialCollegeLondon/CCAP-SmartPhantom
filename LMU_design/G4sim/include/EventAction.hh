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

#ifndef EventAction_h
#define EventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "DetectorConstruction.hh"
#include "SciFiHit.hh"
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

/// Event action
// **************************************************************
struct EventData
{
    TTree* tree;
    G4int id;
    SciFiHitsCollection* hColl;
};

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    virtual ~EventAction();
    
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    void CreateSciFiTree();
    void DeleteSciFiTree();
    void ResetTreeWriteStatus();
    
private:
    G4int eventID;

    G4RunManager* fManager;
    DetectorConstruction* detector;
    
    EventData phantom;
    
    EventData station1Epoxy;
    EventData station1PlaneH;
    EventData station1PlaneV;
    
    EventData station2Epoxy;
    EventData station2PlaneH;
    EventData station2PlaneV;

    EventData station3Epoxy;
    EventData station3PlaneH;
    EventData station3PlaneV;

    EventData station4Epoxy;
    EventData station4PlaneH;
    EventData station4PlaneV;
    
    G4bool wrotePhantomTree;
    G4bool wroteSciFiTree;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
