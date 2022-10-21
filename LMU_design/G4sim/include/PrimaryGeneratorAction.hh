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

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4RunManager.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "globals.hh"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <cmath>
#include <numeric>
#include <fstream>

class G4ParticleGun;
class G4Event;
class G4ParticleDefinition;
class G4GenericMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();
    virtual void GeneratePrimaries(G4Event*);
    
    G4bool CheckNumeric(G4String& str);
    void SkipNonNumeric(std::ifstream& readFile, G4String& readLine);
    void GetLine(std::ifstream& readFile, G4String& readLine);
    void ReadLine(std::ifstream& readFile, G4String& readLine);
    std::vector<G4double> GetToken(G4String& line);
    void AssignToGun(G4String& line, G4Event* event);
    void ResetStartLine();
    void CheckTime(std::ifstream& readFile);
    
    G4String GetOutputName() { return fOutputFile; }
    
private:
    void DefineCommands();
    
    G4ParticleGun* fParticleGun;
    G4ParticleDefinition* fProton;
    G4RunManager* fManager;
    G4GenericMessenger* fMessenger;
    G4String fInputFile;
    G4String fOutputFile;
    
    G4int startLine;
    G4bool resetSL;
    G4double startTime;
    G4bool checkedTime;    
    G4bool Debug;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
