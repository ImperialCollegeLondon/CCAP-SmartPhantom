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

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

#include "RootIO.hh"

#include <iostream>
#include <fstream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),     
  fParticleGun(0), fProton(0)
{
    fParticleGun  = new G4ParticleGun(1);
 
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();    
    G4String particleName;

    fProton = particleTable->FindParticle(particleName="proton");
    
    fParticleGun->SetParticleDefinition(fProton);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{

    //G4cout << "Preparing to Read File." << G4endl;
    std::ifstream readfile;
    readfile.open("test_file.dat");
    //G4cout << "Read File." << G4endl;
    std::string line;
    line.reserve(100);
  
    //lineskip = 1;
    // skip N lines
    for(int i = 0; i<lineskip; i++)
    {
      std::getline(readfile,line);
    }

    std::getline(readfile,line);
  
    // Loop back to start if EOF
    if(!readfile)
    {
      lineskip=0;
      readfile.clear();
      readfile.seekg(0, std::ios::beg);
      std::getline(readfile,line);
    }
  
    std::istringstream iss(line);
    std::string token;
    std::vector<float> tokens;
  
    while(std::getline(iss,token,' '))
    {
      tokens.push_back(std::stof(token));
    }
    float x0 = tokens[0]*m;  // Specifying unit as meters
    float px = tokens[1];    // Normalised horizontal momentum
    float y0 = tokens[2]*m;  // Specifying unit as meters
    float py = tokens[3];    // Normalised vertical momentum
    float z0 = -300*mm/2;// Starting at start of world volume
    float energy = tokens[4]*MeV; //Specifying unit as MeV
    float dt = tokens[5]*ns; //Specifying unit as nanoseconds
    float pz = sqrt(1-pow(px,2)-pow(py,2));
        
/*
    G4double Ekin = 207; // Specifyig kinetic energy 
     
    // Introducing Energy Spread
    G4double randomNumber = G4UniformRand();
    G4double eSpread = 0;
    if (randomNumber > 0.5)
    {
      eSpread = Ekin*(1 + 0.005*G4UniformRand());
    }
    else
    {
      eSpread = Ekin*(1 - 0.005*G4UniformRand());   
    }
*/
    
    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
    fParticleGun->SetParticleTime(dt);
    fParticleGun->GeneratePrimaryVertex(event);
    lineskip++;
    
    readfile.close();
    //G4cout << "File Closed." << G4endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
