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
#include "Randomize.hh"

#include "RootIO.hh"

#include <iostream>
#include <fstream>
#include <limits>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),     
  fParticleGun(0), fProton(0), fMessenger(nullptr)
{
    // Default IO names
    fInputFile = "test_file.dat";
    fOutputFile = "waterhits";
    
    fParticleGun  = new G4ParticleGun(1);
 
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();    
    G4String particleName;

    fProton = particleTable->FindParticle(particleName="proton");
    
    fParticleGun->SetParticleDefinition(fProton);
    
    fManager = G4RunManager::GetRunManager();
    
    DefineCommands();
    
    startLine = 0;
    startTime = 0;
    
    checkedTime = false;
    resetSL = false;
    Debug = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
    std::ifstream sFile;
    
    if (Debug) G4cout << "PrimaryGeneratorAction begins:" << G4endl;
    if (Debug) G4cout << "    ----> Preparing to Read File:" << G4endl;
    
    try{
      if(fInputFile.isNull()) {        
        throw(fInputFile);  
      }
      sFile.open(fInputFile);
    }
    catch (G4String ex)
    {
      G4ExceptionDescription msg;
      msg << "No input file specified.\n";
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries()","Code001",RunMustBeAborted,msg);
      fManager->AbortRun();
      return;
    }
        
    if (Debug) G4cout << "    <---- File " << fInputFile << " read" << G4endl;
    G4String line;
    line.reserve(100);
  
    if(!checkedTime)
        CheckTime(sFile); // Find time offset so all particles arrive >= 0
    
    /////////////////////////////////////////////////////////////////////////////////
    // We want to loop through file, but we only want to shoot one primary event   //
    // at a time. So use variable startLine to skip previous lines.                //
    /////////////////////////////////////////////////////////////////////////////////
    ReadLine(sFile, line);
    AssignToGun(line, event);
    
    sFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Just for convenience to ease iterating the startLine variable
void PrimaryGeneratorAction::GetLine(std::ifstream& readFile, G4String& readLine)
{
    std::getline(readFile, readLine);
    startLine++;
    
    // To handle looping back to start if EOF reached
    if( readFile.peek() == EOF)
    {
        if (Debug) G4cout << "    <---- EOF detected on next line, preparing to loop to start on next iteration." << G4endl;            
        resetSL = true;
    }
}

// To check if line starts with a numeric value, otherwise it is likely a header
G4bool PrimaryGeneratorAction::CheckNumeric(G4String& str)
{
    try
    {
        size_t size;
        std::stod(str, &size);
        if (Debug) G4cout << "    <---- Line read without issue." << G4endl;
        return true;
    }
    catch(const std::invalid_argument&)
    {
        if (Debug) G4cout << "    <---- Detected a possible header line, skipping to next line." << G4endl;
        return false;
    }
    catch(const std::out_of_range&)
    {
        if (Debug) G4cout << "    <---- Out of range of double type." << G4endl;
        return false;        
    }
}

void PrimaryGeneratorAction::SkipNonNumeric(std::ifstream& readFile, G4String& readLine)
{
    while(!CheckNumeric(readLine))
    {
        if (Debug) G4cout << "        ----> Attempting to read line " << startLine+1 << "." << G4endl;
        if(resetSL) ResetStartLine();
        GetLine(readFile, readLine);        
    }
}

// Obtains the line, skipping as necessary via skipLine
void PrimaryGeneratorAction::ReadLine(std::ifstream& readFile, G4String& readLine)
{
    readFile.clear();
    readFile.seekg(0, std::ios::beg);

    // Note line numbers will be reported using 1-based numbering (i.e. very first line is line one not zero)
    if (Debug) G4cout << "    ----> Preparing to loop through lines; startLine = " << startLine << "." << G4endl;
    
    if ( startLine == 0 )    // First line
    {
        if (Debug) G4cout << "        ----> Attempting to read line " << startLine+1 << "." << G4endl;
        GetLine(readFile,readLine);
        SkipNonNumeric(readFile, readLine);
        if (Debug) G4cout << "        <---- First valid line read; line " << startLine << "." << G4endl;
    }
    else    // Other lines
    {
        // Skipping lines
        for(int currLine=0; currLine<startLine; currLine++)
        {
            if (Debug) G4cout << "        ----> Skipping line " << currLine+1 << "." << G4endl;
            readFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
        // Read Line
        if (Debug) G4cout << "        ----> Attempt to read line " << startLine+1 << "." << G4endl;
        GetLine(readFile, readLine);
        SkipNonNumeric(readFile, readLine);
        if (Debug) G4cout << "        <---- Line " << startLine << " read." << G4endl;
    }
    
    if(resetSL) ResetStartLine();
            
    if (Debug) G4cout << "    <---- End of file handling done." << G4endl;
}

// Reads the line splitting into tokens via tab delimited
std::vector<G4double> PrimaryGeneratorAction::GetToken(G4String& newline)
{
    std::istringstream iss(newline);
    std::string token;
    std::vector<G4double> tokens;
    
    while(std::getline(iss,token,' '))
    {
      tokens.push_back(std::stof(token));
    }
    
    return tokens;
}

// Assigns the values in the obtained line to the particle gun
void PrimaryGeneratorAction::AssignToGun(G4String& newline, G4Event* event)
{
    DetectorConstruction* detector = (DetectorConstruction*)fManager->GetUserDetectorConstruction();
    G4double worldZ = detector->GetWorldZ();

    std::vector<G4double> tokens = GetToken(newline);

    G4double x0 = tokens[0]*m;  // Specifying unit as meters
    G4double xp0 = tokens[1];    // Normalised horizontal momentum
    G4double y0 = tokens[2]*m;  // Specifying unit as meters
    G4double yp0 = tokens[3];    // Normalised vertical momentum
    G4double z0 = -worldZ; // Starting at the start of world volume
    G4double zp0 = sqrt(1-pow(xp0,2)-pow(yp0,2));
    G4double energy = tokens[4]*MeV; //Specifying unit as MeV
    G4double dt = (tokens[5] - startTime)*ns; //Specifying unit as nanoseconds
    
    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xp0,yp0,zp0));
    fParticleGun->SetParticleTime(dt);
    fParticleGun->GeneratePrimaryVertex(event);
}

// Scan of the times of input file to find minimum and shift in order that the earliest particle arrives at time = 0
void PrimaryGeneratorAction::CheckTime(std::ifstream& readFile)
{
    if(Debug)
        G4cout << "    ----> Preparing to loop through file to find the earliest time." << G4endl;
    
    readFile.clear();
    readFile.seekg(0, std::ios::beg);
    G4String fileLine;
   
    // Check header    
    while( !readFile.eof() )
    {
        GetLine(readFile, fileLine);
        SkipNonNumeric(readFile, fileLine);
        std::vector<G4double> tokens = GetToken(fileLine);

        if (tokens[5] < startTime)
        {
            startTime = tokens[5];
        }
    }
    
    if(Debug)
        G4cout << "    <---- Earliest time at " << startTime << ", adjusting times so earliest time begins at time = 0" << G4endl;    
    
    checkedTime = true;
    ResetStartLine();
}

void PrimaryGeneratorAction::ResetStartLine()
{
    startLine = 0;
    resetSL = false;
    if(Debug)
        G4cout << "    <---- Looped to start of file for next event or run." << G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// UI commands
void PrimaryGeneratorAction::DefineCommands()
{
    // Define /SP/IO command directory with generic messenger class
    fMessenger = new G4GenericMessenger(this, "/SP/IO/", "I/O control");
    
    // Input file command
    auto& inputFileCmd = fMessenger->DeclareProperty("inputFile", fInputFile);
    G4String guidance = "Specify input file. Note the z position is set to the start of world volume.\n";
    guidance += "Format is: x [m], x' [rad], y [m], y' [rad], E [MeV], t [ns] \n";
    guidance += "(i.e. horizontal position, normalised horizontal momentum, vertical position, normalised vertical momentum, kinetic energy, time)";
    inputFileCmd.SetGuidance(guidance);
    
    // Output file command
    auto& outputFileCmd = fMessenger->DeclareProperty("outputFile", fOutputFile);
    guidance = "Specify output file name.";
    outputFileCmd.SetGuidance(guidance);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
