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

#include "RootIO.hh"
//
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

static RootIO* instance = 0;
G4String RootIO::fOutputFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootIO::RootIO():fNevents(0)
{    
    /*
        Sets up ROOT file
    */
    
    // Initialize ROOT
    // **************************************************************
    if (instance == 0 )
    {
        G4String rootFileName = GetRootName();
        
        fFile = new TFile(rootFileName.c_str(),"RECREATE");

        // Store some parameters of phantom geometry
        // **************************************************************
        TTree* phantomGeo = new TTree("Geometry","Geometry");
        double geomData[3];
        phantomGeo->Branch("Size", &geomData, "Width/D:Height/D:Depth/D");
        
        // Get phantom geometry from DetectorConstruction
        // **************************************************************
        G4RunManager* fManager = G4RunManager::GetRunManager();
        DetectorConstruction* detector = (DetectorConstruction*)fManager->GetUserDetectorConstruction();
        
        geomData[0] = 2*(detector->GetPhantomX());          // Total width
        geomData[1] = 2*(detector->GetPhantomY());          // Total height
        geomData[2] = 2*(detector->GetPhantomZ());          // Total depth
        phantomGeo->Fill();        
    }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootIO::~RootIO()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RootIO* RootIO::GetInstance(G4String name)
{
    /*
        Retrieve current instance or create a new one
        (I don't remember why I did it this way...)
    */
    
    if (instance == 0 )
    {
        G4cout << "Creating new RootIO instance." << G4endl;
        fOutputFile = name;
        instance = new RootIO();
    }
    return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TTree* RootIO::CreateTree(TString treeName)
{
    /*
        Create new TTree in ROOT file with the preset branches:     
                Edep -- Energy deposition
                HitTime -- Global time of hit event (since start)
                EventN -- Which primary event hit originates from
                PosX -- Horizontal position
                PosY -- Vertical position
                PosZ -- Longitudinal position
                Steplength -- Steplength
                DeltaT -- Duration of hit event
    */
    
    TTree* tree = new TTree(treeName,treeName);
    tree->Branch("Data", &data, "Edep/D:HitTime/D:EventN/D:PosX/D:PosY/D:PosZ/D:Steplength/D:DeltaT/D"); // Branch to store numerical data
    tree->Branch("Name", &name);                                                                         // Branch to store particle names
    return tree;
}

void RootIO::WriteToRoot(SciFiHitsCollection* hsf, double evtID, TTree* &tree)
{        
    /*
        Function to write fill TTree with hits
    */
    
    int n_hit = hsf->entries();

    for(int i=0;i<n_hit;i++){
        SciFiHit* hit = (*hsf)[i];
        double eDep = hit->GetEdep();
        if(eDep>0){
            
            data[0] = eDep;
            data[1] = hit->GetTime();
            data[2] = evtID;
            data[3] = hit->GetPos().x();
            data[4] = hit->GetPos().y();
            data[5] = hit->GetPos().z();
            data[6] = hit->GetStepLength();
            data[7] = hit->GetDeltaT();
            name = hit->GetName();
            
            tree->Fill();
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RootIO::Close()
{
    /*
        Close file and reset instance
    */
    
    fFile->Write();
    fFile->Close();
    instance = 0;
    delete fFile;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4String> RootIO::GetFileExt(const G4String& str)
{
    /*
        Check for file extension if any to output ROOT file
    */
    
    std::vector<G4String> file;
    size_t sz = str.rfind('.', str.length());
    if(sz != std::string::npos)
    {
        G4String filename = str.substr(0, sz);
        G4String ext = str.substr(sz+1, str.length() -1);
        file.push_back(filename);
        file.push_back(ext);
    }
    else
    {
        file.push_back(str);
        file.push_back("");
    }
    return file;
}

G4String RootIO::GetRootName()
{
    /*
        Gets the specified name of output file and checks if it has a ".root" extension,
        if not it adds it to the end of filename
    */

    PrimaryGeneratorAction* fPGA = (PrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
    fOutputFile = fPGA->GetOutputName();

    // Check name specified
    // **************************************************************
    if(fOutputFile == "")   // If empty 
    {  
        fOutputFile = "waterhits.root";
        G4cout << "Defaulting to waterhits.root" << G4endl;
    }
    else
    {
        std::vector<G4String> splitFileName = GetFileExt(fOutputFile);
        G4String fName = splitFileName[0];        // Filename
        G4String ext = splitFileName[1];          // Extension
        
        // Check if filename has root extension add one if needed
        // **************************************************************
        if( ext == "" || ext == "root") 
        {
            checkedName = fName + ".root";
        }
        else // If there is already an extension, don't add root extension and leave name as is, 
             // but supply a message to inform user
        {
            checkedName = fName + "." + ext;
            G4cout << "Please note the filename extension is not .root" << G4endl;
        }
    }
    return checkedName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
