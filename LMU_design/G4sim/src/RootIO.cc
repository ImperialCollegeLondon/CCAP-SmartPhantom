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
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TBenchmark.h"

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
                
                fibreN -- Fibre number of event in scifi
                LightYield -- light yield from event in scifi
    */
    
    TTree* tree = new TTree(treeName,treeName);
    tree->Branch("Data", &data, "Edep/D:HitTime/D:EventN/D:PosX/D:PosY/D:PosZ/D:Steplength/D:DeltaT/D:FibreN/D:LightYield/D"); // Branch to store numerical data
    tree->Branch("Name", &name); 
    //")                                                                        // Branch to store particle names
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
            data[8] = hit->GetFibreN();
            data[9] = hit->GetLightYield();
            name = hit->GetName();
            
            tree->Fill();
        }
    }

}

void RootIO::GetGraph(SciFiHitsCollection* hsf){
    int n_hit = hsf->entries();

    int fibs[33] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33};
    TCanvas *c1 = new TCanvas("c1", "Light Yield in Fibres", 10,10,900,500);
    
    gBenchmark->Start("canvas");

    c1->SetGrid();
    c1->SetTopMargin(0.15);
    TH1F *h = new TH1F("h","test",3,0,3);
    h->SetStats(0);
    h->SetFillColor(38);
    h->SetCanExtend(TH1::kAllAxes);

    for (Int_t i=0;i<n_hit;i++) {
        SciFiHit* hit = (*hsf)[i];
        Int_t fibre = hit->GetFibreN();
        h->Fill(fibs[fibre]);
    }
    h->LabelsDeflate();
    h->Draw();
    TObject *c2 = TObject::Clone(c1);

    gBenchmark->Show("canvas");

}




// void RootIO::GetGraph(){
// 
//    int array1[] = {3, 2, 1, 4}; //array of all nums with same fibre number
//    int array2[] = {3, 2, 1, 4}; //array of all nums with same fibre number
//    int array2[] = {3, 2, 1, 4}; //array of all nums with same fibre number
//    int array2[] = {3, 2, 1, 4}; //array of all nums with same fibre number
//    int array2[] = {3, 2, 1, 4}; //array of all nums with same fibre number
//    int array2[] = {3, 2, 1, 4}; //array of all nums with same fibre number
//    int array2[] = {3, 2, 1, 4}; //array of all nums with same fibre number
// 
//    int arrayFibres[] = {0,1,2,3,4,5,6,7};
//    int sum = 0;
// 
//    num = 0;
//    while num <= 33:
//        if hit->GetFibreN() = num
//        
// 
//        num = num + 1
//    h->Draw();
// 


//void RootIO::GetGraph(TString treeName)
//{
//
//
//    //TTree *waterBox;
//    fFile->GetObject("tree", waterBox);
//    TCanvas *myCanvas = new TCanvas();
//    myCanvas->Divide(2,2);
//    myCanvas->cd(1);
//    waterBox->Draw("fibreN:lightYield");
// }
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

    //void make_graph(){
        //const int n_points= sizeof(eDep);
        //TGraphErrors graph(n_points);
        //double y_vals[n_points] = lightYield;
        //double x_vals[n_points] = fibreN;

        //TGraph graph(n_points, x_vals, y_vals);
        //graph.SetTitle("Light yield per fibre");

        //auto mycanvas = newCanvas();

        //graph.DrawClone("APE");

        //mycanvas->Print("graph.pdf")
    //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
