/*    Analysis class implementation file */

#include <iostream>
//#include <filesystem>
#include <vector>
#include "dirent.h"

#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"

//#include "TVector.h"
#include "TLorentzVector.h"

#include "RunControl.hpp"
#include "Analysis.hpp"

Analysis::Analysis(bool Dbg) {

  Debug = Dbg;
  if ( Debug ) {
    std::cout << " Analysis: create instance, start of Debug print:" << std::endl;
  }

  RunControl* RC = RunControl::getInstance();
  if ( Debug ) {
    std::cout << "     ----> RunControl: instance with parameters:" << std::endl;
    std::cout << "           Debug         : "
	      << RC->getDebug() << std::endl;
    std::cout << "           ROOT file name: "
	      << RC->getROOTfilename() << std::endl;
    std::cout << "           Directory with nTuples to chain: "
	      << RC->getCHAINdirname() << std::endl;
    std::cout << "           Output directory: "
	      << RC->getOUTPUTdirname() << std::endl;
  }
  
  // Initialise files for reading:  Event
  // Event1_ch = new TChain("Event;1");
  Event_ch  = new TChain("Event");
  if ( RC->getChainFlag() ){
    struct dirent *entry = nullptr;
    DIR *dp = nullptr;
    dp = opendir(RC->getCHAINdirname().c_str());
    if (dp != nullptr) {
      while ((entry = readdir(dp)))
	printf ("%s\n", entry->d_name);
    }
    closedir(dp);
  }
  if ( RC->getFileFlag() ){
    //Event1_ch->AddFile(RC->getROOTfilename().c_str());
    Event_ch->AddFile(RC->getROOTfilename().c_str());
  }
  if ( Debug ) {
    std::cout << "     ----> Analysis: print TChains: " << std::endl;
    // std::cout << "           Title: " << Event1_ch->GetName() << std::endl;
    std::cout << "           Title: " << Event_ch->GetName()    << std::endl;
  }
  //Event1_ch->Print();
  Event_ch->Print();
  
}

void Analysis::PreEventLoop( bool Dbg ) {
  
  if ( Debug ) {
    std::cout << " ----> Analysis: Pre-event loop method entered:"
	      << std::endl;
  }
  
  if ( Debug ) {
    std::cout << " <---- Leaving Pre-event loop method." << std::endl;
  }
    
}

void Analysis::EventLoop( bool Dbg ) {

  if ( Debug ) {
    std::cout << " ----> Analysis: Event loop method entered:" << std::endl;
  }

  // Loop over entries in ntuple:
  int nEvt = Event_ch->GetEntries();
  if ( Debug ) {
    std::cout << "     ----> Beam ntuple has "
  	      << nEvt << " entries" << std::endl;
  }
    
  for (int i=0 ; i<nEvt ; i++) {
    Event_ch->GetEntry(i);
    if ( Debug and i<10) {
      std::cout << "          ----> Event: " << i << std::endl;
    }
    
  }
  
  if ( Debug ) {
    std::cout << " <---- Leaving event loop method." << std::endl;
  }
 
}

void Analysis::PostEventLoop( bool Dbg ) {

  if ( Debug ) {
    std::cout << " ----> Analysis: Post event loop method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();
  
  if ( Debug ) {
    std::cout << " <---- Leaving post event loop method." << std::endl;
  }
 
}

void Analysis::HistFitDo() {

  if ( Debug ) {
    std::cout << " ----> Analysis: HisFitDo method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();
  
  std::string OutFile = RC->getOUTPUTdirname() + "Analysis.root";

  TFile oRf(OutFile.c_str(),"RECREATE");
  for (int iH = 0; iH < TH1Flist.size(); iH++)
      TH1Flist[iH]->Write();
  for (int iF = 0; iF < TF1list.size(); iF++)
      TF1list[iF]->Write();
  
  if ( Debug ) {
    std::cout << " <---- Leaving HistFitDo method." << std::endl;
  }
 
}

void Analysis::setDebug(bool Dbg) {
  Analysis::Debug = Dbg;
}

BraggPeak::BraggPeak(bool Dbg) {

  Analysis::setDebug(Dbg);
  if ( Analysis::getDebug() ) {
    std::cout << " BraggPeak::BraggPeak: create instance, start of" \
	      << " Debug print:" << std::endl;
  }

}

void BraggPeak::PreEventLoop( bool Dbg ) {
  
  if ( BraggPeak::getDebug() ) {
    std::cout << " ----> BraggPeak: Pre-event loop method entered:"
	      << std::endl;
  }
  
  // Set up histograms etc. prior to event loop:
  TH1F *Bragg = new TH1I("Bragg", "Bragg curve", 100, 0., 3.);
  nuAnalysis::TH1Flist.push_back(Bragg);
    
  if ( BraggPeak::getDebug() ) {
    std::cout << " <---- Leaving Pre-event loop method." << std::endl;
  }
    
}

void BraggPeak::EventLoop( bool Dbg ) {

  if ( BraggPeak::getDebug() ) {
    std::cout << " ----> nuAnalysis: Event loop method entered:" << std::endl;
  }

  // Loop over neutrinos in flux ntuple:
  int nEvt = Analysis::getEvent_ch()->GetEntries();
  if ( Analysis::getDebug() ) {
    std::cout << "     ----> Beam ntuple has "
	      << nEvt << " entries" << std::endl;
  }

  Double_t Edep, HitTime, EventN, PosX, PosY, PosZ, StepLength, DeltaT;
  Analysis::getEvent_ch()->SetBranchAddress("Edep",        &Edep);
  Analysis::getEvent_ch()->SetBranchAddress("HitTime",     &HitTime);
  Analysis::getEvent_ch()->SetBranchAddress("EventN",      &EventN);
  Analysis::getEvent_ch()->SetBranchAddress("PosX",        &PosX);
  Analysis::getEvent_ch()->SetBranchAddress("PosY",        &PosY);
  Analysis::getEvent_ch()->SetBranchAddress("PosZ",        &PosZ);
  Analysis::getEvent_ch()->SetBranchAddress("StepLength",  &StepLength);
  Analysis::getEvent_ch()->SetBranchAddress("DeltaT",      &DeltaT);

  TH1F *Brag = Analysis::TH1Flist[0];

  Double_t z = 0.;
  for (int i=0 ; i<nEvt ; i++) {
    Analysis::getEvent_ch()->GetEntry(i);
    if ( nuAnalysis::getDebug() and i<10) {
      std::cout << "          ----> Event: " << i << std::endl;
    }

    z = PosZ + 150.;
    hmumass->Fill(z, Edep);
    
  }
  
  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving event loop method." << std::endl;
  }
 
}

void BraggPeak::PostEventLoop( bool Dbg ) {

  if ( nuAnalysis::getDebug() ) {
    std::cout << " ----> BraggPeak: Post event loop method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();

  TH1F *Bragg    = Analysis::TH1Flist[0];

  TCanvas *c = new TCanvas();
  gErrorIgnoreLevel = kWarning;
  
  std::string PltFile = RC->getOUTPUTdirname() + "BraggPeak.png";
  Bragg->Draw();
  c->Print(PltFile.c_str());
  
  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving post event loop method." << std::endl;
  }
 
}
