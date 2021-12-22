/*    nuAnalysis class implementation file */

#include <iostream>
#include <filesystem>
#include <vector>

#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"

//#include "TVector.h"
#include "TLorentzVector.h"

#include "RunControl.hpp"
#include "nuAnalysis.hpp"

nuAnalysis::nuAnalysis(bool Dbg) {

  Debug = Dbg;
  if ( Debug ) {
    std::cout << " nuAnalysis: create instance, start of Debug print:" << std::endl;
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
  
  // Initialise files for reading:
  runInfo_ch = new TChain("runInfo");
  beam_ch    = new TChain("beam");
  flux_ch    = new TChain("flux");
  if ( RC->getChainFlag() ){
    for ( const auto & entry :
	    std::filesystem::directory_iterator(RC->getCHAINdirname()) ){
      runInfo_ch->AddFile(entry.path().c_str());
      beam_ch->AddFile(entry.path().c_str());
      flux_ch->AddFile(entry.path().c_str());
    }
  }
  if ( RC->getFileFlag() ){
    runInfo_ch->AddFile(RC->getROOTfilename().c_str());
    beam_ch->AddFile(RC->getROOTfilename().c_str());
    flux_ch->AddFile(RC->getROOTfilename().c_str());
  }
  if ( Debug ) {
    std::cout << "     ----> nuAnalysis: print TChains: " << std::endl;
    std::cout << "           Title: " << runInfo_ch->GetName() << std::endl;
    std::cout << "           Title: " << beam_ch->GetName()    << std::endl;
    std::cout << "           Title: " << flux_ch->GetName()    << std::endl;
  }
  
}

void nuAnalysis::PreEventLoop( bool Dbg ) {
  
  if ( Debug ) {
    std::cout << " ----> nuAnalysis: Pre-event loop method entered:"
	      << std::endl;
  }
  
  if ( Debug ) {
    std::cout << " <---- Leaving Pre-event loop method." << std::endl;
  }
    
}

void nuAnalysis::EventLoop( bool Dbg ) {

  if ( Debug ) {
    std::cout << " ----> nuAnalysis: Event loop method entered:" << std::endl;
  }

  // Loop over neutrinos in flux ntuple:
  int nEvt = beam_ch->GetEntries();
  if ( Debug ) {
    std::cout << "     ----> Beam ntuple has "
	      << nEvt << " entries" << std::endl;
  }
    
  for (int i=0 ; i<nEvt ; i++) {
    beam_ch->GetEntry(i);
    if ( Debug and i<10) {
      std::cout << "          ----> Event: " << i << std::endl;
    }
    
  }
  
  if ( Debug ) {
    std::cout << " <---- Leaving event loop method." << std::endl;
  }
 
}

void nuAnalysis::PostEventLoop( bool Dbg ) {

  if ( Debug ) {
    std::cout << " ----> nuAnalysis: Post event loop method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();
  
  if ( Debug ) {
    std::cout << " <---- Leaving post event loop method." << std::endl;
  }
 
}

void nuAnalysis::HistFitDo() {

  if ( Debug ) {
    std::cout << " ----> nuAnalysis: HisFitDo method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();
  
  std::string OutFile = RC->getOUTPUTdirname() + "nuAnalysis.root";

  TFile oRf(OutFile.c_str(),"RECREATE");
  for (int iH = 0; iH < TH1Flist.size(); iH++)
      TH1Flist[iH]->Write();
  for (int iF = 0; iF < TF1list.size(); iF++)
      TF1list[iF]->Write();
  
  if ( Debug ) {
    std::cout << " <---- Leaving HistFitDo method." << std::endl;
  }
 
}

void nuAnalysis::setrunInfo_ch( TChain* rInfo_ch ) {
  nuAnalysis::runInfo_ch = rInfo_ch;
}

void nuAnalysis::setbeam_ch( TChain* bm_ch ) {
  nuAnalysis::beam_ch = bm_ch;
}

void nuAnalysis::setflux_ch( TChain* flx_ch ) {
  nuAnalysis::flux_ch = flx_ch;
}

void nuAnalysis::setDebug(bool Dbg) {
  nuAnalysis::Debug = Dbg;
}

Double_t nuAnalysis::nueErest(Double_t *xin, Double_t *par) {
  Double_t mmu  = 0.1056583745;
  Double_t s    = par[0];
  Double_t x    = (2.*xin[0]/mmu);
  Double_t scl  = (12.*s);
  Double_t dx   = (2.*par[1]/mmu);
  Double_t y    = scl*x*x*(1.-x)*dx;
  if ( y < 0. or x > 1.) y = 0.; 
  return y;
}

Double_t nuAnalysis::numuErest(Double_t *xin, Double_t *par) {
  Double_t mmu  = 0.1056583745;
  Double_t s    = par[0];
  Double_t x    = (2.*xin[0]/mmu);
  Double_t scl  = (2.*s);
  Double_t dx   = (2.*par[1]/mmu);
  Double_t y    = scl*x*x*(3. - 2.*x)*dx;
  if ( y < 0. or x > 1.) y = 0.; 
  return y;
}

nuSIMtstRestFrame::nuSIMtstRestFrame(bool Dbg) {

  nuAnalysis::setDebug(Dbg);
  if ( nuAnalysis::getDebug() ) {
    std::cout << " nuAnalysis::nuAnalysis: create instance, start of Debug print:" << std::endl;
  }

}

void nuSIMtstRestFrame::PreEventLoop( bool Dbg ) {
  
  if ( nuAnalysis::getDebug() ) {
    std::cout << " ----> nuAnalysis: Pre-event loop method entered:"
	      << std::endl;
  }
  
  // Set up histograms etc. prior to event loop:
  TH1F *hmumass    = new TH1F("hmumass",   "muon mass", 100, 0.105, 0.107);
  nuAnalysis::TH1Flist.push_back(hmumass);
    
  TH1F *hnueErest  = new TH1F("hnueErest", "nue energy in rest frame",
			      100, 0., 0.07);
  nuAnalysis::TH1Flist.push_back(hnueErest);
    
  TH1F *hnumuErest = new TH1F("hnumuErest", "numu energy in rest frame",
			      100, 0., 0.07);
  nuAnalysis::TH1Flist.push_back(hnumuErest);

  int nEvt = nuAnalysis::getbeam_ch()->GetEntries();
  Double_t dEvt = nEvt;
  Double_t dE   = 0.07/100.;

  TF1 *fnueErest  = new TF1("nueE",   nueErest,  0., 0.07, 2);
  fnueErest->SetParameters(dEvt,dE);
  fnueErest->SetLineColor(kRed); fnueErest->SetLineStyle(2);
  nuAnalysis::TF1list.push_back(fnueErest);

  TF1 *fnumuErest = new TF1("nuemu", numuErest, 0., 0.07, 2);
  fnumuErest->SetParameters(dEvt,dE);
  fnumuErest->SetLineColor(kRed); fnumuErest->SetLineStyle(2);
  nuAnalysis::TF1list.push_back(fnumuErest);

  if ( nuAnalysis::getDebug() ) {
    std::cout << "     ----> Booked "
	      << TH1Flist.size() << " histos:" << std::endl;
    for(int iH = 0; iH < TH1Flist.size(); iH++) {
      std::cout << "           Printing histo " << iH
		<< ": Title: " << TH1Flist[iH]->GetTitle() << std::endl;
    }
    std::cout << "     ----> Booked "
	      << TF1list.size() << " fits:" << std::endl;
    for(int iF = 0; iF < TF1list.size(); iF++) {
      std::cout << "           Printing fit " << iF
		<< ": Title: " << TF1list[iF]->GetTitle() << std::endl;
    }
  }
  
  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving Pre-event loop method." << std::endl;
  }
    
}

void nuSIMtstRestFrame::EventLoop( bool Dbg ) {

  if ( nuAnalysis::getDebug() ) {
    std::cout << " ----> nuAnalysis: Event loop method entered:" << std::endl;
  }

  // Loop over neutrinos in flux ntuple:
  int nEvt = nuAnalysis::getbeam_ch()->GetEntries();
  if ( nuAnalysis::getDebug() ) {
    std::cout << "     ----> Beam ntuple has "
	      << nEvt << " entries" << std::endl;
  }

  Float_t nueP[4], numuP[4], muP[4];
  nuAnalysis::getbeam_ch()->SetBranchAddress("Nue",  &nueP);
  nuAnalysis::getbeam_ch()->SetBranchAddress("NuMu", &numuP);
  nuAnalysis::getbeam_ch()->SetBranchAddress("muon", &muP);

  TH1F *hmumass    = nuAnalysis::TH1Flist[0];
  TH1F *hnueErest  = nuAnalysis::TH1Flist[1];
  TH1F *hnumuErest = nuAnalysis::TH1Flist[2];

  TLorentzVector Pmu;
  TLorentzVector PmuRest;
  TLorentzVector Pnue, Pnumu;
  TLorentzVector PnueRest;
  TLorentzVector PnumuRest;
  double Mmu;
  TVector3 b;
    
  for (int i=0 ; i<nEvt ; i++) {
    nuAnalysis::getbeam_ch()->GetEntry(i);
    if ( nuAnalysis::getDebug() and i<10) {
      std::cout << "          ----> Event: " << i << std::endl;
      std::cout << "              ----> Nue.P: ";
      for(int ii=0;ii<4;ii++){
	std::cout << nueP[ii] << ", ";
      }
      std::cout << std::endl;
      std::cout << "              ----> NueMu.P: ";
      for(int ii=0;ii<4;ii++){
	std::cout << numuP[ii] << ", ";
      }
      std::cout << std::endl;
      std::cout << "              ----> muon.P: ";
      for(int ii=0;ii<4;ii++){
	std::cout << muP[ii] << ", ";
      }
      std::cout << std::endl;
    }

    Pmu.SetPxPyPzE(muP[1], muP[2], muP[3], muP[0]);
    Mmu = sqrt(Pmu*Pmu);
    if ( nuAnalysis::getDebug() and i<10) {
      std::cout << "                  ----> Mass of muon check: "
		<< Mmu << std::endl;
    }
    hmumass->Fill(Mmu);
    
    b = -Pmu.BoostVector();
    if ( nuAnalysis::getDebug() and i<10) {
      std::cout <<
	"                  ----> Setting up boost parameters check: "
		<< std::endl;
      std::cout << "                      Boost vector: ";
      for(int ii=0;ii<3;ii++){
	std::cout << b[ii] << ", ";
      }
      std::cout << std::endl;
    }
    
    PmuRest = Pmu;
    PmuRest.Boost(b);
    if ( nuAnalysis::getDebug() and i<10) {
      std::cout <<
	"                  ----> Muon 4-mmtm in rest frame check: ";
      for(int ii=0;ii<4;ii++){
	std::cout << PmuRest[ii] << ", ";
      }
      std::cout << std::endl;
    }
  
    Pnue.SetPxPyPzE( nueP[1],  nueP[2],  nueP[3],  nueP[0]);
    Pnumu.SetPxPyPzE(numuP[1], numuP[2], numuP[3], numuP[0]);
    PnueRest  = Pnue;
    PnumuRest = Pnumu;
    PnueRest.Boost(b);
    PnumuRest.Boost(b);

    hnueErest->Fill(PnueRest[3]);
    hnumuErest->Fill(PnumuRest[3]);
    
  }
  
  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving event loop method." << std::endl;
  }
 
}

void nuSIMtstRestFrame::PostEventLoop( bool Dbg ) {

  if ( nuAnalysis::getDebug() ) {
    std::cout << " ----> nuAnalysis: Post event loop method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();

  TH1F *hmumass    = nuAnalysis::TH1Flist[0];
  TH1F *hnueErest  = nuAnalysis::TH1Flist[1];
  TH1F *hnumuErest = nuAnalysis::TH1Flist[2];

  TF1 *fnueErest  = nuAnalysis::TF1list[0];
  TF1 *fnumuErest = nuAnalysis::TF1list[1];

  TCanvas *c = new TCanvas();
  gErrorIgnoreLevel = kWarning;
  
  std::string PltFile = RC->getOUTPUTdirname() + "hnueErest.png";
  hnueErest->Draw();
  fnueErest->Draw("SAME");
  c->Print(PltFile.c_str());
  
  PltFile = RC->getOUTPUTdirname() + "hnumuErest.png";
  hnumuErest->Draw();
  fnumuErest->Draw("SAME");
  c->Print(PltFile.c_str());
  
  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving post event loop method." << std::endl;
  }
 
}
