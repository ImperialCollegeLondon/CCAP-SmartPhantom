/*    User Analysis class implementation file */

#include <iostream>
#include <filesystem>
#include <vector>

#include "RunControl.hpp"
#include "nuAnalysis.hpp"

nuSIMUserAnal::nuSIMUserAnal(bool Dbg) {

  nuAnalysis::setDebug(Dbg);
  if ( nuAnalysis::getDebug() ) {
    std::cout << " nuSIMUserAnal: create instance, start of Debug print:" << std::endl;
  }

}

void nuSIMUserAnal::PreEventLoop( bool Dbg ) {
  
  if ( nuAnalysis::getDebug() ) {
    std::cout << " ----> nuSIMUserAnal: Pre-event loop method entered:"
	      << std::endl;
  }
    
  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving Pre-event loop method." << std::endl;
  }
    
}

void nuSIMUserAnal::EventLoop( bool Dbg ) {

  if ( nuAnalysis::getDebug() ) {
    std::cout << " ----> nuSIMUserAnal: Event loop method entered:" << std::endl;
  }

  // Loop over neutrinos in flux ntuple:
  int nEvt = nuAnalysis::getbeam_ch()->GetEntries();
  if ( nuAnalysis::getDebug() ) {
    std::cout << "     ----> Beam ntuple has "
	      << nEvt << " entries" << std::endl;
  }
    
  for (int i=0 ; i<nEvt ; i++) {
    nuAnalysis::getbeam_ch()->GetEntry(i);
    if ( nuAnalysis::getDebug() and i<10) {
      std::cout << "          ----> Event: " << i << std::endl;
    }
    
  }
  
  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving event loop method." << std::endl;
  }
 
}

void nuSIMUserAnal::PostEventLoop( bool Dbg ) {

  if ( nuAnalysis::getDebug() ) {
    std::cout << " ----> nuSIMUserAnal: Post event loop method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();
  
  if ( nuAnalysis::getDebug() ) {
    std::cout << " <---- Leaving post event loop method." << std::endl;
  }
 
}

int main(int nArgs, char **ArgV){

  // Initialse and parse input arguments:
  RunControl* RC = RunControl::getInstance();
  RC->ParseArgs(nArgs, ArgV);
  if ( RC->getDebug() )
    RC->print();

  nuSIMUserAnal* nuA = new nuSIMUserAnal(RC->getDebug());
  nuA->PreEventLoop();
  nuA->EventLoop();
  nuA->PostEventLoop();
  nuA->HistFitDo();

  }
