/*    User Analysis class implementation file */

#include <iostream>
//#include <filesystem>
#include <vector>

#include "RunControl.hpp"
#include "Analysis.hpp"

UserAnal::UserAnal(bool Dbg) {

  Analysis::setDebug(Dbg);
  if ( Analysis::getDebug() ) {
    std::cout << " UserAnal: create instance, start of Debug print:" << std::endl;
  }

}

void UserAnal::PreEventLoop( bool Dbg ) {
  
  if ( Analysis::getDebug() ) {
    std::cout << " ----> UserAnal: Pre-event loop method entered:"
	      << std::endl;
  }
    
  if ( Analysis::getDebug() ) {
    std::cout << " <---- Leaving Pre-event loop method." << std::endl;
  }
    
}

void UserAnal::EventLoop( bool Dbg ) {

  if ( Analysis::getDebug() ) {
    std::cout << " ----> UserAnal: Event loop method entered:" << std::endl;
  }

  // Loop over neutrinos in flux ntuple:
  int nEvt = Analysis::getEvent_ch()->GetEntries();
  if ( Analysis::getDebug() ) {
    std::cout << "     ----> Event ntuple has "
	      << nEvt << " entries" << std::endl;
  }
    
  for (int i=0 ; i<nEvt ; i++) {
    Analysis::getEvent_ch()->GetEntry(i);
    if ( Analysis::getDebug() and i<10) {
      std::cout << "          ----> Event: " << i << std::endl;
    }
    
  }
  
  if ( Analysis::getDebug() ) {
    std::cout << " <---- Leaving event loop method." << std::endl;
  }
 
}

void UserAnal::PostEventLoop( bool Dbg ) {

  if ( Analysis::getDebug() ) {
    std::cout << " ----> UserAnal: Post event loop method entered:"
	      << std::endl;
  }
  RunControl* RC = RunControl::getInstance();
  
  if ( Analysis::getDebug() ) {
    std::cout << " <---- Leaving post event loop method." << std::endl;
  }
 
}

int main(int nArgs, char **ArgV){

  // Initialse and parse input arguments:
  RunControl* RC = RunControl::getInstance();
  RC->ParseArgs(nArgs, ArgV);
  if ( RC->getDebug() )
    RC->print();

  UserAnal* A = new UserAnal(RC->getDebug());
  A->PreEventLoop();
  A->EventLoop();
  A->PostEventLoop();
  A->HistFitDo();

  }
