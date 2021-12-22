// Read the previously produced N-Tuple and print on screen
// its content

#include <iostream>

#include "RunControl.hpp"
#include "Analysis.hpp"

int main(int nArgs, char **ArgV){

  // Initialse and parse input arguments:
  RunControl* RC = RunControl::getInstance();
  RC->ParseArgs(nArgs, ArgV);
  if ( RC->getDebug() )
    RC->print();

  // Initialse skeleton analysis:
  Analysis* A = new Analysis(RC->getDebug());
  A->PreEventLoop();
  A->EventLoop();
  A->PostEventLoop();
  A->HistFitDo();
  
  // Execute built-in Bragg peak analysis test:
  BraggPeak* B = new BraggPeak(RC->getDebug());
  B->PreEventLoop();
  B->EventLoop();
  B->PostEventLoop();
  B->HistFitDo();
  
}
