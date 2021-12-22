// Read the previously produced N-Tuple and print on screen
// its content

#include <iostream>

#include "RunControl.hpp"
#include "nuAnalysis.hpp"

int main(int nArgs, char **ArgV){

  // Initialse and parse input arguments:
  RunControl* RC = RunControl::getInstance();
  RC->ParseArgs(nArgs, ArgV);
  if ( RC->getDebug() )
    RC->print();

  // Initialse analysis:
  nuAnalysis* nuA = new nuAnalysis(RC->getDebug());
  nuA->PreEventLoop();
  nuA->EventLoop();
  nuA->PostEventLoop();
  nuA->HistFitDo();
  
  nuSIMtstRestFrame* nuB = new nuSIMtstRestFrame(RC->getDebug());
  nuB->PreEventLoop();
  nuB->EventLoop();
  nuB->PostEventLoop();
  nuB->HistFitDo();
  
}
