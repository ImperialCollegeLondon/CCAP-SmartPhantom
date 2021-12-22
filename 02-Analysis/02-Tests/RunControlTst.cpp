/* Tests of RunControl class */

#include <iostream>
#include <string.h>

#include "RunControl.hpp"

int main(int nArgs, char *ArgV[]){

  bool Debug = false;
  char *Arg;
  for (int i = 0 ; i < nArgs ; i++) {
    Arg = ArgV[i];
    if ( strcmp(Arg, "-d") == 0 )
      Debug = true;
  }
  
  if ( Debug ) std::cout  << "Start tests of RunControl class:" << std::endl;

  int i = 0;
  // Test of instantiation of singleton class:
  if ( Debug ) std::cout  << "    Test " << i << ": test singleton nature of class:"
			  << std::endl;

  // Initialise run control singleton class:
  RunControl* RC = RunControl::getInstance();
  if ( Debug ) std::cout  << "    ----> Initial RunControl instance: "
			  << RC << std::endl;
  if ( Debug ) RC->print();
  RunControl* RC1 = RunControl::getInstance();
  if ( Debug ) std::cout  << "    ----> Second RunControl instance: " << RC1 << std::endl;
  if (RC != RC1) {
    std::cout  << "    ----> FAILED!  Not a singleton!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  else {
    if ( Debug ) std::cout  << "        ----> SUCCESS!  A singleton!" << std::endl;
  }
  i++;
  // Test parsing of input arguments:
  if ( Debug ) std::cout  << "    Test " << i << ": test parsing of input arguments:"
			  << std::endl;
  if ( Debug ) std::cout  << "    ----> Number of input arguments: "
			  << nArgs << std::endl;
  std::string Hlp = "-h";
  const char* cHlp = Hlp.c_str();
  char *ArgVPlus[2];
  int nArgsPlus;
  if ( nArgs == 1 ){
    if ( Debug ) std::cout  << "    ----> No arguments given, set -h, help flag."
	      << std::endl;
    ArgVPlus[0] = ArgV[0]; 
    ArgVPlus[1] = const_cast<char *>(cHlp);
    if ( Debug ) std::cout  << "; ArgPlus: " << ArgVPlus[nArgs] << std::endl;
    nArgsPlus = nArgs + 1;
    if ( Debug ) std::cout  << "        ----> Updated: nArgs, ArgV: "
			    << nArgsPlus << ": " << ArgVPlus << std::endl;
    RC->ParseArgs(nArgsPlus, ArgVPlus);
  }
  else {
    RC->ParseArgs(nArgs, ArgV);
  }
  if ( Debug ) std::cout  << "    ----> RunControl parameters:" << std::endl;
  if ( Debug ) RC->print();
  //--> Test getters:
  i++;
  if ( Debug ) std::cout  << "    Test " << i << ": test getters:"
	    << std::endl;
  bool        Dbg = RC->getDebug();
  bool        FlF = RC->getFileFlag();
  bool        ChF = RC->getChainFlag();
  bool        OPF = RC->getOutPutFlag();
  std::string FlN = RC->getROOTfilename();
  std::string ChN = RC->getCHAINdirname();
  std::string OPN = RC->getOUTPUTdirname();
  if ( Debug ) std::cout  << "RunControl initialised with parameters:" << std::endl;
  if ( Debug ) std::cout  << "                Debug: " << Dbg << std::endl;
  if ( Debug ) std::cout  << "             FileFlag: " << FlF << std::endl;
  if ( Debug ) std::cout  << "            ChainFlag: " << ChF << std::endl;
  if ( Debug ) std::cout  << "           OutPutFlag: " << OPF << std::endl;
  if ( Debug ) std::cout  << "       ROOT file name: " << FlN << std::endl;
  if ( Debug ) std::cout  << " ROOT chain directory: " << ChN << std::endl;
  if ( Debug ) std::cout  << "     OutPut directory: " << OPN << std::endl;
  
}
