/*  Run control class implementation file */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string.h>
#include "sys/stat.h"

#include "RunControl.hpp"

RunControl* RunControl::instance = NULL;

RunControl* RunControl::getInstance() {
    if (instance == NULL) {
      
        instance     = new RunControl();
    }
    return instance;
}
  
RunControl::RunControl(bool Dbg,
		       std::string Rfn,
		       std::string Cdn,
		       std::string Odn){
  Debug         = Dbg;
  ROOTfilename  = Rfn;
  CHAINdirname  = Cdn;
  OUTPUTdirname = Odn;
}

void RunControl::print() {
  std::cout << " RunControl::print: RunControl initialised with parameters:"
	    << std::endl;
  std::cout << "                 Debug: " << Debug         << std::endl;
  std::cout << "        ROOT file name: " << ROOTfilename  << std::endl;
  std::cout << "  ROOT chain directory: " << CHAINdirname  << std::endl;
  std::cout << " OUTPUT directory name: " << OUTPUTdirname << std::endl;
}

void DumpHelp() {
  std::cout << " RunControl help: explanation of arguments:"
	    << std::endl;
  std::cout << "\t -h \t\t Generates this help printout to show use of flags"
            << " etc." << std::endl;
  std::cout << "\t -d \t\t Sets debug flag: RunControl::Debug = true"
	    << std::endl;
  std::cout << "\t -f <filename> \t Sets ROOT <filename> to be read "
	    << "(single file)." << std::endl;
  std::cout << "\t\t\t If -f is specified, <filename> must exist or "
	    << "execution is" << std::endl;
  std::cout << "\t\t\t terminated" << std::endl;
  std::cout << "\t -c <dir name> \t Sets directory containing ROOT "
	    << "files to be chained."
	    << std::endl;
  std::cout << "\t\t\t If -c is specified, <dir name> must exist or "
	    << "execution is" << std::endl;
  std::cout << "\t\t\t terminated" << std::endl;
  std::cout << "\t -o <dir name> \t Sets OUTPUT directory.  Defaul=./ "
	    << std::endl;
  std::cout << "\t\t\t If -o is specified, <dir name> must exist or "
	    << "execution is" << std::endl;
  std::cout << "\t\t\t terminated" << std::endl;
  std::cout << "\t Note that -f and -c can be used together, all files "
	    << "will be read." << std::endl;
}

void RunControl::ParseArgs(int nArgs, char *ArgV[]){
  char *Arg;

  //--> Defaults:
  Debug      = false;
  FileFlag   = false;
  ChainFlag  = false;
  OutPutFlag = false;

  //--> Scan input arguments:
  for (int i = 0 ; i < nArgs ; i++) {
    Arg = ArgV[i];
    if ( strcmp(Arg, "-d") == 0 )
      Debug = true;
    else if ( strcmp(Arg, "-f") == 0 ) {
      FileFlag     = true;
      if ( i+1 < nArgs)
	ROOTfilename = ArgV[i+1];
      else {
	std::cout << " RunControl::ParseArgs:Error! Void file name; STOP!"
		  << std::endl;
	std::exit(EXIT_FAILURE);
      }
    }
    else if ( strcmp(Arg, "-c") == 0 ) {
      ChainFlag    = true;
      if ( i+1 < nArgs)
	CHAINdirname = ArgV[i+1];
      else {
	std::cout <<
	  " RunControl::ParseArgs:Error! Void chain directory name; STOP!"
		  << std::endl;
	std::exit(EXIT_FAILURE);
      }
    }
    else if ( strcmp(Arg, "-o") == 0 ) {
      OutPutFlag    = true;
      if ( i+1 < nArgs)
	OUTPUTdirname = ArgV[i+1];
      else {
	std::cout <<
	  " RunControl::ParseArgs:Error! Void output directory name; STOP!"
		  << std::endl;
	std::exit(EXIT_FAILURE);
      }
    }
    else if ( strcmp(Arg, "-h") == 0 ) 
      DumpHelp();
  }

  //--> Print em if Debug
  if ( Debug == true) {
    std::cout << " RunControl::ParseArgs:Debug: parsed "
	      << nArgs << " arguments \n" ;
    for (int i = 0 ; i < nArgs ; i++) {
      std::cout << "     ----> Argument: " << i 
		<< " value: " << ArgV[i] << "\n" ;
    }
    std::cout << "     ----> Flags: "
	      << "Debug=" << Debug << ", "
	      << "FileFlag=" << FileFlag << ", "
	      << "ChainFlag=" << ChainFlag << ", "
	      << "OutPutFlag=" << OutPutFlag << ", "
	      << std::endl;
  }

  //--> Validate file and directories:
  std::ifstream f(ROOTfilename);
  bool File = f.good();
  struct stat db;
  bool Dir  = (stat(CHAINdirname.c_str(), &db) == 0);
  bool Out  = (stat(OUTPUTdirname.c_str(), &db) == 0);
  if ( FileFlag and !File ) {
    std::cout << " RunControl::ParseArgs:Error! File "
	      << ROOTfilename << " does not exist, STOP!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ( ChainFlag and !Dir ) {
    std::cout << " RunControl::ParseArgs:Error! Directory "
	      << CHAINdirname << " does not exist. STOP!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ( OutPutFlag and !Out ) {
    std::cout << " RunControl::ParseArgs:Error! Directory "
	      << OUTPUTdirname << " does not exist. STOP!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}
