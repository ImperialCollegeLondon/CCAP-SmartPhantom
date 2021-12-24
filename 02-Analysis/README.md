# nuSTORM -- nuAnalis

The code in this directory tree provides a c++ analysis skeleton that can be used to read the SmartPhantom ntuples (see [CCAP/LhARA/IonAcoustic](https://ccap.hep.ph.ic.ac.uk/trac/wiki/Research/Instrumentation/IonAcoustic).  Documentation for the Analysis package will be posted on the CCAP/LhARA/IonAcoustic wiki.

## To set up and run:
Execute "startup.bash" from this directory (i.e. run the bash command "source startup.bash").  This will:
  * Set the "AnalysisPATH" to this director; and
  * Set "AnalysisCPATH" so that g++ can see the 01-Code directory.  The scripts in "02-Tests" may then be run using the compile scripts included in 02-Tests.  A sample user analysis programme is included in "03-Skeleton"
  * The actions of the Analysis code are controlled using the "RunControl" class.  The options defined by this class and how to control the input and output ROOT files may be found in the RunControl documentation.

## Directories:
 * C++ classes and "library" code stored in "01-Code"
 * Test scripts stored in "02-Tests"
 * Integration tests are stored in "03-Integration-Test"
 * Sample data ntuples are stored in "31-Data"
 * Scratch space for output files etc can be found in 99-Scratch.
   Files in this directory are not included in the repository.

## Dependencies:
 * g++11 and assumes ROOT is installed.

## History
 * 22 December 2021:  First version.
