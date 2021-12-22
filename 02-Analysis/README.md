# nuSTORM -- nuAnalis

The code in this directory tree provides a c++ analysis skeleton that can be used to read the nuSIM ntuples (see [nuSTORM](https://www.nustorm.org/trac/).  Documentation for the nuAnalysis package will be posted on the nuSTORM wiki.

## To set up and run:
Execute "startup.bash" from this directory (i.e. run the bash command "source startup.bash").  This will:
  * Set the "nuAnalysisPATH" to this director; and
  * Set "nuAnalysisCPATH" so that g++ can see the 01-Code directory.  The scripts in "02-Tests" may then be run using the compile scriptis included in 02-Tests.  A sample user analysis programme is included in "03-Skeleton"

## Directories:
 * C++ classes and "library" code stored in "01-Code"
 * Test scripts stored in "02-Tests"
 * Integration tests are stored in "03-Integration-Test"
 * Sample data ntuples are stored in "31-Data"

## Dependencies:
 * g++11 and assumes ROOT is installed.

## History
 * 06 December 2021:  First version.
