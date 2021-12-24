Documentation for Analysis package
==================================

Cheep and cheeful documentation for the Analysis package.

RunControl:
-----------
The RunControl class is conceived as providing the funtionality
necessary to control the user Analysis.  In particular it provides the
file-handling code to open and write ROOT files.  It will also be
developed and maintained to handle machine dependences at the i/o
level.  The idea is to make it easier for people to get going on the
Analysis.

The action of the RunControl class is to pass options to the Analysis
package.  The following options are implemented:

              -h : Generates "help" printout to define flags and
                   options etc.
		   
              -d : Sets the debug flag: RunControl::Debug = true
	      
   -f <filename> : Sets ROOT <filename> to be read (single file).  If -f
                   is specified, <filename> must exist or execution is
                   terminated
		   
   -c <dir name> : Sets directory containing ROOT files to be chained. If -c
                   is specified, <dir name> must exist or execution is
                   terminated
		   
   -o <dir name> : Sets output directory for files to be put in. If -o
                   is specified, <dir name> must exist or execution is
                   terminated.  If -o is not specified, files will be
                   written to the present working directory.
		   
Note that -f and -c can be used together, all files will be read.

