/*      -------- nuAnalysis:  Run control class header file --------      */
// RunControl call is a singleton.  Its role is to set flags and to allow
// various file names to be defined from the command line.
// The various flags can be printed with teh option "-h".  The flags are
// defined as follows:
//              -h : Generates "help" printout to show use of flags etc.
//              -d : Sets debug flag: RunControl::Debug = true
//   -f <filename> : Sets ROOT <filename> to be read (single file).  If -f
//                   is specified, <filename> must exist or execution is
//                   terminated
//   -c <dir name> : Sets directory containing ROOT files to be chained. If -c
//                   is specified, <dir name> must exist or execution is
//                   terminated
//   -o <dir name> : Sets output directory for files to be put in. If -o
//                   is specified, <dir name> must exist or execution is
//                   terminated
// Note that -f and -c can be used together, all files will be read.
class RunControl {

private:
  static RunControl* instance;

  bool                Debug;    // Global debug flag
  bool             FileFlag;    // Global flag to say read file
  bool            ChainFlag;    // Global flag to say chain files in dir
  bool           OutPutFlag;    // Global flag to say output file specified
  std::string  ROOTfilename;    // ROOT input file name
  std::string  CHAINdirname;    // Director for chain of ROOT files
  std::string OUTPUTdirname;    // Output directory path
  
  RunControl(bool Dbg = true,
	     std::string Rfn = "Initaliation dummy",
	     std::string Cdn = "Initaliation dummy",
     	     std::string Odn = "./");

  ~RunControl() { };
  
public:
  //--> Getters:
  bool                   getDebug(){ return Debug; };
  bool                getFileFlag(){ return FileFlag; };
  bool               getChainFlag(){ return ChainFlag; };
  bool              getOutPutFlag(){ return OutPutFlag; };
  std::string     getROOTfilename(){ return ROOTfilename; }
  std::string     getCHAINdirname(){ return CHAINdirname; }
  std::string    getOUTPUTdirname(){ return OUTPUTdirname; }
  static RunControl*  getInstance();

  //--> Dumpers:
  void print();

  //--> Dumpers:
  void ParseArgs(int nArgs, char *ArgV[]);

};
