#include <cmath>
#include <fstream>
#include <vector>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TRandom.h"
#include <numeric>

void energy(TString filename, TString rootTree, TString name, TString outfile){
  TFile *myfile = TFile::Open(filename);
  
  TTreeReader myReader(rootTree, myfile);
  
  TTreeReaderValue<double> myE(myReader, name+".Edep");
  TTreeReaderValue<double> myTime(myReader, name+".HitTime");
  TTreeReaderValue<double> myEventN(myReader, name+".EventN");
  TTreeReaderValue<double> myX(myReader, name+".PosX");
  TTreeReaderValue<double> myY(myReader, name+".PosY");
  TTreeReaderValue<double> myZ(myReader, name+".PosZ");
  TTreeReaderValue<double> myDeltaT(myReader, name+".DeltaT");
  //TTreeReaderValue<std::string> myName(myReader, "Name.c_str()");
  TTreeReaderValue<std::string> myName(myReader,"Name");
  
  double count = 0;
  
  std::vector<double> eventNum;
  std::vector<double> elossVec;
  std::vector<double> fibNVec;
  
  std::cout << "Reading '" << filename << "'." << std::endl;

  // If first station then generate new file and write header, else append to existing file
  ofstream file;
  file.open(outfile);
  file << "StationName" << '\t' << "EventNum" << '\t' << "Total Energy Deposit [MeV]" << '\t' << "True PosX [mm]" << '\t' << "True PosY [mm]" << '\t' << "True PosZ [mm]" << '\t' << "Global Time [ns]" << '\t' << "Delta Time [ns]" << '\t' << "Particle Name" << endl;      
 
  while(myReader.Next()){
    double edep = *myE;
    double time = *myTime;
    double eventN = *myEventN;
    
    int letters = myName.Get()->size(); // Number of letters in particle name
    std::string pname = "";
    for(int i=0; i<letters; i++)
    {
      pname += myName->at(i);
    }

    file << name << '\t' << eventN << '\t' << edep << '\t' << *myX << '\t' << *myY << '\t' << *myZ << '\t' << time << '\t' << *myDeltaT << '\t' << pname <<endl;
  }
  
  file.close();
  
  std::cout << "Output to '" << outfile << "'." << std::endl;
}


void root2Ascii(){
  TString rootfilename = "waterhits.root";
  TString outfilename = "waterhits.dat";
  energy(rootfilename, "Event", "WaterBox", outfilename);
}
