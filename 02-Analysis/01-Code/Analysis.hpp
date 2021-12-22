/*    Analysis class header file */

#include "TChain.h"
#include "TH1.h"

class Analysis {

private:
  
  bool        Debug;  // Global debug flag

  TChain* Event1_ch;  // Initalise attribute of type TChain for Event1 ntuple
  TChain*  Event_ch;  // Initalise attribute of type TChain for Event ntuple

public:
  std::vector<TH1F*> TH1Flist;
  std::vector<TF1*>  TF1list;

  Analysis( bool Dbg=true );

  ~Analysis() {  }

  void PreEventLoop( bool Dbg=true );
  
  void EventLoop( bool Dbg=true );

  void PostEventLoop( bool Dbg=true );

  void HistFitDo( );

  //--> Setters:
  void                 setDebug( bool Dbg );
  void             setEvent1_ch( TChain* rInfo_ch );
  void              setEvent_ch( TChain* bm_ch );
  
  //--> Getters:
  bool                 getDebug(){ return Debug; };
  TChain*          getEvent1_ch(){ return Event1_ch; };
  TChain*           getEvent_ch(){ return Event_ch; };

};
