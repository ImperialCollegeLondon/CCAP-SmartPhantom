/*    Analysis class header file */

#include "TChain.h"
#include "TH1.h"

class Analysis {

private:
  
  bool        Debug;  // Global debug flag

  TChain* Event_ch;  // Initalise attribute of type TChain for Event ntuple

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
  void              setEvent_ch( TChain* Evt_ch );
  
  //--> Getters:
  bool                 getDebug(){ return Debug; };
  TChain*           getEvent_ch(){ return Event_ch; };

};

class BraggPeak : public Analysis{

public: 
  BraggPeak( bool Dbg=true );
  ~BraggPeak() {  }

  void PreEventLoop( bool Dbg=true ); 
  void EventLoop( bool Dbg=true );
  void PostEventLoop( bool Dbg=true );

};

class UserAnal : public Analysis{

public: 
  UserAnal( bool Dbg=true );
  ~UserAnal() {  }

  void PreEventLoop( bool Dbg=true ); 
  void EventLoop( bool Dbg=true );
  void PostEventLoop( bool Dbg=true );

};
