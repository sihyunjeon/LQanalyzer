#ifndef HNSSSFMuMuE_h
#define HNSSSFMuMuE_h

#include "AnalyzerCore.h"

class HNSSSFMuMuE : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNSSSFMuMuE();
  ~HNSSSFMuMuE();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();


  snu::KParticle GENmu[2], GENel, GENnu, GENHN;
  snu::KParticle RAWmu[2], RAWel, RAWnu[2];
  snu::KParticle RECOmu[2], RECOel, RECOnu_lowmass, RECOnu_highmass, RECOW_pri_lowmass, RECOW_sec_lowmass, RECOW_pri_highmass, RECOW_sec_highmass, RECOHN[4];
  snu::KParticle MET;
  void GENSignalStudy( bool doGENEventSelection );
  void GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index );
  void GENEventSelectionStudy( snu::KParticle GENmu[], snu::KParticle GENel, snu::KParticle GENnu, snu::KParticle GENHN );
  void EventSelectionStudy( snu::KParticle RAWmu[], snu::KParticle RAWel, int signal_class );
  void DrawHistograms( TString suffix, double weight );
  int DefineClass();
    

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HNSSSFMuMuE, 1);
};
#endif
