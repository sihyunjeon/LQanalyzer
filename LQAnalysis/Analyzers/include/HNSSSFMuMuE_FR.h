#ifndef HNSSSFMuMuE_FR_h
#define HNSSSFMuMuE_FR_h

#include "AnalyzerCore.h"

class HNSSSFMuMuE_FR : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNSSSFMuMuE_FR();
  ~HNSSSFMuMuE_FR();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();


  snu::KParticle RAWmu[2], RAWel, RAWnu[2];
  snu::KParticle RECOmu[2], RECOel, RECOnu_lowmass, RECOnu_highmass, RECOW_pri_lowmass, RECOW_sec_lowmass, RECOW_pri_highmass, RECOW_sec_highmass, RECOHN[4];
  void EventSelectionStudy( snu::KParticle RAWmu[], snu::KParticle RAWel, int signal_class );
  void DrawHistograms( TString suffix, double weight );
  int DefineClass();

  double weight_err;    

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HNSSSFMuMuE_FR, 1);
};
#endif
