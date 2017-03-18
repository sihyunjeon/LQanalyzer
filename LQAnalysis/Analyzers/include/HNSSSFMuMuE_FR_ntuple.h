#ifndef HNSSSFMuMuE_FR_ntuple_h
#define HNSSSFMuMuE_FR_ntuple_h

#include "AnalyzerCore.h"

class HNSSSFMuMuE_FR_ntuple : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNSSSFMuMuE_FR_ntuple();
  ~HNSSSFMuMuE_FR_ntuple();

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
  snu::KParticle MET;
  void EventSelectionStudy( snu::KParticle RAWmu[], snu::KParticle RAWel, int signal_class );
  void DrawHistograms( TString suffix, double weight );
  TString GetSystematicString( int it_sys );

  int DefineClass();
  int n_bjets;
   int GetPeriodIndex(void);

    

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HNSSSFMuMuE_FR_ntuple, 1);
};
#endif
