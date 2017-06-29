#ifndef HNEMuAnalyzer_h
#define HNEMuAnalyzer_h

#include "AnalyzerCore.h"


class HNEMuAnalyzer : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNEMuAnalyzer();
  ~HNEMuAnalyzer();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
  void FillEventCutFlow(TString cut, TString label, float w);
  void FillIsoCutFlow(TString cut, float w);
  void CheckJetsCloseToLeptons(std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets,  TString name);
  

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;

  ClassDef ( HNEMuAnalyzer, 1);
};
#endif
