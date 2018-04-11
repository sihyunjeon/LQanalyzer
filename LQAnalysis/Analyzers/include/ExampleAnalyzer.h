#ifndef ExampleAnalyzer_h
#define ExampleAnalyzer_h

#include "AnalyzerCore.h"
class ExampleAnalyzer : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  ExampleAnalyzer();
  ~ExampleAnalyzer();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );

bool PassEMuTriggerPt(std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons);
void GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index );
  std::vector<TString> triggerlist_emBG1;
  std::vector<TString> triggerlist_emBG2;
  std::vector<TString> triggerlist_emH1;
  std::vector<TString> triggerlist_emH2;  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( ExampleAnalyzer, 1);
};
#endif
