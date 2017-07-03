#ifndef CutflowCheck_h
#define CutflowCheck_h

#include "AnalyzerCore.h"
class CutflowCheck : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  CutflowCheck();
  ~CutflowCheck();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( CutflowCheck, 1);
};
#endif
