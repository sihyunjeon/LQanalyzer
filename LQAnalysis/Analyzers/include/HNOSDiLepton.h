#ifndef HNOSDiLepton_h
#define HNOSDiLepton_h

#include "AnalyzerCore.h"
class HNOSDiLepton : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNOSDiLepton();
  ~HNOSDiLepton();

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


  ClassDef ( HNOSDiLepton, 1);
};
#endif
