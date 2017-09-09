#ifndef CFRateCalculator2_h
#define CFRateCalculator2_h

#include "AnalyzerCore.h"
class CFRateCalculator2 : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  CFRateCalculator2();
  ~CFRateCalculator2();

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


  ClassDef ( CFRateCalculator2, 1);
};
#endif
