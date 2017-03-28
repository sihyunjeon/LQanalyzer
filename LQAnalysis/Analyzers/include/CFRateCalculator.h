#ifndef CFRateCalculator_h
#define CFRateCalculator_h

#include "AnalyzerCore.h"
class CFRateCalculator : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  CFRateCalculator();
  ~CFRateCalculator();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );

  void DrawHistograms( snu::KElectron lep, int eta_region, double weight, bool is_CF, bool is_Z );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( CFRateCalculator, 1);
};
#endif
