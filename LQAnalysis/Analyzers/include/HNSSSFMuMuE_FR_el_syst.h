#ifndef HNSSSFMuMuE_FR_el_syst_h
#define HNSSSFMuMuE_FR_el_syst_h

#include "AnalyzerCore.h"

class HNSSSFMuMuE_FR_el_syst : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNSSSFMuMuE_FR_el_syst();
  ~HNSSSFMuMuE_FR_el_syst();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  int GetPeriodIndex(void);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HNSSSFMuMuE_FR_el_syst, 1);
};
#endif
