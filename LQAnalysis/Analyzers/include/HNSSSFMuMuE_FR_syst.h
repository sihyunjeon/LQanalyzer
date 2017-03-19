#ifndef HNSSSFMuMuE_FR_syst_h
#define HNSSSFMuMuE_FR_syst_h

#include "AnalyzerCore.h"

class HNSSSFMuMuE_FR_syst : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNSSSFMuMuE_FR_syst();
  ~HNSSSFMuMuE_FR_syst();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  TString GetSuffix(std::vector<double> RelIso_array, int RelIso_bin, std::vector<double> dXYSig_array, int dXYSig_bin);
  TString GetSuffix(std::vector<double> awayJetPt_array, int awayJetPt_bin);

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HNSSSFMuMuE_FR_syst, 1);
};
#endif
