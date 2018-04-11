#ifndef CFRateCalculator_Final_h
#define CFRateCalculator_Final_h

#include "AnalyzerCore.h"
class CFRateCalculator_Final : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  CFRateCalculator_Final();
  ~CFRateCalculator_Final();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );

  void GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index );
  float GetCFweight(std::vector<snu::KElectron> electrons, bool apply_sf, TString el_ID, bool do_halftest);
  float GetCFRates(double el_pt, double el_eta, TString el_ID, bool do_halftest);
  snu::KElectron ShiftEnergy( snu::KElectron old_lep, double shift_rate );
  bool JHsConv( snu::KElectron lepton );

  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( CFRateCalculator_Final, 1);
};
#endif
