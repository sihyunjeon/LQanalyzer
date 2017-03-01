#ifndef FakeRateCalculator_Mu_h
#define FakeRateCalculator_Mu_h

#include "AnalyzerCore.h"
class FakeRateCalculator_Mu : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  FakeRateCalculator_Mu();
  ~FakeRateCalculator_Mu();

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
  double GetPrescale(std::vector<snu::KMuon> muonColl, bool passlow, bool passhigh);
  void StudyDijetMuon(TString suffix, snu::KMuon muon, double weight);
  void StudyHighdXYMuon(TString suffix, snu::KMuon muon, double weight);


  ClassDef ( FakeRateCalculator_Mu, 1);
};
#endif
