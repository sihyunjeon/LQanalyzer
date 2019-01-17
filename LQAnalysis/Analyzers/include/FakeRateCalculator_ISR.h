#ifndef FakeRateCalculator_ISR_h
#define FakeRateCalculator_ISR_h

#include "AnalyzerCore.h"
class FakeRateCalculator_ISR : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  FakeRateCalculator_ISR();
  ~FakeRateCalculator_ISR();

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

  void GetFakeRatesAndCorrectionFactors(std::vector<snu::KMuon> muons, double this_weight);
  void GetFakeRatesAndCorrectionFactors(std::vector<snu::KElectron> electrons, double this_weight);

  void GetFakeRates(TString flavour, KLepton Probe, std::vector<snu::KJet> jets, bool Probe_passing_Tight, double this_weight);
  void GetFakeRates_FillHistograms(KLepton Probe, snu::KJet Tag, snu::KParticle MET, TString cut, TString workingpoint, TString flavour, double this_weight);

  void GetCorrectionFactors(TString flavour, std::vector<KLepton> leptons, std::vector<snu::KJet> jets, double this_weight, bool both_pass_Tight);
  void GetCorrectionFactors_FillHistograms(TString flavour, std::vector<KLepton> leptons, std::vector<snu::KJet> jets, TString workingpoint, TString trigger, double this_weight);

  void DoMCClosureTests(std::vector<snu::KMuon> muons, std::vector<bool> is_Prompt, double this_weight);
  void DoMCClosureTests(std::vector<snu::KElectron> electrons, std::vector<bool> is_Prompt, double this_weight);
  void DoMCClosureTests(TString flavour, std::vector<KLepton> leptons, std::vector<bool> is_Tight, std::vector<bool> is_Prompt, double this_weight);

  bool SelectZpeak(std::vector<KLepton> TightLep, double zmass, double window);
  double GetTransverseMass(KLepton, snu::KParticle);

  double GetWeightFromFakeRate(TString flavour, std::vector<KLepton> leptons, std::vector<bool> is_Tight);
  double GetFakeRateForWeight(TString flavour, KLepton Lepton);


  ClassDef ( FakeRateCalculator_ISR, 1);
};
#endif
