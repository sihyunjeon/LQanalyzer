#ifndef DiEl_SS_h
#define DiEl_SS_h

#include "AnalyzerCore.h"
class DiEl_SS : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  DiEl_SS();
  ~DiEl_SS();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();

  double Get2DCFRates(bool apply_sf, double el_pt, double el_eta, TString el_ID, TString halfsample, TString CFsample);
  snu::KParticle ShiftEnergy( snu::KParticle old_lep, double shift_rate );

  TH2F* HNTIGHT_CF_hist_madgraph;
  TH2F* HNTIGHT_CF_sampleA_hist_madgraph;
  TH2F* HNTIGHT_CF_sampleB_hist_madgraph;
  TH2F* MVATIGHT_CF_hist_madgraph;
  TH2F* MVATIGHT_CF_sampleA_hist_madgraph;
  TH2F* MVATIGHT_CF_sampleB_hist_madgraph;

  TH2F* HNTIGHT_CF_hist_powheg;
  TH2F* HNTIGHT_CF_sampleA_hist_powheg;
  TH2F* HNTIGHT_CF_sampleB_hist_powheg;
  TH2F* MVATIGHT_CF_hist_powheg;
  TH2F* MVATIGHT_CF_sampleA_hist_powheg;
  TH2F* MVATIGHT_CF_sampleB_hist_powheg;



 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( DiEl_SS, 1);
};
#endif
