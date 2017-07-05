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
  void GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index );
  float GetPrescale( std::vector<snu::KElectron> electronColl, bool pass_low, bool pass_high );
  snu::KElectron ShiftEnergy( snu::KElectron old_lep, double shift_rate );
  void CFvalidation(void);
  double GetCFRates(int sys, double el_pt, double el_eta, TString el_ID);
  double GetCFweight(int sys, std::vector<snu::KElectron> electron, bool apply_sf, TString Zwidth);
  void CF_MCClosure(std::vector<snu::KElectron> electron);

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
