#ifndef CFRateCalculator_FR_h
#define CFRateCalculator_FR_h

#include "AnalyzerCore.h"
class CFRateCalculator_FR : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  CFRateCalculator_FR();
  ~CFRateCalculator_FR();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );

  void DrawHistograms( snu::KElectron lep, int eta_region, double weight, bool is_CF, bool is_Z );
  void GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index );
  float GetPrescale( std::vector<snu::KElectron> electronColl, bool pass_low, bool pass_high );
  snu::KParticle ShiftEnergy( snu::KParticle old_lep, double shift_rate );
  void CFvalidation(void);
  double Get2DCFRates(bool apply_sf, double el_pt, double el_eta);

  TH2F* TIGHT1_CF_hist;
  TH2F* TIGHT2_CF_hist;

  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( CFRateCalculator_FR, 1);
};
#endif
