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
  void CFvalidation(void);
  snu::KParticle ShiftEnergy( snu::KParticle old_lep, double shift_rate );


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
