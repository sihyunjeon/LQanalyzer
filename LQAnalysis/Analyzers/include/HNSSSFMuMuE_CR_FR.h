#ifndef HNSSSFMuMuE_CR_FR_h
#define HNSSSFMuMuE_CR_FR_h

#include "AnalyzerCore.h"
class HNSSSFMuMuE_CR_FR : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNSSSFMuMuE_CR_FR();
  ~HNSSSFMuMuE_CR_FR();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();

  bool CR_WZ_mumue, CR_Zjet_mumue, CR_ttW_mumue, CR_ttZ_mumue;
  snu::KParticle SF[2], OF;
  void DrawHistograms(TString suffix, snu::KParticle SF[], snu::KParticle OF, snu::KParticle MET,  std::vector<snu::KJet> jetTightColl, double weight);

  TH2F* hist_single_highdxy_trilep;
  TH2F* hist_single_dijet;
  double GetMuonFakeRate( TString method, snu::KMuon muon );
  int GetPeriodIndex(void);

  double weight_err;

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HNSSSFMuMuE_CR_FR, 1);
};
#endif
