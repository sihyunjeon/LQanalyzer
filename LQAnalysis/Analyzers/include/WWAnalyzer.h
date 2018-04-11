#ifndef WWAnalyzer_h
#define WWAnalyzer_h

#include "AnalyzerCore.h"
class WWAnalyzer : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  WWAnalyzer();
  ~WWAnalyzer();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  void DoTruthMCStudy( void );

  void DrawHistColl( TString this_string,
                     std::vector<KLepton> leptons,
                     std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets,
                     snu::KParticle MET,
                     double this_weight );
  void DrawMAOSDifference(TString this_string, snu::KParticle MT2assistedMETreco_p, snu::KParticle MT2assistedMETreco_m, double this_weight);

  void DrawAsymmetryHist( TString this_string,
                          snu::KParticle wp, snu::KParticle wm,
                          double this_weight );

  std::vector<KLepton> SortByPtOrder( std::vector<KLepton> leptons );
  double GetProjectedMET( snu::KParticle MET, std::vector<KLepton> leptons );

  void GetDijetMassClosest( std::vector<snu::KJet> jets, double target_mass, int& m, int& n);
  void GetNuPzFromWMass( KLepton lepton, snu::KParticle MET, std::vector<snu::KParticle>& reconstructedMET, bool& img );
  void GENFindDecayIndex( std::vector<snu::KTruth> truthColl, int it, std::vector<int>& index );
  double GetMT( KLepton this_lepton, double this_MET );
  double GetMT( snu::KParticle this_lepton, double this_MET );
  double GetMT2( std::vector<KLepton> leptons, snu::KParticle MET, std::vector<snu::KParticle>& MT2_assisted_nu );
  double GetMT2( std::vector<snu::KParticle> leptons, snu::KParticle MET, std::vector<snu::KParticle>& MT2_assisted_nu );

  //
  // The output variables 
  //
  /// Vectors for output objetcs
  snu::KParticle gen_nup, gen_num;
  
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( WWAnalyzer, 1);
};
#endif
