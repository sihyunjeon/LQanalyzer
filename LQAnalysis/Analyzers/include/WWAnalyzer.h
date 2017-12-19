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
  std::vector<KLepton> SortByPtOrder( std::vector<KLepton> leptons );
  double GetProjectedMET( snu::KParticle MET, std::vector<KLepton> leptons );

  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( WWAnalyzer, 1);
};
#endif
