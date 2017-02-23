#ifndef HNSSSFMuMuE_CR_h
#define HNSSSFMuMuE_CR_h

#include "AnalyzerCore.h"
class HNSSSFMuMuE_CR : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNSSSFMuMuE_CR();
  ~HNSSSFMuMuE_CR();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();

  snu::KParticle SF[2], OF;
  void DrawHistograms(TString suffix, snu::KParticle SF[], snu::KParticle OF, snu::KParticle MET,  std::vector<snu::KJet> jetTightColl, double weight);
 

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HNSSSFMuMuE_CR, 1);
};
#endif
