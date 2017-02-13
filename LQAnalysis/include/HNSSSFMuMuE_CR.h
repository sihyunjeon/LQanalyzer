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

  snu::KParticle RAWmu[3], RAWel, RAWnu[2];
  snu::KParticle RECOnu;
  void EventSelectionStudy( snu::KParticle RAWmu[], snu::KParticle RAWel, std::vector<snu::KJet> jetColl, snu::KParticle RECOnu, int control_region, double weight );
  void StudyControlRegion( TString suffix, snu::KParticle RAWmu[], snu::KParticle RAWel, std::vector<snu::KJet> jetColl, snu::KParticle RECOnu, int control_region, double weight );
    

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
