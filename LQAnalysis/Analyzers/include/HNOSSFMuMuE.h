#ifndef HNOSSFMuMuE_h
#define HNOSSFMuMuE_h

#include "AnalyzerCore.h"

class HNOSSFMuMuE : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNOSSFMuMuE();
  ~HNOSSFMuMuE();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();


  snu::KParticle GENmu[2], GENel, GENnu, GENHN;
  snu::KParticle RAWmu[2], RAWel, RAWnu[2];
  snu::KParticle RECOmu[2], RECOel, RECOnu_lowmass, RECOnu_highmass, RECOW_pri_lowmass, RECOW_sec_lowmass, RECOW_pri_highmass, RECOW_sec_highmass, RECOHN[4];
  snu::KParticle MET;
  void GENSignalStudy( bool doGENEventSelection );
  void GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index );
  void GENEventSelectionStudy( snu::KParticle GENmu[], snu::KParticle GENel, snu::KParticle GENnu, snu::KParticle GENHN );
  void EventSelectionStudy( snu::KParticle RAWmu[], snu::KParticle RAWel, int signal_class );
  void DrawHistograms( TString suffix, double weight );
  int DefineClass();
  int n_bjets;
   int GetPeriodIndex(void);
  TString PutString_PassOptimizedCuts(int sig_mass);
  bool PassOptimizedCuts(double first_pt, double second_pt, double third_pt, double METPt, double RECOW_pri_highmass, double RECOW_pri_lowmass, double RECOW_sec_highmass, double RECOW_sec_lowmass, int sig_mass);
    

 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HNOSSFMuMuE, 1);
};
#endif
