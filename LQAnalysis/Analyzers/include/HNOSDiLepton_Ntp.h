#ifndef HNOSDiLepton_Ntp_h
#define HNOSDiLepton_Ntp_h

#include "AnalyzerCore.h"
class HNOSDiLepton_Ntp : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNOSDiLepton_Ntp();
  ~HNOSDiLepton_Ntp();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  bool Pass_Preselection;
  bool GetCuts(TString region, TString cut, std::vector<KLepton> lep, std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets, std::vector<snu::KFatJet> fatjets, double MET, double LT, double HT, double ST, bool Is_SR, bool is_Tchannel);
  TString GetCuts_name(TString region, int cut, bool Is_SR);
  void DrawHistograms(TString region, std::vector<KLepton> lep, std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets, std::vector<snu::KJet> bjetsloose, std::vector<snu::KFatJet> fatjets, double MET, double LT, double HT, double ST, bool Draw_SR, bool Draw_CR, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, bool is_Tchannel);
  double GetWeight(bool geterr, TString region, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons);
  void FillLeptonHist(TString hist_prefix, TString hist_suffix, KLepton this_lep, double this_weight);
  std::vector<snu::KElectron> ShiftElectronEnergy(std::vector<snu::KElectron> beforeshift, TString el_ID, bool applyshift);
  void GENSignalStudy( bool Is_Signal );
  void GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index );

  std::vector<TString> triggerlist_mm;
  std::vector<TString> triggerlist_ee;
  std::vector<TString> triggerlist_em1;
  std::vector<TString> triggerlist_em2;
  double weight_err;


  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HNOSDiLepton_Ntp, 1);
};
#endif
