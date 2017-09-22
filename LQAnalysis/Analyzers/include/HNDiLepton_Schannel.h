#ifndef HNDiLepton_Schannel_h
#define HNDiLepton_Schannel_h

#include "AnalyzerCore.h"
class HNDiLepton_Schannel : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  HNDiLepton_Schannel();
  ~HNDiLepton_Schannel();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
 private:
  
  int jet_lowindex[2], jet_highindex[2];
  bool flip;
  bool Pass_Preselection, Pass_LowPreselection, Pass_HighPreselection;
  bool GetCuts(TString region, TString cut, std::vector<KLepton> lep, std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets, std::vector<snu::KFatJet> fatjets, double MET, double LT, double HT, double ST, bool Is_SR, bool is_Tchannel);
  TString GetCuts_name(TString region, int cut, bool Is_SR);
  void DrawHistograms(TString region, std::vector<KLepton> lep, std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets, std::vector<snu::KJet> bjetsloose, std::vector<snu::KFatJet> fatjets, double MET, double LT, double HT, double ST, bool Draw_SR, bool Draw_CR, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, bool is_Tchannel);
  double GetWeight(bool geterr, TString region, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons);
  void FillLeptonHist(TString hist_prefix, TString hist_suffix, KLepton this_lep, double this_weight);
  std::vector<snu::KElectron> ShiftElectronEnergy(std::vector<snu::KElectron> beforeshift, TString el_ID, bool applyshift);
  void GENSignalStudy( bool Is_Signal );
  void GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index );
  double MCT(snu::KJet jet1, snu::KJet jet2);
  bool PassEMuTriggerPt(std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons);

  std::vector<TString> triggerlist_mm;
  std::vector<TString> triggerlist_ee;
  std::vector<TString> triggerlist_emNoDZ1;
  std::vector<TString> triggerlist_emNoDZ2;
  std::vector<TString> triggerlist_emDZ1;
  std::vector<TString> triggerlist_emDZ2; 
  double weight_err;

  TH2D *hist_Muon_FR;
  double CorrPt(KLepton lep, double T_iso);
  double CorrPt(snu::KMuon lep, double T_iso);
  double CorrPt(snu::KElectron lep, double T_iso);
  double GetMuonFR(bool geterr, float pt, float eta);
  double GetMuonPR(bool geterr, float pt, float eta);

  double get_eventweight(bool geterr, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons);

  bool run_fake, run_cf;

  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( HNDiLepton_Schannel, 1);
};
#endif
