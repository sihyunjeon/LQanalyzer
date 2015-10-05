#ifndef AnalyzerCore_H
#define AnalyzerCore_H

//forward declarations                                                                                                                                            
class Reweight;
class EventBase;
class MuonPlots;
class ElectronPlots;
class JetPlots;
class SignalPlotsEE;
class SignalPlotsMM;
class SignalPlotsEM;
class TriLeptonPlots;
class EventBase;

#include "LQCycleBase.h"
#include "HNCommonLeptonFakes/HNCommonLeptonFakes/HNCommonLeptonFakes.h"

class AnalyzerCore : public LQCycleBase {
  
 public:
  
  // Default constructor
  AnalyzerCore();

  //destructor
  virtual ~AnalyzerCore();

  // SetUpEvent CORE function: accesses event in ntuple
  virtual void SetUpEvent(Long64_t entry, float ev_weight)throw( LQError );
  virtual void EndEvent()throw( LQError );
  virtual void WriteHistograms()throw( LQError );


  TDirectory*   getTemporaryDirectory(void) const;

  std::vector<snu::KJet>  GetJets(TString label);
  std::vector<snu::KMuon> GetMuons(TString label);
  std::vector<snu::KMuon> GetMuons(TString label, bool keepfakes);
  std::vector<snu::KElectron> GetElectrons(bool keepcf, bool keepfake, TString label);
  std::vector<snu::KElectron> GetElectrons( TString label);

  bool HasCloseBJet(snu::KElectron el);
  bool HasCloseLBJet(snu::KElectron el);


  float WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip);

<<<<<<< HEAD
  float WeightCFEvent(std::vector<snu::KElectron> electrons, bool runchargeflip, bool useoldrates);

  void GetIDEfficiency( std::vector<snu::KElectron> electronNoCutColl, std::vector<snu::KJet> jets, TString label, float w);
  void GetIDEfficiency( std::vector<snu::KElectron> electronNoCutColl, std::vector<snu::KMuon> muons,  std::vector<snu::KJet> jets, float w);
    
  void PlotFakeLeptons( std::vector<snu::KJet> jets, float w);
  double TriggerScaleFactor( vector<snu::KElectron> el, int sys=0);;
  double TriggerScaleFactor( vector<snu::KMuon> mu, int sys=0);;
  double TriggerScaleFactorEMu(int sys=0);
  float GetEEMass(std::vector<snu::KElectron> electrons);
=======
  double TriggerScaleFactor( vector<snu::KElectron> el);;
  double TriggerScaleFactor( vector<snu::KMuon> mu);;
  double TriggerScaleFactorEMu();
  float GetZMass(std::vector<snu::KMuon> muons);
  float GetZMass(std::vector<snu::KElectron> electrons);
>>>>>>> cat

  double ElectronScaleFactor( double eta, double pt, bool tight, int sys=0);
  double MuonScaleFactor(double eta, double pt, int sys=0);
  float  JetResCorr(snu::KJet jet, std::vector<snu::KGenJet> genjets);
  float SumPt( std::vector<snu::KJet> particles);
  bool isPrompt(long pdgid);
  bool IsTight(snu::KElectron electron);
  bool IsTight(snu::KMuon muon);
  std::vector<snu::KElectron> GetTruePrompt(vector<snu::KElectron> electrons,  bool keep_chargeflip, bool keepfake);
  std::vector<snu::KMuon> GetTruePrompt(vector<snu::KMuon> muons,   bool keepfake);
  int NBJet(std::vector<snu::KJet> jets);
  bool Zcandidate(vector<snu::KMuon> muons, float interval, bool require_os=true);
  bool Zcandidate(vector<snu::KElectron> electrons, float interval, bool require_os=true);
  bool SameCharge(std::vector<snu::KMuon> muons);
  bool SameCharge(std::vector<snu::KElectron> electrons);
  
  float CFRate(snu::KElectron el);
  std::vector<snu::KElectron>  ShiftElectronEnergy(std::vector<snu::KElectron> el, bool applyshift);

  float Get_DataDrivenWeight_EE(vector<snu::KElectron> k_electrons);
  float Get_DataDrivenWeight_MM(vector<snu::KMuon> k_muons);
  float Get_DataDrivenWeight_EM(vector<snu::KMuon> k_muons, vector<snu::KElectron> k_electrons);
 
  double MuonDYMassCorrection(std::vector<snu::KMuon> mu, double w);

  
  vector<TLorentzVector> MakeTLorentz( vector<snu::KElectron> el);
  vector<TLorentzVector> MakeTLorentz( vector<snu::KMuon> mu);
  vector<TLorentzVector> MakeTLorentz( vector<snu::KJet> jet);
  // enum for plotting functions/classes
  enum histtype {muhist, elhist, jethist, sighist_ee, sighist_mm, sighist_em, trilephist};
  
  
  //
  // Useful message function 
  //
  void Message(TString message, LQMsgType type=INFO);


  //
  //  Specify which triggers will be avaiable in KTrigger
  //
  void AddTriggerToList(TString triggername);
  
  /// Pileup Reweighting class
  static const Bool_t MC_pu = true;
  Reweight *reweightPU;

  //// Event base pointer. Used to get all objects for analysis
  EventBase* eventbase;
  
  UInt_t numberVertices;
  Bool_t *goodVerticiesB;

  TDirectory *Dir;
  map<TString, TH1*> maphist;
  map<TString, TH2*> maphist2D;
  TH2F* FRHist;
  TH2F* MuonSF;
  HNCommonLeptonFakes* m_fakeobj;
  
  /// Event weights
  Double_t MCweight, weight;

  // used to get trigger prescale
  Int_t prescale;
  
  std::vector<TString> triggerlist;

  //// Making cleaver hist maps
  map<TString, TriLeptonPlots*> mapCLhistTriLep;
  map<TString, SignalPlotsEE*> mapCLhistSigEE;
  map<TString, SignalPlotsMM*> mapCLhistSigMM;
  map<TString, SignalPlotsEM*> mapCLhistSigEM;
  map<TString, ElectronPlots*> mapCLhistEl;
  map<TString, MuonPlots*> mapCLhistMu;
  map<TString, JetPlots*> mapCLhistJet;
  
  //
  // Function that closes rootfile
  //
  void CloseFiles();
  

  //
  // Make Histograms and fill maphist
  //
  void MakeHistograms();
  void MakeHistograms(TString hname, int nbins, float xmin, float xmax);
  void MakeHistograms(TString hname, int nbins, float xbins[]);
  void MakeHistograms2D(TString hname, int nbinsx, float xbins[], int nbinsy, float ybins[]);
  void MakeHistograms2D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax);
    //
    // Makes temporary dir
    //
    TDirectory* GetTemporaryDirectory(void) const;                                                                                                                                 
  //
  // Checks if a file exists
  //
  void CheckFile(TFile* file) throw( LQError );
  
  //// Plotting 
  TH1* GetHist(TString hname);
  TH2* GetHist2D(TString hname);

  /// Fills hist in maphist
  void FillHist(TString histname, float value, float w );
  void FillHist(TString histname, float value, float w , float xmin, float xmax, int nbins=0);
  void FillHist(TString histname, float value, float w , float xmin[], int nbins=0);
  void FillHist(TString histname, float value1, float value2, float w , float x[], int nbinsx, float y[], int nbinsy);
  void FillHist(TString histname, float value1,  float value2, float w , float xmin, float xmax, int nbinsx,  float ymin, float ymax, int nbinsy);

  /// Fills clever hists
  void FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double weight);
  void FillCLHist(histtype type, TString hist, vector<snu::KMuon> muons , double weight);
  void FillCLHist(histtype type, TString hist, vector<snu::KElectron> electrons , double weight);
  void FillCLHist(histtype type, TString hist, vector<snu::KJet> jets , double weight);

  // Makes clever histograms
  void MakeCleverHistograms(histtype type, TString clhistname );

  /// File related                                                                                                                                                
  void OpenPutputFile();
  void WriteHists();
  void WriteCLHists();

  //// Event related                                                                                                                                              
  bool PassTrigger(std::vector<TString> list, int& prescale);
  bool PassBasicEventCuts();

};
#endif
