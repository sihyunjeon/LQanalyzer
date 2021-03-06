#include "MuonSelection.h"

using namespace snu;

MuonSelection::MuonSelection(LQEvent ev) :
  BaseSelection()
{
  k_lqevent = ev;  

};


MuonSelection::~MuonSelection() {};

void MuonSelection::BasicSelection( std::vector<KMuon>& leptonColl, bool m_debug) {

  std::vector<KMuon> allmuons = k_lqevent.GetMuons();
  int ilep(0);
  for (std::vector<KMuon>::iterator muit = allmuons.begin(); muit!=allmuons.end(); muit++, ilep++)
    {

      if(muit->Pt() == 0.) continue;

      bool pass_selection = true;

      /// ONLY CUT ON PT/ETA/LOOSE ID
      if( muit->Pt() < pt_cut_min ) {
	pass_selection = false;
	if(m_debug) cout << "BasicSelection:: Muon Fails Pt cut " << endl; 
      }
      if( fabs(muit->Eta()) > eta_cut){
	pass_selection =false;
	if(m_debug) cout << "BasicSelection:: Muon Fails Eta cut " << endl; 
      }
      if(! (PassID("MUON_POG_LOOSE", *muit, m_debug))){
	pass_selection =false;
	if(m_debug) cout << "BasicSelection:: Muon Fails Loose Selection" << endl;
      }
      
      if(pass_selection) leptonColl.push_back(*muit);
    }
  if(m_debug) cout << "BasicSelection::Number of muons = " << leptonColl.size() << endl;
  return;
}


void MuonSelection::SkimSelection( std::vector<KMuon>& leptonColl , bool m_debug) {

  std::vector<KMuon> allmuons = k_lqevent.GetMuons();
  int ilep(0);
  for (std::vector<KMuon>::iterator muit = allmuons.begin(); muit!=allmuons.end(); muit++, ilep++)
    {
      
      if(muit->Pt() == 0.) continue;

      bool pass_selection = true;
      
      /// ONLY CUT ON PT/ETA/LOOSE ID
      if( muit->Pt() < pt_cut_min ) {
        pass_selection = false;
        if(m_debug) cout << "SkimSelection:: Muon Fails Pt cut " << endl;
      }
      if( fabs(muit->Eta()) > eta_cut){
        pass_selection =false;
        if(m_debug) cout << "SkimSelection:: Muon Fails Eta cut " << endl;
      }
      if(! (PassID("MUON_POG_LOOSE", *muit, m_debug))){
        pass_selection =false;
        if(m_debug) cout << "SkimSelection:: Muon Fails Loose Selection" << endl;
      }

      if (pass_selection) leptonColl.push_back(*muit);  
    }
  if(m_debug) cout << "SkimSelection::Number of muons = " << leptonColl.size() << endl;
  return;

}

void MuonSelection::Selection( std::vector<KMuon>& leptonColl, bool m_debug) {
  

  std::vector<KMuon> allmuons = k_lqevent.GetMuons();
  int ilep(0);
  for (std::vector<KMuon>::iterator muit = allmuons.begin(); muit!=allmuons.end(); muit++, ilep++){ 

      bool pass_selection(true);      
      if(muit->Pt() == 0.) continue;

      TString MuID=GetString(k_id);
      if(apply_ID && !PassID(MuID, *muit)) pass_selection =false;
      if(m_debug&& apply_ID && !PassID(MuID, *muit)) cout << "Fails Selection::ID cut " << endl;


      if(apply_ptcut && ! ( muit->Pt() > pt_cut_min )) pass_selection = false;
      if(m_debug&&apply_ptcut && ! (muit->Pt() >= pt_cut_min && muit->Pt() < pt_cut_max)) cout << "Fails Selection::pt cut " << endl;
      
      if(apply_etacut && !(fabs(muit->Eta()) < eta_cut)) pass_selection =false;
      if(m_debug&&apply_etacut && !(fabs(muit->Eta()) < eta_cut))  cout << "Fails Selection::eta cut " << endl;


      //// Calculate PF isolation
      /// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
      //LeptonRelIso = (muit->RelIso04());
      if     (apply_relisocut && RelIsoType.Contains("Default"))    LeptonRelIso=muit->RelIso04();
      else if(apply_relisocut && RelIsoType.Contains("PFRelIso03")) LeptonRelIso=muit->RelIso03();
      else   LeptonRelIso=muit->RelIso04();
      
      if(apply_relisocut && !( LeptonRelIso < relIso_cut)) pass_selection = false;
      if(m_debug&&apply_relisocut && !( LeptonRelIso < relIso_cut))  cout << "Fails Selection::reliso cut " << endl;
      
  
      /// impact parameter cuts
      // Uses fabs(recoMu.muonBestTrack()->dxy(vertex->position())) as described in https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon 
      // Also stores dB() on pat::Muon  as dxy_pat
      // Also stores D0 
      if(apply_dzcut && !(fabs(muit->dZ())<  dz_cut )) pass_selection = false;
      if(m_debug&&apply_dzcut && !(fabs(muit->dZ())<  dz_cut ))  cout << "Fails Selection::dz cut " << endl;
      if(apply_dxycut && !(fabs(muit->dXY())< dxy_cut )) pass_selection = false;
      if(m_debug&&apply_dxycut && !(fabs(muit->dXY())< dxy_cut ))cout << "Fails Selection::dxy cut " << endl;
      
      
      if(apply_chi2cut && !( muit->GlobalChi2() < chiNdof_cut && muit->GlobalChi2() >= chiNdofMIN_cut )) pass_selection = false;
      if(m_debug&&apply_chi2cut && !( muit->GlobalChi2() <chiNdof_cut && muit->GlobalChi2()  >=chiNdofMIN_cut)) cout << "Fails chi2 cut " << endl;

      
      //// ADD EXTRA  cut on D0sig? or same vertex?     
      if(pass_selection) leptonColl.push_back(*muit);    

  }/// muon loop end

  BaseSelection::reset();
  return;
}  



void MuonSelection::SelectMuons(std::vector<KMuon>& leptonColl, ID muid, float ptcut, float etacut){
  
  return SelectMuons(leptonColl, GetString(muid), ptcut, etacut);
}

void MuonSelection::SelectMuons(std::vector<KMuon>& leptonColl, TString muid, float ptcut, float etacut){
  std::vector<KMuon> allmuons = k_lqevent.GetMuons();

  if (ptcut == -999.) ptcut = AccessFloatMap("ptmin",muid);
  if (etacut == -999.) etacut = AccessFloatMap("|etamax|",muid);
  //cout << "cuts " << ptcut  << " " << etacut << endl;
  for (std::vector<KMuon>::iterator muit = allmuons.begin(); muit!=allmuons.end(); muit++){

    bool pass_selection(true);
    if(muit->Pt() == 0.)   continue;

    MuonID = PassUserID(muid, *muit);
    if(!MuonID)  pass_selection = false;

 
   if(( muit->Pt() < ptcut )) pass_selection = false;
    if(!(fabs(muit->Eta()) < etacut)) pass_selection = false;
    if(pass_selection)  leptonColl.push_back(*muit);
  }
  return;
}


bool MuonSelection::PassUserID(TString id, snu::KMuon mu){

  if(mu.Pt() == 0.) return false;

  float isomax = AccessFloatMap("isomax04",id);
  //float isomin = AccessFloatMap("isomin04",id);
  float dxymax = AccessFloatMap("|dxymax|",id);
  //float dxymin = AccessFloatMap("|dxymin|",id);
  float dxysigmax = AccessFloatMap("|dxysigmax|",id);
  float dxysigmin = AccessFloatMap("|dxysigmin|",id);

  float dzmax = AccessFloatMap("|dzmax|",id);
  //float dzmin = AccessFloatMap("|dzmin|",id);
  float chi2max = AccessFloatMap("chi2max",id);
  //float chi2min = AccessFloatMap("chi2min",id);

  bool checkisomax     = CheckCutFloat("isomax04",id);
  //bool checkisomin     = CheckCutFloat("isomin04",id);
  bool checkdxymax      = CheckCutFloat("|dxymax|",id);
  //bool checkdxymin      = CheckCutFloat("|dxymin|",id);
  bool checkdzmax       = CheckCutFloat("|dzmax|",id);
  bool checkdxysigmin  = CheckCutFloat("|dxysigmin|",id);
  bool checkdxysigmax  = CheckCutFloat("|dxysigmax|",id);
  //bool checkdzmin       = CheckCutFloat("|dzmin|",id);
  bool checkchi2max     = CheckCutFloat("chi2max",id);
  //bool checkchi2min     = CheckCutFloat("chi2min",id);
  bool checkisloose  = (CheckCutString("IsLoose(POG)",id));
  bool checkismedium = (CheckCutString("IsMedium(POG)",id));
  bool checkistight  = (CheckCutString("IsTight(POG)",id));
  

  bool debug=false;
  //if(id.Contains("VETO")) debug=true;
  LeptonRelIso = (mu.RelIso04());

  bool pass_selection=true;
  if(checkisloose && ! mu.IsLoose ()) { pass_selection = false;if(debug){ cout << "Fail isloose" << endl;}}
  
  if(checkismedium && ! mu.IsMedium ()) { pass_selection = false;if(debug){ cout << "Fail ismedium"  << endl;}}
  if(checkdxymax || checkchi2max || checkdzmax) {
    if(checkistight && ! PassID("MUON_POG_TIGHT",mu, !checkdxymax,!checkdzmax,!checkchi2max)) { pass_selection = false;if(debug){ cout << "Fail pogtight"  << endl;}}
    if(checkdxymax && (fabs(mu.dXY()) > dxymax)){ pass_selection = false;if(debug){ cout << "Fail dxy "  << endl;}}
    if(checkdzmax && (fabs(mu.dZ()) > dzmax)){ pass_selection = false;if(debug){ cout << "Fail dZ"  << endl;}}
    if(checkchi2max && (fabs(mu.GlobalChi2()) > chi2max)){ pass_selection = false;if(debug){ cout << "Fail GlobalChi2"  << endl;}}
  }
  else  if(checkistight &&  ! mu.IsTight ()) { pass_selection = false;if(debug){ cout << "Fail tight"  << endl;}}


  if(checkisomax && (LeptonRelIso > isomax)) { pass_selection = false;if(debug){ cout << "Fail iso"  << endl;}}
  //if(checkisomin && (LeptonRelIso < isomin)) { pass_selection = false;if(debug){ cout << "Fail "  << endl;}}

  if(checkdxysigmin &&(fabs(mu.dXYSig()) < dxysigmin)) { pass_selection = false;if(debug){ cout << "Fail dsximin"  << endl;}}
  if(checkdxysigmax &&(fabs(mu.dXYSig()) > dxysigmax)) { pass_selection = false;if(debug){ cout << "Fail dsigmax"  << endl;}}
  //if(checkdxymin && (fabs(mu.dXY()) < dxymin)){ pass_selection = false;if(debug){ cout << "Fail "  << endl;}}

  //if(checkdzmin && (fabs(mu.dZ()) < dzmin)){ pass_selection = false;if(debug){ cout << "Fail "  << endl;}}

  return pass_selection;
  
}


bool MuonSelection::HNTightMuonSelection(KMuon mu) {
  
  if(HNIsTight(mu))  return true;
  else return false;
}



////////// PREDEFINED MUON SELECTIONS FOR TOP ANALYSIS


bool MuonSelection::TopVetoMuonSelection(KMuon mu) {
  
  bool pass_selection(true);
  LeptonRelIso = (mu.RelIso03());
  
  if(!( LeptonRelIso < 0.2))     pass_selection = false;
  if(!PassID("MUON_POG_LOOSE", mu))      pass_selection =false;
  //// Make Loose selection
  
  return pass_selection;

}

// Intended use for di-lepton channel
bool MuonSelection::TopLooseMuonSelection(KMuon mu) {

  bool pass_selection(true);
  
  if(mu.Pt() == 0.) return false;
  LeptonRelIso = (mu.RelIso03());

  if(!PassID("MUON_POG_LOOSE", mu)) pass_selection =false;
  if(!( LeptonRelIso < 0.2)) pass_selection = false;

  
  return pass_selection;

}

bool MuonSelection::TopTightMuonSelection(KMuon mu) {

  bool pass_selection(true);
  
  if(mu.Pt() == 0.) return false;
  LeptonRelIso = (mu.RelIso03());
  
  if(!( LeptonRelIso < 0.1))    pass_selection = false;
  if(!(fabs(mu.dZ())< 0.10  ))  pass_selection = false;
  if(!(fabs(mu.dXY())< 0.005 ))   pass_selection = false;
  
  /// TIGHT MUON from muon POG
  if(!PassID("MUON_POG_TIGHT", mu)) pass_selection =false;
  if(mu.Pt() == 0.) return false;
  return pass_selection;

}


/// NO LONGER NEEDED
bool MuonSelection::PassID(TString id, snu::KMuon mu, bool cutondxy, bool cutondz, bool cutonchi2, bool m_debug){


  /// Taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonIdRun2
  bool passID(true);
  if (id == "MUON_POG_LOOSE") {
    if(!(mu.IsPF() == 1)) {
      passID = false;
      if(m_debug)cout << "PassID: Fail isPF" << endl;
    }
    if(!(mu.IsGlobal()==1 || mu.IsTracker() == 1 )){
      passID = false;
      if(m_debug){
        cout << "PassID: Fail isGlobal||isTracker" << endl;
        cout << "PassID: mu.IsGlobal()=  " << mu.IsGlobal() << endl;
        cout << "PassID: mu.IsTracker()= " << mu.IsTracker() << endl;
      }
    }
  }


  else if (id == "MUON_POG_TIGHT") {
    if(!(mu.IsPF() == 1        )){
      passID = false;
      if(m_debug)cout << "PassID: Fail isPF" << endl;
    }
    if(!(mu.IsGlobal()==1      )) {
      passID = false;
      if(m_debug)cout << "PassID: Fail isGlobal" << endl;
    }
    if( mu.validHits() == 0     ) {
      passID = false;
      if(m_debug)cout << "PassID: Fail validHit" << endl;
    }
    if( mu.validPixHits() == 0) {
      passID = false;
      if(m_debug)cout << "PassID: Fail validPixelHit" << endl;
    }
    if( mu.validStations() <= 1 ) {
      passID = false;
      if(m_debug)cout << "PassID: Fail validStations" << endl;
    }
    if( mu.ActiveLayer() <= 5   ) {
      passID = false;
      if(m_debug)cout << "PassID: Fail ActiveLayer " << endl;
    }
    if (cutondxy){
      if( fabs(mu.dXY())    >= 0.2) {
        passID = false;
        if(m_debug)cout << "PassID: Fail dXY" << endl;
      }
    }
    if(cutondz){
      if( fabs(mu.dZ())    >= 0.5) {
        passID = false;
        if(m_debug)cout << "PassID: Fail dZ" << endl;
      }
    }
    if(cutonchi2){
      if( mu.GlobalChi2() >=  10.){
        passID = false;
        if(m_debug) cout << "PassID: Fail  Chi2" << endl;
      }
    }
  }

  else{
    cout << "Invalid ID set for muon selection" << endl;
  }
  return passID;
}




void MuonSelection::SetDeposits(Double_t ECalDeposit , Double_t HCalDeposit) {
    ECalDeposit ? ECalDeposit_max = ECalDeposit : ECalDeposit_max=4.0;
    HCalDeposit ? HCalDeposit_max = HCalDeposit : HCalDeposit_max=6.0;
    ECalDeposit_min = 0.0;
    HCalDeposit_min = 0.0;
}

void MuonSelection::SetDeposits(Double_t ECalDeposit1 , Double_t HCalDeposit1, Double_t ECalDeposit2 , Double_t HCalDeposit2) {
    ECalDeposit1 ? ECalDeposit_min = ECalDeposit1 : ECalDeposit_min=0.0;
    HCalDeposit1 ? HCalDeposit_min = HCalDeposit1 : HCalDeposit_min=0.0;
    ECalDeposit2 ? ECalDeposit_max = ECalDeposit2 : ECalDeposit_max=4.0;
    HCalDeposit2 ? HCalDeposit_max = HCalDeposit2 : HCalDeposit_max=6.0;
}


MuonSelection& MuonSelection::operator= (const MuonSelection& ms) {
  if(this != &ms){    
    BaseSelection::operator = (ms);
    k_lqevent = ms.k_lqevent;  

  }
  return *this;
};

MuonSelection::MuonSelection(const MuonSelection& ms):
  BaseSelection(ms)
{
  k_lqevent = ms.k_lqevent; 
};

