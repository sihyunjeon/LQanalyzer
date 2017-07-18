// $Id: HNOSDiLepton_Ntp.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNOSDiLepton_Ntp Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNOSDiLepton_Ntp.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNOSDiLepton_Ntp);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNOSDiLepton_Ntp::HNOSDiLepton_Ntp() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNOSDiLepton_Ntp");
  
  Message("In HNOSDiLepton_Ntp constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
cout<<"nowrong"<<endl;

}


void HNOSDiLepton_Ntp::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
 cout<<"nowrong"<<endl; 
  Message("Making clever hists for Z ->ll test code", INFO);
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:

}


void HNOSDiLepton_Ntp::ExecuteEvents()throw( LQError ){
cout<<"what'swrong"<<endl;
  // ========== Apply the gen weight ====================
  if(!isData) weight*=MCweight;
  // ================================================================================

  // ========== No cut ====================
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  FillCutFlow("NoCut", weight);
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  // ================================================================================


  // ========== MET filter cut ====================
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", weight);
  // ================================================================================


  // ========== Primary vertex cut ====================
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex      
  // ================================================================================

  triggerlist_mm.clear(); triggerlist_em1.clear(); triggerlist_em2.clear(); triggerlist_ee.clear();
  // ========== Trigger cut ====================
  triggerlist_mm.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_mm.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  triggerlist_em1.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_em2.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
//  triggerlist_em2.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");

//  triggerlist_me1.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
//  triggerlist_me2.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_ee.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  if(!(PassTriggerOR(triggerlist_mm)||PassTriggerOR(triggerlist_em1)||PassTriggerOR(triggerlist_em2)||PassTriggerOR(triggerlist_ee))) return;
  bool Pass_Trigger_mm = PassTriggerOR(triggerlist_mm);
  bool Pass_Trigger_em = (PassTriggerOR(triggerlist_em1) || PassTriggerOR(triggerlist_em2));
  bool Pass_Trigger_ee = PassTriggerOR(triggerlist_ee);

  if(k_channel.Contains("DoubleMuon")) if(!Pass_Trigger_mm && (Pass_Trigger_ee||Pass_Trigger_em)) return;
  if(k_channel.Contains("DoubleEG")) if(!Pass_Trigger_ee && (Pass_Trigger_mm||Pass_Trigger_em)) return;
  if(k_channel.Contains("MuonEG")) if(!Pass_Trigger_em && (Pass_Trigger_ee||Pass_Trigger_mm)) return;
  // ================================================================================

  // ========== Get Objects (muon, electron, jet) ====================
  std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_HN_VETO", false);
  std::vector<snu::KMuon> muons;
  std::vector<snu::KMuon> muonTightColl;
  muonTightColl.clear(); muons.clear();
  int muonVetoN = muonVetoColl.size(), muonTightN = 0, muonsN = 0;
  for(unsigned int i=0;i<muonVetoN;i++){
    if(PassID(muonVetoColl.at(i), "MUON_HN_LOOSE")){
      muons.push_back(muonVetoColl.at(i));
      muonsN ++;
      if(PassID(muonVetoColl.at(i), "MUON_HN_TIGHT")){
        muonTightColl.push_back(muonVetoColl.at(i));
        muonTightN++;
      }
    }
  }

  std::vector<snu::KElectron> electronVetoColl = GetElectrons(false,false,"ELECTRON_HN_VETO");
  std::vector<snu::KElectron> electrons;
  std::vector<snu::KElectron> electronTightColl;
  electronTightColl.clear(); electrons.clear();
  int electronVetoN = electronVetoColl.size(), electronTightN = 0, electronsN = 0;
  for(unsigned int i=0;i<electronVetoN;i++){
    if(PassID(electronVetoColl.at(i), "ELECTRON_HN_FAKELOOSE")){
      electrons.push_back(electronVetoColl.at(i));
      electronsN ++;
      if(PassID(electronVetoColl.at(i), "ELECTRON_HN_TIGHTv4")){
        electronTightColl.push_back(electronVetoColl.at(i));
        electronTightN++;
      }
    }
  }
cout<< electronTightN << muonTightN <<endl;
  if(!k_running_nonprompt){
    if(!(muonsN == muonTightN && electronsN == electronTightN)) return;
    if(!(muonsN == muonVetoN && electronsN == electronVetoN)) return;
  }
  if(k_running_nonprompt){
    if(muonsN == muonTightN && electronsN == electronTightN) return;
    if(!(muonsN == muonVetoN && electronsN == electronVetoN)) return;
  }
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  std::vector<snu::KFatJet> fatjets = GetFatJets("FATJET_HN");
//  std::vector<snu::KFatJet> fatjets = GetFatJets("FATJET_NOCUT");
  std::vector<snu::KJet> jets_for_bjet = GetJets("JET_NOLEPTONVETO");
  std::vector<snu::KJet> Tjets = GetJets("JET_HN_TChannel", 20.);
  std::vector<snu::KJet> bjets, bjetsloose;
  std::vector<snu::KJet> frontTjets;
  std::vector<snu::KJet> backTjets;
  bjets.clear(); bjetsloose.clear(); frontTjets.clear(); backTjets.clear();
  for(unsigned int j=0; j<jets_for_bjet.size(); j++){
    if(jets_for_bjet.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Loose)){
      bjetsloose.push_back(jets_for_bjet.at(j));
      if(jets_for_bjet.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) bjets.push_back(jets_for_bjet.at(j));
    }
  }
  for(unsigned int j=0; j<Tjets.size(); j++){
    if(Tjets.at(j).Eta() > 2.5) frontTjets.push_back(Tjets.at(j));
    if(Tjets.at(j).Eta() < -2.5) backTjets.push_back(Tjets.at(j));
  }
  bool is_Tchannel = (frontTjets.size() != 0 && backTjets.size() != 0);

  double MET = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  // ================================================================================

  // ========== Momentum Correction ===================
  if(k_running_chargeflip){
    double old_sum=0., new_sum=0.;
    for(int i=0;i<electronsN; i++) old_sum+=electrons.at(i).Pt();

    electrons = ShiftElectronEnergy(electronTightColl, "ELECTRON_HN_TIGHTv4", true);
    for(int i=0;i<electronsN; i++) new_sum+=electrons.at(i).Pt();

    MET += fabs(new_sum-old_sum);
  }

  CorrectMuonMomentum(muons);
  MET = CorrectedMETRochester(muons, true);
  // ================================================================================

  TString LeptonConfig_mm = "NULL", ChargeConfig_mm = "NULL";
  bool Pass_Pt_mm = false;
  if(Pass_Trigger_mm ){
    if(muonsN == 2 && electronsN == 0){

      LeptonConfig_mm = "2mu0el";

      Pass_Pt_mm = (muons.at(0).Pt() > 20 && muons.at(1).Pt() >10);

      if(muons.at(0).Charge() == muons.at(1).Charge()) ChargeConfig_mm = "SSSF";
      else ChargeConfig_mm = "OSSF";

    }
  }

  TString LeptonConfig_ee = "NULL", ChargeConfig_ee = "NULL";
  bool Pass_Pt_ee = false;
  if(Pass_Trigger_ee){
    if(muonsN == 0 && electronsN == 2){

       LeptonConfig_ee = "0mu2el";

       Pass_Pt_ee = (electrons.at(0).Pt() > 25 && electrons.at(1).Pt() >15);

      if(electrons.at(0).Charge() == electrons.at(1).Charge()) ChargeConfig_ee = "SSSF";
      else ChargeConfig_ee = "OSSF";

    }
  }

  TString LeptonConfig_em = "NULL", ChargeConfig_em = "NULL";
  bool Pass_Pt_em = false;
  if(Pass_Trigger_em){
    if(muonsN == 1 && electronsN == 1){

      LeptonConfig_em = "1mu1el";
 
      Pass_Pt_em = (muons.at(0).Pt() > 10 && electrons.at(0).Pt() > 25);

      if(muons.at(0).Charge() != electrons.at(0).Charge()) ChargeConfig_em = "OSOF";
      else ChargeConfig_em = "SSOF";

    }
  }
  if(LeptonConfig_mm == "NULL" && LeptonConfig_ee == "NULL" && LeptonConfig_em == "NULL") return;
  std::vector<KLepton> lep;
cout<<"here1"<<endl;
  if(muonsN==1 && electronsN==1){
    if(muons.at(0).Pt() > electrons.at(0).Pt()){
      lep.push_back(muons.at(0));
      lep.push_back(electrons.at(0));
    }
    else{
      lep.push_back(electrons.at(0));
      lep.push_back(muons.at(0));
    }
  }
  else{
    for(int i=0; i<muonsN; i++){
      lep.push_back(muons.at(i));
    }
    for(int i=0; i<electronsN; i++){
      lep.push_back(electrons.at(i));
    }
  }

  TString region = "NULL";
  if(Pass_Trigger_mm && LeptonConfig_mm == "2mu0el" && Pass_Pt_mm && ChargeConfig_mm == "SSSF") region = "DiMu_SS";
  if(Pass_Trigger_mm && LeptonConfig_mm == "2mu0el" && Pass_Pt_mm && ChargeConfig_mm == "OSSF") region = "DiMu_OS";
  if(Pass_Trigger_ee && LeptonConfig_ee == "0mu2el" && Pass_Pt_ee && ChargeConfig_ee == "SSSF") region = "DiEl_SS";
  if(Pass_Trigger_ee && LeptonConfig_ee == "0mu2el" && Pass_Pt_ee && ChargeConfig_ee == "OSSF") region = "DiEl_OS";
  if(Pass_Trigger_em && LeptonConfig_em == "1mu1el" && Pass_Pt_em && ChargeConfig_em == "SSOF") region = "MuEl_SS";
  if(Pass_Trigger_em && LeptonConfig_em == "1mu1el" && Pass_Pt_em && ChargeConfig_em == "OSOF") region = "MuEl_OS";

  if(region == "NULL") return;

  if(k_running_chargeflip){
    if(region != "DiEl_OS" && region != "MuEl_OS"){
      return;
    }
  }

  double HT = 0.;
  for(unsigned int i=0; i<jets.size(); i++){
    HT += jets.at(i).Pt();
  }
  double LT = 0.;
  for(unsigned int i=0; i<lep.size(); i++){
    LT += lep.at(i).Pt();
  }
  double ST = LT + HT;


  Pass_Preselection = false;
  DrawHistograms(region, lep, jets, bjets, bjetsloose, fatjets, MET, LT, HT, ST, true, true, muons, electrons, is_Tchannel);

   return;
}// End of execute event loop

void HNOSDiLepton_Ntp::DrawHistograms(TString region, std::vector<KLepton> lep, std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets, std::vector<snu::KJet> bjetsloose, std::vector<snu::KFatJet> fatjets, double MET, double LT, double HT, double ST, bool Draw_SR, bool Draw_CR, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, bool is_Tchannel){

  if(jets.size() < 2 && fatjets.size() < 1) return;
  if(bjets.size() != 0) return;
  if(region != "MuEl_OS" && region != "MuEl_SS") return;
cout<<"HERE?2"<<endl;
  int region_index = 999;
  if(region == "MuEl_OS") region_index=1;
  if(region == "MuEl_SS") region_index =2; 

  if(k_running_chargeflip) region_index=2;

  double temp_weight = GetWeight(false, region, muons, electrons);
  double temp_weight_err = GetWeight(true, region, muons, electrons);

  double leading_pt(0.), subleading_pt(0.);
  snu::KParticle W_candidate, W_selection;
  int jetselection[2] = {0,1};
  if(jets.size() > 1){
    W_selection = jets.at(0)+jets.at(1);
    for(unsigned int i=0; i<jets.size(); i++){
      for(unsigned int j=i+1; j<jets.size(); j++){
        W_candidate = jets.at(i)+jets.at(j);
        if(fabs(W_candidate.M() - 80.4) < fabs(W_selection.M() - 80.4)){
          W_selection = W_candidate;
          jetselection[0] = i;
          jetselection[1] = j;
        }
      }
    }
  }
cout<<"here3"<<endl;
  double cutop[100] = {-999.,};
  leading_pt = lep.at(0).Pt();
  subleading_pt = lep.at(1).Pt();

  cutop[0] = region_index;
  cutop[1] = temp_weight;
  cutop[2] = temp_weight_err;

  cutop[3] = (lep.at(0).Pt());
  cutop[4] = (lep.at(1).Pt());
  cutop[5] = MET;
  cutop[6] = LT;
  cutop[7] = ST;
  cutop[8] = HT;
  cutop[9] = (MET*MET/ST);
  cutop[10] = (lep.at(0)+lep.at(1)).M();
  cutop[11] = (lep.at(0).DeltaR(lep.at(1)));
  

  if(jets.size() > 1){
    cutop[12] = (jets.size());
    cutop[13] = (jets.at((jetselection[0])).Pt());
    cutop[14] = (jets.at((jetselection[1])).Pt());
    cutop[15] = (jets.at(jetselection[0]).DeltaR(jets.at(jetselection[1])));
    cutop[16] = (W_selection.M());
    cutop[17] = ((lep.at(0)+W_selection).M());
    cutop[18] = ((lep.at(1)+W_selection).M());
  }

  if(fatjets.size() > 0){
    cutop[19] = (fatjets.size());
    cutop[20] = (fatjets.at(0).M());
    cutop[21] = (fatjets.at(0).Pt());
    cutop[22] = (fatjets.at(0).PrunedMass());
    cutop[23] = ((lep.at(0)+fatjets.at(0)).M());
    cutop[24] = ((lep.at(1)+fatjets.at(0)).M());
    cutop[25] = (lep.at(0).DeltaR(fatjets.at(0)));
    cutop[26] = (lep.at(1).DeltaR(fatjets.at(0)));
  }  


  FillNtp("Ntp_Preselection",cutop);

  return;

}
 
double HNOSDiLepton_Ntp::GetWeight(bool geterr, TString region, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons){

  bool flip = false;
  if(std::find(k_flags.begin(), k_flags.end(), "flip") !=k_flags.end()) flip = true;

  if(!k_running_nonprompt && !flip){
    double muon_id_iso_sf = 1.;
    double MuTrkEffSF = 1.;
    double trigger_sf = 1.;
    double pileup_reweight = 1.;
    double electron_idsf = 1.;
    double electron_recosf = 1.;
    if(!k_isdata){
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TIGHT", muons, 0);
      MuTrkEffSF = mcdata_correction->MuonTrackingEffScaleFactor(muons);
      double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "ELECTRON_HN_TIGHTv4", muons, "MUON_HN_TIGHT", 0, 0, 0);
      double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "ELECTRON_HN_TIGHTv4", muons, "MUON_HN_TIGHT", 0, 1, 0);
      trigger_sf = trigger_eff_Data/trigger_eff_MC;
      pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
      electron_idsf = mcdata_correction->ElectronScaleFactor("ELECTRON_HN_TIGHTv4", electrons);
      electron_recosf = mcdata_correction->ElectronRecoScaleFactor(electrons);
    }

    double trigger_ps_weight = 0.;
    if(region == "DiMu_OS" || region == "DiMu_SS"){
      trigger_ps_weight = WeightByTrigger(triggerlist_mm, TargetLumi);
    }
    if(region == "DiEl_OS" || region == "DiEl_SS"){
      trigger_ps_weight = WeightByTrigger(triggerlist_ee, TargetLumi);
    }
    if(region == "MuEl_OS"){
      if(PassTriggerOR(triggerlist_em1)) trigger_ps_weight += WeightByTrigger(triggerlist_em1, TargetLumi);
      if(PassTriggerOR(triggerlist_em2)) trigger_ps_weight += WeightByTrigger(triggerlist_em2, TargetLumi);
    }

    double cutflow_weight = weight;
    if(!k_isdata){
      weight *= muon_id_iso_sf;
      weight *= MuTrkEffSF;
      weight *= trigger_sf;
      weight *= pileup_reweight;
      weight *= electron_idsf;
      weight *= electron_recosf;
      weight *= trigger_ps_weight;
      weight_err = 0.;
    }
    if(std::find(k_flags.begin(), k_flags.end(), "cutflow") !=k_flags.end()){
      cutflow_weight *= trigger_ps_weight;
      cutflow_weight *= pileup_reweight;
      weight = cutflow_weight;
    }

    if(k_isdata){ weight = 1.; weight_err = 0.; }
  }
  if(k_running_nonprompt && !flip){
    weight = m_datadriven_bkg->Get_DataDrivenWeight(false,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");
    weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");
  }

  if(!k_running_nonprompt && flip){
    weight = GetCFweight(electrons, true, "ELECTRON_HN_TIGHTv4");
    weight_err = 0.;
  }

  if(k_running_nonprompt && flip){
    weight = m_datadriven_bkg->Get_DataDrivenWeight(false,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");
    weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSE", "mva");
    weight *= GetCFweight(electrons, true, "ELECTRON_HN_TIGHTv4");
    weight = weight*(-1.);    
  }


  if(geterr) return weight_err;
  if(!geterr) return weight;
}

void HNOSDiLepton_Ntp::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNOSDiLepton_Ntp::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

HNOSDiLepton_Ntp::~HNOSDiLepton_Ntp() {
  
  Message("In HNOSDiLepton_Ntp Destructor" , INFO);
  
}


void HNOSDiLepton_Ntp::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNOSDiLepton_Ntp::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  MakeNtp("Ntp_Preselection", "chargeconfig:weight:weight_err:lep1pt:lep2pt:MET:lt:st:ht:metsqdivst:lep1lep2mass:lep1lep2deltar:jetn:jet1pt:jet2pt:jet1jet2deltar:jet1jet2mass:lep1jet1jet2mass:lep2jet1jet2mass:fatjetn:fatjetmass:fatjetpt:fatjetprunedmass:lep1fatjetmass:lep2fatjetmass:lep1deltar:lep2fatjetdeltar");

  /**
   *  Remove//Overide this HNOSDiLepton_NtpCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNOSDiLepton_Ntp::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


std::vector<snu::KElectron> HNOSDiLepton_Ntp::ShiftElectronEnergy(std::vector<snu::KElectron> beforeshift, TString el_ID, bool applyshift){

  if(el_ID != "ELECTRON_HN_TIGHTv4") return beforeshift;
  if(!applyshift) return beforeshift;

  std::vector<snu::KElectron> aftershift;
  double shiftrate = -999.;
  if(beforeshift.size() == 1) shiftrate = (1-0.024);
  if(beforeshift.size() == 2) shiftrate = (1-0.011);
  if(beforeshift.size() > 2) shiftrate = (-999.);


   for(unsigned int i=0; i < beforeshift.size(); i++){
     beforeshift.at(i).SetPtEtaPhiM(beforeshift.at(i).Pt()*shiftrate, beforeshift.at(i).Eta(), beforeshift.at(i).Phi(), 0.511e-3);
     aftershift.push_back(beforeshift.at(i));
   }
   return aftershift;
}
