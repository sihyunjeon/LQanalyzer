// $Id: HNDiLeptonAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNDiLeptonAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNDiLeptonAnalyzer.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNDiLeptonAnalyzer);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNDiLeptonAnalyzer::HNDiLeptonAnalyzer() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNDiLeptonAnalyzer");
  
  Message("In HNDiLeptonAnalyzer constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNDiLeptonAnalyzer::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  Message("Making clever hists for Z ->ll test code", INFO);
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:

  TDirectory* origDir = gDirectory;
  TString lqdir = getenv("LQANALYZER_DIR");

/*  TString MuonFRType_Data = "v7_SIP3_";
  TFile *file_Muon_FR = new TFile( lqdir+"/data/Fake/80X/Muon_Data_v7_SIP3_FR.root");

  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet
    std::stringstream dirname;
    dirname << "HNCommonLeptonFakes_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }
  tempDir->cd();

  hist_Muon_FR = (TH2D*)file_Muon_FR->Get("Muon_Data_v7_SIP3_FR_Awayjet40")->Clone();

  file_Muon_FR->Close();
  delete file_Muon_FR;

  origDir->cd();
  return;*/
}


void HNDiLeptonAnalyzer::ExecuteEvents()throw( LQError ){

  jet_lowindex[0] = 0;  jet_lowindex[1] = 1; jet_highindex[0] = 0;  jet_highindex[1] = 1;
  Pass_Preselection = false; Pass_LowPreselection = false; Pass_HighPreselection = false;
  triggerlist_mm.clear(); triggerlist_emNoDZ1.clear(); triggerlist_emNoDZ2.clear(); triggerlist_emDZ1.clear(); triggerlist_emDZ2.clear(); triggerlist_ee.clear();

  run_fake = false;
  run_cf = false;

  if(std::find(k_flags.begin(), k_flags.end(), "fake") !=k_flags.end()) run_fake = true;
  if(std::find(k_flags.begin(), k_flags.end(), "cf") !=k_flags.end()) run_cf = true;


  // ========== Apply the gen weight ====================
  if(!isData) weight*=MCweight;
  // ================================================================================

  // ========== Trigger cut ====================
  triggerlist_mm.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_mm.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  triggerlist_emNoDZ1.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_emNoDZ2.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");

  triggerlist_emDZ1.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_emDZ2.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_ee.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  bool Pass_Trigger_mm = PassTriggerOR(triggerlist_mm);
  bool Pass_Trigger_em = (PassTriggerOR(triggerlist_emNoDZ1) || PassTriggerOR(triggerlist_emNoDZ2) || PassTriggerOR(triggerlist_emDZ1) || PassTriggerOR(triggerlist_emDZ2));
  bool Pass_Trigger_ee = PassTriggerOR(triggerlist_ee);
  if(!(Pass_Trigger_mm || Pass_Trigger_em || Pass_Trigger_ee)) return;

  if(k_channel.Contains("DoubleMuon")) if(!Pass_Trigger_mm && (Pass_Trigger_ee||Pass_Trigger_em)) return;
  if(k_channel.Contains("DoubleEG")) if(!Pass_Trigger_ee && (Pass_Trigger_mm||Pass_Trigger_em)) return;
  if(k_channel.Contains("MuonEG")) if(!Pass_Trigger_em && (Pass_Trigger_ee||Pass_Trigger_mm)) return;
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

  GENSignalStudy( true ); 

  // ========== Get Objects (muon, electron, jet) ====================
  std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_HN_VETO", true);
  std::vector<snu::KMuon> muons;
  std::vector<snu::KMuon> muonTightColl;
  muonTightColl.clear(); muons.clear();
  int muonVetoN = muonVetoColl.size(), muonTightN = 0, muonsN = 0;
  for(unsigned int i=0;i<muonVetoN;i++){
    if(PassID(muonVetoColl.at(i), "MUON_HN_LOOSEv7_SIP3")){
      muons.push_back(muonVetoColl.at(i));
      muonsN ++;
      if(PassID(muonVetoColl.at(i), "MUON_HN_TIGHT")){
        muonTightColl.push_back(muonVetoColl.at(i));
        muonTightN++;
      }
    }
  }

  std::vector<snu::KElectron> electronVetoColl = GetElectrons(false,true,"ELECTRON_HN_VETO");
  std::vector<snu::KElectron> electrons;
  std::vector<snu::KElectron> electronTightColl;
  electronTightColl.clear(); electrons.clear();
  int electronVetoN = electronVetoColl.size(), electronTightN = 0, electronsN = 0;
  for(unsigned int i=0;i<electronVetoN;i++){
    if(PassID(electronVetoColl.at(i), "ELECTRON_HN_FAKELOOSEv7")){
      electrons.push_back(electronVetoColl.at(i));
      electronsN ++;
      if(PassID(electronVetoColl.at(i), "ELECTRON_HN_TIGHTv4")){
        electronTightColl.push_back(electronVetoColl.at(i));
        electronTightN++;
      }
    }
  }

  if(!run_fake){
    if(!(muonsN == muonTightN && electronsN == electronTightN)) return;
    if(!(muonsN == muonVetoN && electronsN == electronVetoN)) return;
  }
  if(run_fake){
    if(muonsN == muonTightN && electronsN == electronTightN) return;
    if(!(muonsN == muonVetoN && electronsN == electronVetoN)) return;
  }

  std::vector<snu::KJet> jets = GetJets("JET_HN");
  std::vector<snu::KFatJet> fatjets = GetFatJets("FATJET_HN");
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
  bool is_Schannel = (frontTjets.size() != 0 && backTjets.size() != 0);

  // ================================================================================
  // ========== Momentum Correction ===================
  CorrectMuonMomentum(muons);
  CorrectedMETRochester(muons);
  double MET = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();

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

      Pass_Pt_em = PassEMuTriggerPt(electrons, muons);

      if(muons.at(0).Charge() != electrons.at(0).Charge()) ChargeConfig_em = "OSOF";
      else ChargeConfig_em = "SSOF";

    }
  }
  if(LeptonConfig_mm == "NULL" && LeptonConfig_ee == "NULL" && LeptonConfig_em == "NULL") return;

  std::vector<KLepton> lep;
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
  if(run_cf){
    if(region == "DiEl_OS") region = "DiEl_SS";
    else if(region == "MuEl_OS") region = "MuEl_SS";
    else return;
  }

  double HT = 0.;
  for(unsigned int i=0; i<jets.size(); i++){
    HT += jets.at(i).Pt();
  }
  for(unsigned int i=0; i<fatjets.size(); i++){
//    HT += fatjets.at(i).Pt();
  }
  double LT = 0.;
  for(unsigned int i=0; i<lep.size(); i++){
    LT += lep.at(i).Pt();
  }
  double ST = LT + HT + MET;

  bool Draw_SR = true, Draw_CR = true;
  DrawHistograms(region, lep, jets, bjets, bjetsloose, fatjets, MET, LT, HT, ST, Draw_SR, Draw_CR, muons, electrons, is_Schannel);

  return;
}// End of execute event loop

void HNDiLeptonAnalyzer::DrawHistograms(TString region, std::vector<KLepton> lep, std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets, std::vector<snu::KJet> bjetsloose, std::vector<snu::KFatJet> fatjets, double MET, double LT, double HT, double ST, bool Draw_SR, bool Draw_CR, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, bool is_Schannel){

  int cutN_SR = 0, cutN_CR = 0;
  if(Draw_SR){
    if(region == "DiMu_SS") cutN_SR = 0;
    else if(region == "MuEl_SS") cutN_SR = 0;
    else if(region == "DiEl_SS") cutN_SR = 0;

  }
  if(Draw_CR){
    if(region == "DiMu_SS") cutN_CR = 0;
    else if(region == "MuEl_SS") cutN_CR = 1;
    else if(region == "DiEl_SS") cutN_CR = 0;
  }

  double temp_weight = GetWeight(false, region, muons, electrons);
  double temp_weight_err = GetWeight(true, region, muons, electrons);

  std::vector<double> weight_updown;
  int weightN=1;
  if(run_fake) weightN=3;
  weight_updown.push_back(temp_weight);
  weight_updown.push_back(temp_weight+temp_weight_err);
  weight_updown.push_back(temp_weight-temp_weight_err);

  if(Draw_SR){
    for(unsigned int cut_it=0; cut_it<cutN_SR; cut_it++){
      if(GetCuts(region, GetCuts_name(region, cut_it, true), lep, jets, bjets, fatjets, MET, LT, HT, ST, true, is_Schannel)){
        TString cut_suffix = "_"+ GetCuts_name(region, cut_it, true);
        TString hist_prefix = "SR_"+region+cut_suffix;
        TString hist_suffix = "";
        for(unsigned int weight_it=0; weight_it<weightN; weight_it++){
          double this_weight = weight_updown.at(weight_it);
          if(weight_it == 1) hist_suffix ="_Up";
          if(weight_it == 2) hist_suffix ="_Down";

          TString this_lepton_hist_prefix ="";
          for(unsigned int lep_it=0; lep_it<lep.size(); lep_it++){
            if(lep_it == 0) this_lepton_hist_prefix = "SR_"+region+cut_suffix+"_LeadingLepton_";
            if(lep_it == 1) this_lepton_hist_prefix = "SR_"+region+cut_suffix+"_SubLeadingLepton_";
            if(lep_it == 2) this_lepton_hist_prefix = "SR_"+region+cut_suffix+"_TrailingLepton_";
            KLepton this_lepton;
            this_lepton =lep.at(lep_it);
            FillLeptonHist(this_lepton_hist_prefix, hist_suffix, this_lepton, this_weight);
          }
          FillHist(hist_prefix+"_MET"+hist_suffix, MET, this_weight, 0., 1000., 1000);
          FillHist(hist_prefix+"_METsqdivST"+hist_suffix, MET*MET/ST, this_weight, 0., 200., 200);
          FillHist(hist_prefix+"_DiLepton_Mass"+hist_suffix, (lep.at(0)+lep.at(1)).M(), this_weight, 0., 1500., 1500);
          FillHist(hist_prefix+"_DiLepton_DeltaR"+hist_suffix, lep.at(0).DeltaR(lep.at(1)), this_weight, 0., 5., 50);

          if(jets.size() > 1){

            FillHist(hist_prefix+"_LeadingJetL_Pt"+hist_suffix, jets.at(jet_lowindex[0]).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_SubLeadingJetL_Pt"+hist_suffix, jets.at(jet_lowindex[1]).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_DiJetL_Mass"+hist_suffix, (jets.at(jet_lowindex[0])+jets.at(jet_lowindex[1])).M(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_DiLeptonLeadingJetL_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_lowindex[0])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonSubLeadingJetL_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_lowindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_LeadingLeptonDiJetL_Mass"+hist_suffix, (lep.at(0)+jets.at(jet_lowindex[0])+jets.at(jet_lowindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_SubLeadingLeptonDiJetL_Mass"+hist_suffix, (lep.at(1)+jets.at(jet_lowindex[0])+jets.at(jet_lowindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonDiJetL_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_lowindex[0])+jets.at(jet_lowindex[1])).M(), this_weight, 0., 2500., 2500);
            FillHist(hist_prefix+"_DiJetL_DeltaR"+hist_suffix, jets.at(jet_lowindex[0]).DeltaR(jets.at(jet_lowindex[1])), this_weight, 0., 5., 50);

            FillHist(hist_prefix+"_LeadingJetH_Pt"+hist_suffix, jets.at(jet_highindex[0]).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_SubLeadingJetH_Pt"+hist_suffix, jets.at(jet_highindex[1]).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_DiJetH_Mass"+hist_suffix, (jets.at(jet_highindex[0])+jets.at(jet_highindex[1])).M(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_DiLeptonLeadingJetH_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_highindex[0])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonSubLeadingJetH_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_highindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_LeadingLeptonDiJetH_Mass"+hist_suffix, (lep.at(0)+jets.at(jet_highindex[0])+jets.at(jet_highindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_SubLeadingLeptonDiJetH_Mass"+hist_suffix, (lep.at(1)+jets.at(jet_highindex[0])+jets.at(jet_highindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonDiJetH_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_highindex[0])+jets.at(jet_highindex[1])).M(), this_weight, 0., 2500., 2500);
            FillHist(hist_prefix+"_DiJetH_DeltaR"+hist_suffix, jets.at(jet_highindex[0]).DeltaR(jets.at(jet_highindex[1])), this_weight, 0., 5., 50);

          }
          if(fatjets.size() > 0){
/*            FillHist(hist_prefix+"_FatJet_Pt"+hist_suffix, fatjets.at(0).Pt(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_FatJet_Mass"+hist_suffix, fatjets.at(0).M(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_LeadingLeptonFatJet_Mass"+hist_suffix, (lep.at(0)+fatjets.at(0)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_SubLeadingLeptonFatJet_Mass"+hist_suffix, (lep.at(1)+fatjets.at(0)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonFatJet_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+fatjets.at(0)).M(), this_weight, 0., 2500., 2500);
            FillHist(hist_prefix+"_FatJet_PrunedMass"+hist_suffix, fatjets.at(0).PrunedMass(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_FatJet_SoftDropMass"+hist_suffix, fatjets.at(0).SoftDropMass(), this_weight, 0., 1500., 1500);*/
          }
          FillHist(hist_prefix+"_N_of_Leptons"+hist_suffix, lep.size(), this_weight, 0., 5., 5);
          FillHist(hist_prefix+"_N_of_Jets"+hist_suffix, jets.size(), this_weight, 0., 10., 10);
          FillHist(hist_prefix+"_N_of_FatJets"+hist_suffix, fatjets.size(), this_weight, 0., 5., 5);
          FillHist(hist_prefix+"_N_of_BJets"+hist_suffix, bjets.size(), this_weight, 0., 5., 5);
          FillHist(hist_prefix+"_N_of_LooseBJets"+hist_suffix, bjetsloose.size(), this_weight, 0., 5., 5);
          FillHist(hist_prefix+"_N_of_Events"+hist_suffix, 0., this_weight, 0., 1., 1);
          FillHist(hist_prefix+"_Charge_Asymmetry"+hist_suffix, lep.at(1).Charge(), this_weight, -2., 3., 5);
        }
      }
    }
  }

  if(Draw_CR){
    for(unsigned int cut_it=0; cut_it<cutN_CR; cut_it++){
      if(GetCuts(region, GetCuts_name(region, cut_it, false), lep, jets, bjets, fatjets, MET, LT, HT, ST, false, is_Schannel)){
        TString cut_suffix = "_"+ GetCuts_name(region, cut_it, false);
        TString hist_prefix = "CR_"+region+cut_suffix;
        TString hist_suffix = "";
        for(unsigned int weight_it=0; weight_it<weightN; weight_it++){
          double this_weight = weight_updown.at(weight_it);
          if(weight_it == 1) hist_suffix ="_Up";
          if(weight_it == 2) hist_suffix ="_Down";
          TString this_lepton_hist_prefix ="";
          for(unsigned int lep_it=0; lep_it<lep.size(); lep_it++){
            if(lep_it == 0) this_lepton_hist_prefix = "CR_"+region+cut_suffix+"_LeadingLepton_";
            if(lep_it == 1) this_lepton_hist_prefix = "CR_"+region+cut_suffix+"_SubLeadingLepton_";
            if(lep_it == 2) this_lepton_hist_prefix = "CR_"+region+cut_suffix+"_TrailingLepton_";
            KLepton this_lepton;
            this_lepton =lep.at(lep_it);
            FillLeptonHist(this_lepton_hist_prefix, hist_suffix, this_lepton, this_weight);
          }
cout<<cut_suffix<<endl;
          FillHist(hist_prefix+"_MET"+hist_suffix, MET, this_weight, 0., 1000., 1000);
          FillHist(hist_prefix+"_METsqdivST"+hist_suffix, MET*MET/ST, this_weight, 0., 200., 200);
          FillHist(hist_prefix+"_DiLepton_Mass"+hist_suffix, (lep.at(0)+lep.at(1)).M(), this_weight, 0., 1500., 1500);
          FillHist(hist_prefix+"_DiLepton_DeltaR"+hist_suffix, lep.at(0).DeltaR(lep.at(1)), this_weight, 0., 5., 50);
          if(jets.size() > 1){
            FillHist(hist_prefix+"_LeadingJetL_Pt"+hist_suffix, jets.at(jet_lowindex[0]).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_SubLeadingJetL_Pt"+hist_suffix, jets.at(jet_lowindex[1]).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_DiJetL_Mass"+hist_suffix, (jets.at(jet_lowindex[0])+jets.at(jet_lowindex[1])).M(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_DiLeptonLeadingJetL_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_lowindex[0])).M(), this_weight, 0., 2000., 2000);

            FillHist(hist_prefix+"_DiLeptonSubLeadingJetL_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_lowindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_LeadingLeptonDiJetL_Mass"+hist_suffix, (lep.at(0)+jets.at(jet_lowindex[0])+jets.at(jet_lowindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_SubLeadingLeptonDiJetL_Mass"+hist_suffix, (lep.at(1)+jets.at(jet_lowindex[0])+jets.at(jet_lowindex[1])).M(), this_weight, 0., 2000., 2000);

            FillHist(hist_prefix+"_DiLeptonDiJetL_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_lowindex[0])+jets.at(jet_lowindex[1])).M(), this_weight, 0., 2500., 2500);
            FillHist(hist_prefix+"_DiJetL_DeltaR"+hist_suffix, jets.at(jet_lowindex[0]).DeltaR(jets.at(jet_lowindex[1])), this_weight, 0., 5., 50);

            FillHist(hist_prefix+"_LeadingJetH_Pt"+hist_suffix, jets.at(jet_highindex[0]).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_SubLeadingJetH_Pt"+hist_suffix, jets.at(jet_highindex[1]).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_DiJetH_Mass"+hist_suffix, (jets.at(jet_highindex[0])+jets.at(jet_highindex[1])).M(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_DiLeptonLeadingJetH_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_highindex[0])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonSubLeadingJetH_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_highindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_LeadingLeptonDiJetH_Mass"+hist_suffix, (lep.at(0)+jets.at(jet_highindex[0])+jets.at(jet_highindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_SubLeadingLeptonDiJetH_Mass"+hist_suffix, (lep.at(1)+jets.at(jet_highindex[0])+jets.at(jet_highindex[1])).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonDiJetH_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(jet_highindex[0])+jets.at(jet_highindex[1])).M(), this_weight, 0., 2500., 2500);
            FillHist(hist_prefix+"_DiJetH_DeltaR"+hist_suffix, jets.at(jet_highindex[0]).DeltaR(jets.at(jet_highindex[1])), this_weight, 0., 5., 50);

          }  
          FillHist(hist_prefix+"_N_of_Leptons"+hist_suffix, lep.size(), this_weight, 0., 5., 5);
          FillHist(hist_prefix+"_N_of_Jets"+hist_suffix, jets.size(), this_weight, 0., 10., 10);
          FillHist(hist_prefix+"_N_of_FatJets"+hist_suffix, fatjets.size(), this_weight, 0., 5., 5);
          FillHist(hist_prefix+"_N_of_BJets"+hist_suffix, bjets.size(), this_weight, 0., 5., 5);
          FillHist(hist_prefix+"_N_of_LooseBJets"+hist_suffix, bjetsloose.size(), this_weight, 0., 5., 5);
          FillHist(hist_prefix+"_N_of_Events"+hist_suffix, 0., this_weight, 0., 1., 1);
        }
      }
    }
  }


  return;

}
 
bool HNDiLeptonAnalyzer::GetCuts(TString region, TString cut, std::vector<KLepton> lep, std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets, std::vector<snu::KFatJet> fatjets, double MET, double LT, double HT, double ST, bool Is_SR, bool is_Schannel){

  if(Is_SR){
    if(cut == "Preselection"){

      if((lep.at(0)+lep.at(1)).M() < 10) return false;
      if(jets.size() < 2) return false;
      if(region == "DiEl_SS") if(fabs((lep.at(0)+lep.at(1)).M() - 91.1876) < 10) return false;


      snu::KParticle W_candidate, W_selection;
      W_selection = lep.at(0)+lep.at(1)+jets.at(0)+jets.at(1);
      for(unsigned int i=0; i<jets.size(); i++){
        for(unsigned int j=i+1; j<jets.size(); j++){
          W_candidate = lep.at(0)+lep.at(1)+jets.at(i)+jets.at(j);
          if(fabs(W_candidate.M() - 80.4) < fabs(W_selection.M() - 80.4)){
            W_selection = W_candidate;
            jet_lowindex[0] = i;
            jet_lowindex[1] = j;
          }
        }
      }

      W_selection = jets.at(0)+jets.at(1);
      for(unsigned int i=0; i<jets.size(); i++){
        for(unsigned int j=i+1; j<jets.size(); j++){
          W_candidate = jets.at(i)+jets.at(j);
          if(fabs(W_candidate.M() - 80.4) < fabs(W_selection.M() - 80.4)){
            W_selection = W_candidate;
            jet_highindex[0] = i;
            jet_highindex[1] = j;
          }
        }
      }

      Pass_Preselection = true;


      return true;
    }

    if(cut == "LowMass"){
      return true;
    }

    if(cut == "HighMass"){

      return true;
    }

  }

  if(!Is_SR){
    if(cut == "LowMass"){
      if(!Pass_Preselection) return false;
      if((lep.at(0)+lep.at(1)+jets.at(jet_lowindex[0])+jets.at(jet_lowindex[1])).M() > 300) return false;
      if(bjets.size() == 0 && MET < 100.) return false; 

      return true;
    }
    if(cut == "HighMass"){
      if(!Pass_Preselection) return false;
      if((jets.at(jet_highindex[0])+jets.at(jet_highindex[1])).M() > 150) return false;
      if(bjets.size() == 0 && (MET*MET/ST) < 20.) return false;

      return true;
    }
    if(cut == "1Jet0bJet"){
      if((lep.at(0)+lep.at(1)).M() < 100) return false;
      if(jets.size() != 1) return false;
      if(bjets.size() != 0) return false;
      
      return true;
    }
    if(cut == "1Jet0bJetDiLepMass10"){
      if((lep.at(0)+lep.at(1)).M() < 10) return false;
      if(jets.size() != 1) return false;
      if(bjets.size() != 0) return false;

      return true;
    }

    if(cut == "SSEMu"){
      if((lep.at(0)+lep.at(1)).M() < 10) return false;
      if((lep.at(0).Charge() != lep.at(1).Charge())) return false;
 
      return true;
    }

  }

  return false;
}

TString HNDiLeptonAnalyzer::GetCuts_name(TString region, int cut, bool Is_SR){

  if(Is_SR){
    if(cut == 0) return "Preselection";
    if(cut == 1) return "LowMass";
    if(cut == 2) return "HighMass";

    return "NULL";
  }
  if(!Is_SR){
    if(cut == 0) return "SSEMu";
/*    if(cut == 0) return "LowMass";
    if(cut == 1) return "HighMass";
    if(cut == 2) return "1Jet0bJet";
    if(cut == 3) return "1Jet0bJetDiLepMass10";*/

    return "NULL";
  }

  return "NULL";
}

void HNDiLeptonAnalyzer::FillLeptonHist(TString hist_prefix, TString hist_suffix, KLepton this_lep, double this_weight){

  FillHist(hist_prefix+"Pt"+hist_suffix, this_lep.Pt(), this_weight, 0., 1000., 1000);
  FillHist(hist_prefix+"Eta"+hist_suffix, this_lep.Eta(), this_weight, -3., 3., 60);
  FillHist(hist_prefix+"RelIso"+hist_suffix, this_lep.RelIso(), this_weight, 0., 1., 100);
  FillHist(hist_prefix+"dXY"+hist_suffix, fabs(this_lep.dXY()), this_weight, 0., 0.01, 100);
  FillHist(hist_prefix+"dZ"+hist_suffix, fabs(this_lep.dZ()), this_weight, 0., 0.05, 100);
  FillHist(hist_prefix+"Flavour"+hist_suffix, this_lep.LeptonFlavour(), this_weight, 0., 3., 3);

}

double HNDiLeptonAnalyzer::GetWeight(bool geterr, TString region, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons){

  if(!run_fake && !run_cf){
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
    if(region == "MuEl_OS" || region == "MuEl_SS"){
      bool Pass_NoDZ = (PassTriggerOR(triggerlist_emNoDZ1) || PassTriggerOR(triggerlist_emNoDZ2));
      bool Pass_DZ = (PassTriggerOR(triggerlist_emDZ1) || PassTriggerOR(triggerlist_emDZ2));
      if(Pass_NoDZ && Pass_DZ) trigger_ps_weight = (-2041.112) + (-7540.488);

      if(Pass_NoDZ) trigger_ps_weight += WeightByTrigger(triggerlist_emNoDZ1, TargetLumi);
      if(Pass_DZ) trigger_ps_weight += WeightByTrigger(triggerlist_emDZ1, TargetLumi);
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
  if(run_fake && !run_cf){
    weight = get_eventweight(false, muons, electrons);
    weight_err = get_eventweight(true, muons, electrons);

//    weight = m_datadriven_bkg->Get_DataDrivenWeight(false,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSEv7", "mva");
//    weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSEv7", "mva");
  }

  if(!run_fake && run_cf){
    weight = GetCFweight(0, electrons, true, "ELECTRON_HN_TIGHTv4");
    weight_err = 0.;
  }
  if(run_fake && run_cf){

//    weight = m_datadriven_bkg->Get_DataDrivenWeight(false,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSEv7", "mva");
//    weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true,  muons, "MUON_HN_TIGHT", muons.size(), electrons, "ELECTRON_HN_TIGHTv4", electrons.size(), "ELECTRON_HN_FAKELOOSEv7", "mva");
//    weight *= GetCFweight(0, electrons, true, "ELECTRON_HN_TIGHTv4");
//    weight = weight*(-1.);    
  }

  if(geterr) return weight_err;
  if(!geterr) return weight;

}

bool HNDiLeptonAnalyzer::PassEMuTriggerPt(std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons){
  
  bool pass =false;
  snu::KParticle el,mu;
  el = electrons.at(0);
  mu = muons.at(0);

  if(PassTriggerOR(triggerlist_emNoDZ1)){
    if((mu.Pt() >10 && el.Pt() >25)) return true;
  }
  if(PassTriggerOR(triggerlist_emNoDZ2)){
    if((mu.Pt() >25 && el.Pt() >10)) return true;
  }
  if(PassTriggerOR(triggerlist_emDZ1)){
    if((mu.Pt() >10 && el.Pt() >25)) return true;
  }
  if(PassTriggerOR(triggerlist_emDZ2)){
    if((mu.Pt() >25 && el.Pt() >10)) return true;
  }
  
  return pass;
}

void HNDiLeptonAnalyzer::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNDiLeptonAnalyzer::BeginCycle() throw( LQError ){
  
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

HNDiLeptonAnalyzer::~HNDiLeptonAnalyzer() {
  
  Message("In HNDiLeptonAnalyzer Destructor" , INFO);
  
}


void HNDiLeptonAnalyzer::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNDiLeptonAnalyzer::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNDiLeptonAnalyzerCore::MakeHistograms() to make new hists for your analysis
   **/
 
}


void HNDiLeptonAnalyzer::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


void HNDiLeptonAnalyzer::GENSignalStudy( bool Is_Signal ){

  //if( !(k_sample_name.Contains("Tchannel"))) return;
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

//  TruthPrintOut();
//  return;

  int electron_ID = 11, muon_ID = 13, W_ID = 24, HN_ID = 9900012;
  int parton_ID = 1;

  //TruthPrintOut();

  int max = truthColl.size();
  vector<int> HN_index, onshell_W_index, lep1_index, lep2_index, quark1_index, quark2_index, vbf_quark_index;
  HN_index.clear(); onshell_W_index.clear(); lep1_index.clear(); lep2_index.clear(); quark1_index.clear(); quark2_index.clear(); vbf_quark_index.clear();

  // Look for Heavy Neutrino using PdgId
  for( int i = 2 ; i < max ; i++ ){
    if( abs(truthColl.at(i).PdgId()) == HN_ID ){
      HN_index.push_back(i);
      GENFindDecayIndex( truthColl, i, HN_index );
      break;
    }
  }
  if( HN_index.size() == 0 ){
    FillHist("GEN_HN_not_found", 0., 1., 0., 1., 1);
    return;
  }
  int HN_mother_index = truthColl.at( HN_index.at(0) ).IndexMother(); // Save HN mother index for lep1

  // Look for muon1 using PdgId and mother index (sister : HN)
  for( int i = 2 ; i < max ; i++ ){
    if( ((abs(truthColl.at(i).PdgId()) == muon_ID) || (abs(truthColl.at(i).PdgId()) == electron_ID)) && truthColl.at(i).IndexMother() == HN_mother_index ){
      lep1_index.push_back(i);
      GENFindDecayIndex( truthColl, i, lep1_index);
      break;
    }
  }
  if( lep1_index.size () == 0 ){
    FillHist("GEN_lep1_not_found", 0., 1., 0., 1., 1);
    return;
  }

  // Look for muon2 using PdgId and mother index (mother : HN)
  int HN_decay_index = 0;
  for( int index = 0 ; index < HN_index.size() ; index++ ){ // One of the HN index decays into muon
    for( int i = 2 ; i < max ; i++ ){
      if( ((abs(truthColl.at(i).PdgId()) == muon_ID) || (abs(truthColl.at(i).PdgId()) == electron_ID)) && truthColl.at(i).IndexMother() == HN_index.at(index) ){
        lep2_index.push_back(i);
        GENFindDecayIndex( truthColl, i, lep2_index);
        HN_decay_index = index;
        break;
      }
    }
  }
  if( lep2_index.size () == 0 ){
    FillHist("GEN_lep2_not_found", 0., 1., 0., 1., 1);
    return;
  }

  for( int i=2; i < max ; i++){
    if( abs(truthColl.at(i).PdgId()) < 6 && abs(truthColl.at(truthColl.at(i).IndexMother()).PdgId()) == HN_ID ){
      quark1_index.push_back(i);
      GENFindDecayIndex( truthColl, i, quark1_index);
      break;
    }
  }
  if( quark1_index.size () == 0 ){
    FillHist("GEN_quark1_not_found", 0., 1., 0., 1., 1);
    return;
  }

  for( int i=(quark1_index.at(0)+1); i < max ; i++){
    if( abs(truthColl.at(i).PdgId()) < 6 && abs(truthColl.at(truthColl.at(i).IndexMother()).PdgId()) == HN_ID ){
      quark2_index.push_back(i);
      GENFindDecayIndex( truthColl, i, quark2_index);
      break;
    }
  }
  if( quark2_index.size () == 0 ){
    FillHist("GEN_quark2_not_found", 0., 1., 0., 1., 1);
    return;
  }

  for( int i=0; i < max; i++){
    if( (abs(truthColl.at(i).PdgId()) < 6) && (abs(truthColl.at(truthColl.at(i).IndexMother()).PdgId()) != W_ID) && (truthColl.at(i).GenStatus() == 23)){
      vbf_quark_index.push_back(i);
      GENFindDecayIndex( truthColl, i, vbf_quark_index );
    }
  }
  snu::KParticle GEN_lep[2], GEN_quark[2], GEN_onshellW, GEN_HN, GEN_vbfquark;
  bool Is_Tchannel = false;
  if( vbf_quark_index.size() != 0 ) Is_Tchannel = true;
  GEN_lep[0] = truthColl.at( lep1_index.back() );
  GEN_lep[1] = truthColl.at( lep2_index.back() );

  GEN_quark[0] = truthColl.at( quark1_index.back() );
  GEN_quark[1] = truthColl.at( quark2_index.back() );
  if(GEN_quark[0].Pt() < GEN_quark[1].Pt()){
    GEN_quark[1] = truthColl.at( quark1_index.back() );
    GEN_quark[0] = truthColl.at( quark2_index.back() );
  }
  if( Is_Tchannel ) GEN_vbfquark = truthColl.at( vbf_quark_index.back() );

  GEN_onshellW = GEN_quark[0] + GEN_quark[1];
  GEN_HN = GEN_onshellW + GEN_lep[1];
  FillHist("GEN_PrimaryLepton_Pt", GEN_lep[0].Pt(), 1., 0., 1000., 1000);
  FillHist("GEN_SecondaryLepton_Pt", GEN_lep[1].Pt(), 1., 0., 1000., 1000);
  FillHist("GEN_LeadingQuark_Pt", GEN_quark[0].Pt(), 1., 0., 1000., 1000);
  FillHist("GEN_SubLeadingQuark_Pt", GEN_quark[1].Pt(), 1., 0., 1000., 1000);
  FillHist("GEN_OnShellW_Mass", GEN_onshellW.M(), 1., 0., 500., 500);
  FillHist("GEN_HeavyNeutrino_Mass", GEN_HN.M(), 1., 0., 2500., 2500);
  FillHist("GEN_DiJet_DeltaR", GEN_quark[0].DeltaR(GEN_quark[1]), 1., 0., 5., 50);
  FillHist("GEN_PrimaryLeptonAndonshellW_DeltaR", GEN_onshellW.DeltaR(GEN_lep[0]), 1., 0., 5., 50);
  FillHist("GEN_SecondaryLeptonAndonshellW_DeltaR", GEN_onshellW.DeltaR(GEN_lep[1]), 1., 0., 5., 50);
  FillHist("GEN_PrimaryLeptonAndonshellW_DeltaPhi", GEN_onshellW.DeltaPhi(GEN_lep[0]), 1., 0., 5., 50);
  FillHist("GEN_SecondaryLeptonAndonshellW_DeltaPhi", GEN_onshellW.DeltaPhi(GEN_lep[1]), 1., 0., 5., 50);
  FillHist("GEN_PrimaryLeptonAndonshellW_CosDeltaPhi", TMath::Cos(GEN_onshellW.DeltaPhi(GEN_lep[0])), 1., -1.5, 1.5, 30);
  FillHist("GEN_SecondaryLeptonAndonshellW_CosDeltaPhi", TMath::Cos(GEN_onshellW.DeltaPhi(GEN_lep[1])), 1., -1.5, 1.5, 30);
  
  if( Is_Tchannel ){
    FillHist("VBF_GEN_ForwardJet_Pt", GEN_vbfquark.Pt(), 1., 0., 500., 500);
    FillHist("VBF_GEN_ForwardJet_Eta", GEN_vbfquark.Eta(), 1., -5., 5., 100);
    FillHist("VBF_GEN_ForwardJetAndHN_DeltaEta", fabs(GEN_vbfquark.Eta()-GEN_HN.Eta()), 1., 0., 20., 40);
    FillHist("VBF_GEN_ForwardJetAndonshellW_DeltaEta", fabs(GEN_vbfquark.Eta()-GEN_onshellW.Eta()), 1., 0., 20., 40);
    if(fabs(GEN_vbfquark.Eta()) > 2.5){
      FillHist("VBF_GEN_2p5ForwardJet_Pt", GEN_vbfquark.Pt(), 1., 0., 500., 500);
      FillHist("VBF_GEN_2p5ForwardJet_Eta", GEN_vbfquark.Eta(), 1., -5., 5., 100);
      FillHist("VBF_GEN_2p5ForwardJetAndHN_DeltaEta", fabs(GEN_vbfquark.Eta()-GEN_HN.Eta()), 1., 0., 20., 40);
      FillHist("VBF_GEN_2p5ForwardJetAndonshellW_DeltaEta", fabs(GEN_vbfquark.Eta()-GEN_onshellW.Eta()), 1., 0., 20., 40);
    }
  }
  return;  

}


void HNDiLeptonAnalyzer::GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index ){

  for( int i = it+1 ; i < truthColl.size(); i ++ ){
    if( truthColl.at(i).IndexMother() == it && truthColl.at(i).PdgId() == truthColl.at(it).PdgId() ){
      index.push_back(i);
    }
  }
  return;

}

double HNDiLeptonAnalyzer::CorrPt(KLepton lep, double T_iso){

  double ptcorr = lep.Pt()*(1+max(0.,(lep.RelIso()-T_iso)));
  return ptcorr;
}

double HNDiLeptonAnalyzer::CorrPt(snu::KMuon lep, double T_iso){

  double ptcorr = lep.Pt()*(1+max(0.,(lep.RelIso04()-T_iso)));
  return ptcorr;
}

double HNDiLeptonAnalyzer::CorrPt(snu::KElectron lep, double T_iso){

  double ptcorr = lep.Pt()*(1+max(0.,(lep.PFRelIso(0.3)-T_iso)));
  return ptcorr;
}



double HNDiLeptonAnalyzer::GetMuonFR(bool geterr, float pt,  float eta){

  if(pt < 10.) pt = 11.;
  if(pt >= 60.) pt = 59.;
  if(fabs(eta) >= 2.4) eta = 2.3;

  //cout << "[HNDiLeptonAnalyzer::GetMuonFR] pt = " << pt << endl;
  //cout << "[HNDiLeptonAnalyzer::GetMuonFR] eta = " << eta << endl;
  //cout << "[HNDiLeptonAnalyzer::GetMuonFR] MuFR_key = " << MuFR_key << endl;
  //cout << "[HNDiLeptonAnalyzer::GetMuonFR] NearBjet = " << NearBjet << endl;

  TH2D *THISFRHIST = hist_Muon_FR;

  int binx = THISFRHIST->FindBin(pt, abs(eta));
  //cout << "[HNDiLeptonAnalyzer::GetMuonFR] => FR = " << THISFRHIST->GetBinContent(binx) << endl;
  

  if(geterr) return THISFRHIST->GetBinError(binx);
  else return THISFRHIST->GetBinContent(binx);

}

double HNDiLeptonAnalyzer::GetMuonPR(bool geterr, float pt,  float eta){
/*
  if(pt < 10.) pt = 11.;
  if(pt >= 60.) pt = 59.;
  if(fabs(eta) >= 2.5) eta = 2.4;
  int binx = hist_Muon_PR->FindBin(pt, abs(eta));
  if(geterr) return hist_Muon_PR->GetBinError(binx);
  else return hist_Muon_PR->GetBinContent(binx);
*/
  return 1.;
}


double HNDiLeptonAnalyzer::get_eventweight(bool geterr, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons){

  unsigned int n_leptons = muons.size() + electrons.size();
  //cout << "[HNDiLeptonAnalyzer::get_eventweight] muons.size() = " << muons.size() << ", electrons.size() = " << electrons.size() << endl;

  vector<float> lep_pt, lep_eta;
  vector<bool> ismuon;
  vector<bool> isT;
  for(unsigned int i=0; i<muons.size(); i++){
    lep_pt.push_back( CorrPt(muons.at(i), 0.07) );
    lep_eta.push_back(muons.at(i).Eta());
    ismuon.push_back(true);
    if(PassID(muons.at(i), "MUON_HN_TIGHT")) isT.push_back(true);
    else isT.push_back(false);
  }
  for(unsigned int i=0; i<electrons.size(); i++){
    lep_pt.push_back( CorrPt(electrons.at(i), 0.08) );
    lep_eta.push_back(electrons.at(i).Eta());
    ismuon.push_back(false);
    if(PassID(electrons.at(i), "ELECTRON_HN_TIGHTv4")) isT.push_back(true);
    else isT.push_back(false);
  }

  vector<float> fr, pr, fr_err, pr_err;

  for(unsigned int i=0; i<n_leptons; i++){
    //==== Muon
    if(ismuon.at(i)){
      fr.push_back( GetMuonFR(false, lep_pt.at(i), lep_eta.at(i)) );
      pr.push_back( GetMuonPR(false, lep_pt.at(i), lep_eta.at(i)) );
      fr_err.push_back( GetMuonFR(true, lep_pt.at(i), lep_eta.at(i)) );
      pr_err.push_back( GetMuonPR(true, lep_pt.at(i), lep_eta.at(i)) );
    }
    else{
      fr.push_back(0);// GetElectronFR(0, lep_pt.at(i), lep_eta.at(i)) );
      pr.push_back(0);// GetElectronPR(0, lep_pt.at(i), lep_eta.at(i)) );
      fr_err.push_back(0);// GetElectronFR(1, lep_pt.at(i), lep_eta.at(i)) );
      pr_err.push_back(0);// GetElectronPR(1, lep_pt.at(i), lep_eta.at(i)) );
    }
  }

  //==== let a == f/(1-f)

  vector<float> a, fr_onlyLoose;

  for(unsigned int i=0; i<n_leptons; i++) a.push_back( fr.at(i)/(1.-fr.at(i)) );
  for(unsigned int i=0; i<n_leptons; i++){
    if(!isT.at(i)){
      //cout << "[HNDiLeptonAnalyzer::get_eventweight] "<<i<<" th lepton is Loose" << endl;
      fr_onlyLoose.push_back( a.at(i) );
    }
  }

  //==== Initialise weight
  float this_weight=-1.;

  for(unsigned int i=0; i<fr_onlyLoose.size(); i++){
    this_weight *= -fr_onlyLoose.at(i);
  }
  //cout << "[HNDiLeptonAnalyzer::get_eventweight] this_weight = " << this_weight << endl;

  //==== d(a)/a = d(f)/f(1-f)
  //==== so, if w = a1*a2,
  //==== d(w)/w = d(a1)/a1 + d(a2)/a2

  vector<float> da_over_a;
  for(unsigned int i=0; i<n_leptons; i++) da_over_a.push_back( fr_err.at(i) / ( fr.at(i)*(1.-fr.at(i)) ) );
  float this_weight_err = 0.;
  for(unsigned int i=0; i<n_leptons; i++){
    if(!isT.at(i)) this_weight_err += da_over_a.at(i)*da_over_a.at(i);
  }

  this_weight_err = sqrt(this_weight_err);
  this_weight_err = this_weight_err*fabs(this_weight);

  if(!geterr) return this_weight;
  if(geterr) return this_weight_err;


}
