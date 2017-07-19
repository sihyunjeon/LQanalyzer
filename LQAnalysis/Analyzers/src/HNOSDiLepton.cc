// $Id: HNOSDiLepton.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNOSDiLepton Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNOSDiLepton.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNOSDiLepton);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNOSDiLepton::HNOSDiLepton() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNOSDiLepton");
  
  Message("In HNOSDiLepton constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNOSDiLepton::InitialiseAnalysis() throw( LQError ) {
  
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

}


void HNOSDiLepton::ExecuteEvents()throw( LQError ){

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

  GENSignalStudy( (k_sample_name.Contains("HN")) );

  triggerlist_mm.clear(); triggerlist_emBG1.clear(); triggerlist_emBG2.clear(); triggerlist_emH1.clear(); triggerlist_emH2.clear(); triggerlist_ee.clear();
  // ========== Trigger cut ====================
  triggerlist_mm.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_mm.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  triggerlist_emBG1.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_emBG2.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");

  triggerlist_emH1.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_emH2.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

//  triggerlist_me1.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
//  triggerlist_me2.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_ee.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  bool Pass_Trigger_mm = PassTriggerOR(triggerlist_mm);
  bool Pass_Trigger_em = (PassTriggerOR(triggerlist_emBG1) || PassTriggerOR(triggerlist_emBG2) || PassTriggerOR(triggerlist_emH1) || PassTriggerOR(triggerlist_emH2));
  bool Pass_Trigger_ee = PassTriggerOR(triggerlist_ee);

  if(!(Pass_Trigger_mm && Pass_Trigger_em && Pass_Trigger_ee)) return;

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

      if((PassTriggerOR(triggerlist_emBG1) && PassTriggerOR(triggerlist_emBG2)) || (PassTriggerOR(triggerlist_emH1) && PassTriggerOR(triggerlist_emH2)))
        Pass_Pt_em = (muons.at(0).Pt() > 10 && electrons.at(0).Pt() > 10);
      else if((PassTriggerOR(triggerlist_emBG1) && PassTriggerOR(triggerlist_emH2)) || (PassTriggerOR(triggerlist_emH1) && PassTriggerOR(triggerlist_emBG2)))
        Pass_Pt_em = (muons.at(0).Pt() > 10 && electrons.at(0).Pt() > 10);

      else if((PassTriggerOR(triggerlist_emBG1) || PassTriggerOR(triggerlist_emH1)) && !(PassTriggerOR(triggerlist_emBG2) || PassTriggerOR(triggerlist_emH2)))
        Pass_Pt_em = (muons.at(0).Pt() > 10 && electrons.at(0).Pt() > 25);
      else if(!(PassTriggerOR(triggerlist_emBG1) || PassTriggerOR(triggerlist_emH1)) && (PassTriggerOR(triggerlist_emBG2) || PassTriggerOR(triggerlist_emH2)))
        Pass_Pt_em = (muons.at(0).Pt() > 25 && electrons.at(0).Pt() > 10);

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

void HNOSDiLepton::DrawHistograms(TString region, std::vector<KLepton> lep, std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets, std::vector<snu::KJet> bjetsloose, std::vector<snu::KFatJet> fatjets, double MET, double LT, double HT, double ST, bool Draw_SR, bool Draw_CR, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons, bool is_Tchannel){

  int cutN_SR = 0, cutN_CR = 0;
  if(Draw_SR){
    if(region == "DiMu_SS") cutN_SR = 0;
    if(region == "DiMu_OS") cutN_SR = 3;
    if(region == "DiEl_SS") cutN_SR = 0;
    if(region == "DiEl_OS") cutN_SR = 3;
    if(region == "MuEl_SS") cutN_SR = 3;
    if(region == "MuEl_OS") cutN_SR = 3;
  }
  if(Draw_CR){
    if(region == "DiMu_SS") cutN_CR = 0;
    if(region == "DiMu_OS") cutN_CR = 2;
    if(region == "DiEl_SS") cutN_CR = 0;
    if(region == "DiEl_OS") cutN_CR = 0;
    if(region == "MuEl_SS") cutN_SR = 1;
    if(region == "MuEl_OS") cutN_CR = 2;
  }
  double temp_weight = GetWeight(false, region, muons, electrons);
  double temp_weight_err = GetWeight(true, region, muons, electrons);

  if(k_running_chargeflip){
    if(region == "DiEl_OS") region = "DiEl_SS";
    if(region == "MuEl_OS") region = "MuEl_SS";
  }

  std::vector<double> weight_updown;

  int weightN=1;
  if(k_running_nonprompt) weightN=3;
  weight_updown.push_back(temp_weight);
  weight_updown.push_back(temp_weight+temp_weight_err);
  weight_updown.push_back(temp_weight-temp_weight_err);

  if(Draw_SR){
    for(unsigned int cut_it=0; cut_it<cutN_SR; cut_it++){
      if(GetCuts(region, GetCuts_name(region, cut_it, true), lep, jets, bjets, fatjets, MET, LT, HT, ST, true, is_Tchannel)){
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
            FillHist(hist_prefix+"_LeadingJet_Pt"+hist_suffix, jets.at(0).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_SubLeadingJet_Pt"+hist_suffix, jets.at(1).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_DiJet_Mass"+hist_suffix, (jets.at(0)+jets.at(1)).M(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_DiLeptonLeadingJet_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(0)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonSubLeadingJet_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(1)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_LeadingLeptonDiJet_Mass"+hist_suffix, (lep.at(0)+jets.at(0)+jets.at(1)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_SubLeadingLeptonDiJet_Mass"+hist_suffix, (lep.at(1)+jets.at(0)+jets.at(1)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonDiJet_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(0)+jets.at(1)).M(), this_weight, 0., 2500., 2500);
            FillHist(hist_prefix+"_DiJet_DeltaR"+hist_suffix, jets.at(0).DeltaR(jets.at(1)), this_weight, 0., 5., 50);
          }
          if(fatjets.size() > 0){
            FillHist(hist_prefix+"_FatJet_Pt"+hist_suffix, fatjets.at(0).Pt(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_FatJet_Mass"+hist_suffix, fatjets.at(0).M(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_LeadingLeptonFatJet_Mass"+hist_suffix, (lep.at(0)+fatjets.at(0)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_SubLeadingLeptonFatJet_Mass"+hist_suffix, (lep.at(1)+fatjets.at(0)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonFatJet_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+fatjets.at(0)).M(), this_weight, 0., 2500., 2500);
            FillHist(hist_prefix+"_FatJet_PrunedMass"+hist_suffix, fatjets.at(0).PrunedMass(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_FatJet_SoftDropMass"+hist_suffix, fatjets.at(0).SoftDropMass(), this_weight, 0., 1500., 1500);

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

  if(Draw_CR){
    for(unsigned int cut_it=0; cut_it<cutN_CR; cut_it++){
      if(GetCuts(region, GetCuts_name(region, cut_it, false), lep, jets, bjets, fatjets, MET, LT, HT, ST, false, is_Tchannel)){
        TString cut_suffix = "_"+ GetCuts_name(region, cut_it, false);
        TString hist_prefix = "CR_"+region+cut_suffix;
        TString hist_suffix = "";

        for(unsigned int weight_it=0; weight_it<weightN; weight_it++){
          double this_weight = weight_updown.at(weight_it);
          if(weight_it == 1) hist_suffix ="_Up";
          if(weight_it == 2) hist_suffix ="_Down";
          TString lep_hist_prefix = "";
          for(unsigned int lep_it=0; lep_it<lep.size(); lep_it++){
            if(lep_it == 0) lep_hist_prefix = "CR_"+region+cut_suffix+"_LeadingLepton_";
            if(lep_it == 1) lep_hist_prefix = "CR_"+region+cut_suffix+"_SubLeadingLepton_";
            if(lep_it == 2) lep_hist_prefix = "CR_"+region+cut_suffix+"_TrailingLepton_";
            KLepton this_lepton;
            this_lepton =lep.at(lep_it);
            FillLeptonHist(lep_hist_prefix, hist_suffix, this_lepton, this_weight);
          }
          FillHist(hist_prefix+"_MET"+hist_suffix, MET, this_weight, 0., 1000., 1000);
          FillHist(hist_prefix+"_DiLepton_Mass"+hist_suffix, (lep.at(0)+lep.at(1)).M(), this_weight, 0., 1500., 1500);
          FillHist(hist_prefix+"_DiLepton_DeltaR"+hist_suffix, lep.at(0).DeltaR(lep.at(1)), this_weight, 0., 5., 50);
  
          if(jets.size() > 1){
            FillHist(hist_prefix+"_LeadingJet_Pt"+hist_suffix, jets.at(0).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_SubLeadingJet_Pt"+hist_suffix, jets.at(1).Pt(), this_weight, 0., 1000., 1000);
            FillHist(hist_prefix+"_DiJet_Mass"+hist_suffix, (jets.at(0)+jets.at(1)).M(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_DiLeptonLeadingJet_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(0)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonSubLeadingJet_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(1)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_LeadingLeptonDiJet_Mass"+hist_suffix, (lep.at(0)+jets.at(0)+jets.at(1)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_SubLeadingLeptonDiJet_Mass"+hist_suffix, (lep.at(1)+jets.at(0)+jets.at(1)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonDiJet_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+jets.at(0)+jets.at(1)).M(), this_weight, 0., 2500., 2500);
            FillHist(hist_prefix+"_DiJet_DeltaR"+hist_suffix, jets.at(0).DeltaR(jets.at(1)), this_weight, 0., 5., 50);
          }
          if(fatjets.size() > 0){
            FillHist(hist_prefix+"_FatJet_Pt"+hist_suffix, fatjets.at(0).Pt(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_FatJet_Mass"+hist_suffix, fatjets.at(0).M(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_LeadingLeptonFatJet_Mass"+hist_suffix, (lep.at(0)+fatjets.at(0)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_SubLeadingLeptonFatJet_Mass"+hist_suffix, (lep.at(1)+fatjets.at(0)).M(), this_weight, 0., 2000., 2000);
            FillHist(hist_prefix+"_DiLeptonFatJet_Mass"+hist_suffix, (lep.at(0)+lep.at(1)+fatjets.at(0)).M(), this_weight, 0., 2500., 2500);
            FillHist(hist_prefix+"_FatJet_PrunedMass"+hist_suffix, fatjets.at(0).PrunedMass(), this_weight, 0., 1500., 1500);
            FillHist(hist_prefix+"_FatJet_SoftDropMass"+hist_suffix, fatjets.at(0).SoftDropMass(), this_weight, 0., 1500., 1500);

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

  if(std::find(k_flags.begin(), k_flags.end(), "cutflow") !=k_flags.end()){
    double this_weight = weight_updown.at(0);
    if(region == "DiMu_OS"||region == "DiMu_SS"){
      FillHist("CUTFLOWCHECK_STEP3", 0., this_weight, 0., 1., 1);
      if((lep.at(0)+lep.at(1)).M()>10){
        FillHist("CUTFLOWCHECK_STEP4", 0., this_weight, 0., 1., 1);
        if(jets.size()>1){
          FillHist("CUTFLOWCHECK_STEP5", 0., this_weight, 0., 1., 1);
          if(MET<50){
            FillHist("CUTFLOWCHECK_STEP6", 0., this_weight, 0., 1., 1);
            if((jets.at(0)+jets.at(1)).M()<200){
              FillHist("CUTFLOWCHECK_STEP7", 0., this_weight, 0., 1., 1);
            }
          }
        }
      }
    }
  }
  return;

}
 
bool HNOSDiLepton::GetCuts(TString region, TString cut, std::vector<KLepton> lep, std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets, std::vector<snu::KFatJet> fatjets, double MET, double LT, double HT, double ST, bool Is_SR, bool is_Tchannel){

  if(Is_SR){
    if(cut == "Preselection"){

      if(jets.size() < 2 && fatjets.size() < 1) return false;
      if(bjets.size() != 0) return false;
      if(region == "DiMu_OS" || region == "DiEl_OS") if((lep.at(0)+lep.at(1)).M() < 110) return false;

      Pass_Preselection = true;
      return true;
    }

    if(cut == "HighMassNormalJets"){
      if(!Pass_Preselection) return false;
      if(jets.size() < 2) return false;
      snu::KParticle W_candidate, W_selection;
      W_selection = jets.at(0)+jets.at(1);
      for(unsigned int i=0; i<jets.size(); i++){
        for(unsigned int j=i+1; j<jets.size(); j++){
          W_candidate = jets.at(i)+jets.at(j);
          if(fabs(W_candidate.M() - 80.4) < fabs(W_selection.M() - 80.4)) W_selection = W_candidate;
        }
      }
      if((W_selection.M() < 50) || (W_selection.M() > 110)) return false;

      return true;
    }

    if(cut == "HighMassFatJets"){
      if(!Pass_Preselection) return false;
      if(fatjets.size() < 1) return false;

      return true;
    }
  }


  if(!Is_SR){
    if(cut == "DrellYan"){
      if((lep.at(0)+lep.at(1)).M() < 15) return false;
      if(region == "DiMu_OS" || region == "DiEl_OS") if((lep.at(0)+lep.at(1)).M() > 91.1876+10 || (lep.at(0)+lep.at(1)).M() < 91.1876-10) return false;
      if(MET > 40) return false;
      if(bjets.size() != 0) return false;
      return true;
    }
    if(cut == "TTbar"){
      if((lep.at(0)+lep.at(1)).M() < 15) return false;
      if(region == "DiMu_OS" || region == "DiEl_OS") if((lep.at(0)+lep.at(1)).M() < 91.1876+10 && (lep.at(0)+lep.at(1)).M() > 91.1876-10) return false;
      if(MET < 40) return false;
      if(bjets.size() == 0) return false;
      return true;
    }
    if(cut == "ChargeFlip"){
      if((lep.at(0)+lep.at(1)).M() < 15) return false;
      //if((lep.at(0)+lep.at(1)).M() > 91.1876+10 && (lep.at(0)+lep.at(1)).M() < 91.1876-10) return false;
      if(jets.size() > 2) return false;
      return true;
    }
  }

  return false;
}

TString HNOSDiLepton::GetCuts_name(TString region, int cut, bool Is_SR){

  if(Is_SR){
    if(region == "DiMu_OS") if(cut == 0) return "Preselection";
    if(region == "DiMu_OS") if(cut == 1) return "HighMassNormalJets";
    if(region == "DiMu_OS") if(cut == 2) return "HighMassFatJets";
    if(region == "MuEl_OS") if(cut == 0) return "Preselection";
    if(region == "MuEl_OS") if(cut == 1) return "HighMassNormalJets";
    if(region == "MuEl_OS") if(cut == 2) return "HighMassFatJets";
    if(region == "MuEl_SS") if(cut == 0) return "Preselection";
    if(region == "MuEl_SS") if(cut == 1) return "HighMassNormalJets";
    if(region == "MuEl_SS") if(cut == 2) return "HighMassFatJets"; 

    return "NULL";
  }
  if(!Is_SR){
    if(region == "DiMu_OS") if(cut == 0) return "DrellYan";
    if(region == "DiMu_OS") if(cut == 1) return "TTbar";
    if(region == "MuEl_OS") if(cut == 0) return "DrellYan";
    if(region == "MuEl_OS") if(cut == 1) return "TTbar";
    if(region == "MuEL_SS") if(cut == 0) return "ChargeFlip";

    return "NULL";
  }

}

void HNOSDiLepton::FillLeptonHist(TString hist_prefix, TString hist_suffix, KLepton this_lep, double this_weight){

  FillHist(hist_prefix+"Pt"+hist_suffix, this_lep.Pt(), this_weight, 0., 1000., 1000);
  FillHist(hist_prefix+"Eta"+hist_suffix, this_lep.Eta(), this_weight, -3., 3., 60);
  FillHist(hist_prefix+"RelIso"+hist_suffix, this_lep.RelIso(), this_weight, 0., 1., 100);
  FillHist(hist_prefix+"dXY"+hist_suffix, fabs(this_lep.dXY()), this_weight, 0., 0.01, 100);
  FillHist(hist_prefix+"dZ"+hist_suffix, fabs(this_lep.dZ()), this_weight, 0., 0.05, 100);
  FillHist(hist_prefix+"Flavour"+hist_suffix, this_lep.LeptonFlavour(), this_weight, 0., 3., 3);

}

double HNOSDiLepton::GetWeight(bool geterr, TString region, std::vector<snu::KMuon> muons, std::vector<snu::KElectron> electrons){

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
    if(region == "MuEl_OS" || region == "MuEl_SS"){
      if((PassTriggerOR(triggerlist_emBG1) || PassTriggerOR(triggerlist_emBG2))) trigger_ps_weight += WeightByTrigger(triggerlist_emBG1, TargetLumi);
      if(((PassTriggerOR(triggerlist_emH1) || PassTriggerOR(triggerlist_emH2)))) trigger_ps_weight += WeightByTrigger(triggerlist_emH1, TargetLumi);
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

void HNOSDiLepton::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNOSDiLepton::BeginCycle() throw( LQError ){
  
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

HNOSDiLepton::~HNOSDiLepton() {
  
  Message("In HNOSDiLepton Destructor" , INFO);
  
}


void HNOSDiLepton::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNOSDiLepton::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNOSDiLeptonCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNOSDiLepton::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


void HNOSDiLepton::GENSignalStudy( bool Is_Signal ){

  if( !(k_sample_name.Contains("HN"))) return;
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

/*  TruthPrintOut();
  return;*/

  int electron_ID = 11, muon_ID = 13, W_ID = 24, HN_ID = 9900012;
  int parton_ID = 1;

  //TruthPrintOut();

  int max = truthColl.size();
  vector<int> HN_index, onshell_W_index, lep1_index, lep2_index, quark1_index, quark2_index, tchannel_quark1_index, tchannel_quark2_index;
  HN_index.clear(); onshell_W_index.clear(); lep1_index.clear(); lep2_index.clear(); quark1_index.clear(); quark2_index.clear(); tchannel_quark1_index.clear(); tchannel_quark2_index.clear();

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
    if( abs(truthColl.at(i).PdgId()) < 5 && abs(truthColl.at(truthColl.at(i).IndexMother()).PdgId()) == W_ID ){
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
    if( abs(truthColl.at(i).PdgId()) < 5 && abs(truthColl.at(truthColl.at(i).IndexMother()).PdgId()) == W_ID ){
      quark2_index.push_back(i);
      GENFindDecayIndex( truthColl, i, quark2_index);
      break;
    }
  }
  if( quark2_index.size () == 0 ){
    FillHist("GEN_quark2_not_found", 0., 1., 0., 1., 1);
    return;
  }

  GENFindDecayIndex( truthColl, 0, tchannel_quark1_index);
  GENFindDecayIndex( truthColl, 1, tchannel_quark2_index);
  if( (k_sample_name.Contains("Tchannel")) ){
    if(tchannel_quark1_index.size() * tchannel_quark2_index.size() == 0) return;
  }

  snu::KParticle GEN_lep[2], GEN_quark[2], GEN_onshellW, GEN_HN, GEN_Tquark[2];
  GEN_lep[0] = truthColl.at( lep1_index.back() );
  GEN_lep[1] = truthColl.at( lep2_index.back() );
/*  if(GEN_lep[0].Pt() < GEN_lep[1].Pt()){
    GEN_lep[1] = truthColl.at( lep1_index.back() )
    GEN_lep[0] = truthColl.at( lep2_index.back() );
  }*/

  GEN_quark[0] = truthColl.at( quark1_index.back() );
  GEN_quark[1] = truthColl.at( quark2_index.back() );
  if(GEN_quark[0].Pt() < GEN_quark[1].Pt()){
    GEN_quark[1] = truthColl.at( quark1_index.back() );
    GEN_quark[0] = truthColl.at( quark2_index.back() );
  }

  if( (k_sample_name.Contains("Tchannel")) ){
    GEN_Tquark[0] = truthColl.at(tchannel_quark1_index.back());
    GEN_Tquark[1] = truthColl.at(tchannel_quark2_index.back());

    if(GEN_Tquark[0].Pt() < GEN_Tquark[1].Pt()){
      GEN_Tquark[1] = truthColl.at( tchannel_quark1_index.back() );
      GEN_Tquark[0] = truthColl.at( tchannel_quark2_index.back() );
    }

    FillHist("GEN_DiTQuark_Mass", (GEN_Tquark[0]+GEN_Tquark[1]).M(), 1., 0., 2000., 2000);
    FillHist("GEN_LeadingTQuark_Pt", GEN_Tquark[0].Pt(), 1., 0., 2000., 2000);
    FillHist("GEN_SubLeadingTQuark_Pt", GEN_Tquark[1].Pt(), 1., 0., 2000., 2000);
  }

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
  
  double dphi[2] = {0.,};
  dphi[0] = GEN_onshellW.DeltaPhi(GEN_lep[0]);
  dphi[1] = GEN_onshellW.DeltaPhi(GEN_lep[1]);
  if(fabs(dphi[0]-3.) > fabs(dphi[1]-3.)) FillHist("GEN_COMPARE_SecondaryLeptonAndonshellW_DeltaPhi_CloserTo3p0", 1., 1., 0., 2., 2);
  else FillHist("GEN_COMPARE_SecondaryLeptonAndonshellW_DeltaPhi_CloserTo3p0", 0., 1., 0., 2., 2);
  if(GEN_lep[0].Pt() > GEN_lep[1].Pt()) FillHist("GEN_COMPARE_PrimaryLeptonPt_BiggerThan_SecondaryLeptonPt", 1., 1., 0., 2., 2);
  else FillHist("GEN_COMPARE_PrimaryLeptonPt_BiggerThan_SecondaryLeptonPt", 0., 1., 0., 2., 2);

/*
  std::vector<snu::KMuon> muons = GetMuons("MUON_HN_TIGHT", false);
  std::vector<snu::KElectron> electrons = GetElectrons(false, false, "ELECTRON_HN_TIGHTv4");
  if((muons.size() + electrons.size()) != 2) return;
    
  std::vector<snu::KJet> jets = GetJets("JET_HN");
  if(jets.size() < 2) return;

  snu::KJet RECO_jets[2];
  RECO_jets[0] = jets.at(0);
  RECO_jets[1] = jets.at(1);
  int RECO_jets_index[2];
  bool GEN_RECO_Jets_Matched = false;
  for(int i=0; i<jets.size(); i++){
    if(jets.at(i).DeltaR(GEN_quark[0]) < 0.1){
      for(int j=0; j<jets.size(); j++){
        if(i != j){
          if(jets.at(j).DeltaR(GEN_quark[1]) < 0.1){

            if(GEN_RECO_Jets_Matched){
              double chi_2=0., new_chi_2=0.;
              chi_2=(RECO_jets[0].Pt()-GEN_quark[0].Pt())*(RECO_jets[0].Pt()-GEN_quark[0].Pt())/GEN_quark[0].Pt() + (RECO_jets[1].Pt()-GEN_quark[1].Pt())*(RECO_jets[1].Pt()-GEN_quark[1].Pt())/GEN_quark[1].Pt();
              if(i < j){
                new_chi_2=(jets.at(i).Pt()-GEN_quark[0].Pt())*(jets.at(i).Pt()-GEN_quark[0].Pt())/GEN_quark[0].Pt() + (jets.at(j).Pt()-GEN_quark[1].Pt())*(jets.at(j).Pt()-GEN_quark[1].Pt())/GEN_quark[1].Pt();
              }
              else{
                new_chi_2=(jets.at(j).Pt()-GEN_quark[0].Pt())*(jets.at(j).Pt()-GEN_quark[0].Pt())/GEN_quark[0].Pt() + (jets.at(i).Pt()-GEN_quark[1].Pt())*(jets.at(i).Pt()-GEN_quark[1].Pt())/GEN_quark[1].Pt();
              }
              if(chi_2 > new_chi_2){
                if(i < j){
                  RECO_jets[0] = jets.at(i);
                  RECO_jets[1] = jets.at(j);
                }
                else{
                  RECO_jets[1] = jets.at(i);
                  RECO_jets[0] = jets.at(j);
                }
                RECO_jets_index[0] = i;
                RECO_jets_index[1] = j;
              }
            }
            if(!GEN_RECO_Jets_Matched){
              GEN_RECO_Jets_Matched = true;
              if(i < j){
                RECO_jets[0] = jets.at(i);
                RECO_jets[1] = jets.at(j);
              }
              else{
                RECO_jets[1] = jets.at(i);
                RECO_jets[0] = jets.at(j);
              }
              RECO_jets_index[0] = i;
              RECO_jets_index[1] = j;
            }
          }
        }
      }
    }
  }
  if(!GEN_RECO_Jets_Matched) return;

  FillHist("GEN_RECO_MATCHED_Jets_Pt_Order", RECO_jets_index[0], 1., 0., 10., 10); 
  FillHist("GEN_RECO_MATCHED_Jets_Pt_Order", RECO_jets_index[1], 1., 0., 10., 10);

  FillHist("GEN_RECO_MATCHED_DiJet_Mass", (RECO_jets[0]+RECO_jets[1]).M(), 1., 0., 200., 200);
  FillHist("GEN_RECO_MATCHED_DiJet_DeltaR", (RECO_jets[0]).DeltaR(RECO_jets[1]), 1., 0., 5., 50);
  FillHist("GEN_RECO_MATCHED_DiJet_xMass_yDeltaR", (RECO_jets[0]+RECO_jets[1]).M(), (RECO_jets[0]).DeltaR(RECO_jets[1]), 1., 0., 200., 200, 0., 5., 50);
*/


/*
  snu::KParticle RECO_lep[2];

  int GEN_RECO_lep_matched = 1;
  if((GEN_lep[0].DeltaR(muons.at(0)) < 0.1) && (GEN_lep[1].DeltaR(muons.at(1)) < 0.1)){
    GEN_RECO_lep_matched *= -1;
    RECO_lep[0] = muons.at(0);
    RECO_lep[1] = muons.at(1);
  }      
  if((GEN_lep[0].DeltaR(muons.at(1)) < 0.1) && (GEN_lep[1].DeltaR(muons.at(0)) < 0.1)){
    GEN_RECO_lep_matched *= -1;
    RECO_lep[0] = muons.at(1);
    RECO_lep[1] = muons.at(0);
  }
  if( GEN_RECO_lep_matched == 1 ){
    FillHist("GEN_Lepton_Matching_Failed", 0., 1., 0., 1., 1);    
    return;
  }


  int GEN_RECO_jet_matched = -999;
  int RECO_jet_matched[2] = {-999, -999};
  for(int i=0; i<jets.size(); i++){
    if(GEN_quark[0].DeltaR(jets.at(i)) < 0.1){
      GEN_RECO_jet_matched = i;
      RECO_jet_matched[0] = i;
      for(int j=0; j<jets.size(); j++){
        if((GEN_quark[1].DeltaR(jets.at(j)) < 0.1) && GEN_RECO_jet_matched != j){
          RECO_jet_matched[1] = j;
        }
      }
    }
  }

  if(RECO_jet_matched[0] == -999 || RECO_jet_matched[1] == -999){
    FillHist("GEN_Jet_Matching_Failed", 0., 1., 0., 1., 1);
    return;
  }
  if(!(RECO_jet_matched[0] == 0 && RECO_jet_matched[1] == 1) && !(RECO_jet_matched[0] == 1 && RECO_jet_matched[1] == 0)) FillHist("GEN_Jet_Not_Leading", 0., 1., 0., 1., 1);

  snu::KParticle RECO_jet[2];
  RECO_jet[0] = jets.at(RECO_jet_matched[0]);
  RECO_jet[1] = jets.at(RECO_jet_matched[1]);

  if(RECO_jet[0].Pt() < RECO_jet[1].Pt()){
    RECO_jet[1] = jets.at(RECO_jet_matched[0]);
    RECO_jet[0] = jets.at(RECO_jet_matched[1]);
  }

  snu::KParticle RECO_onshellW, RECO_HN;
  RECO_onshellW = RECO_jet[0]+RECO_jet[1];
  RECO_HN = RECO_lep[1] + RECO_onshellW;

  FillHist("RECO_LeadingJet_Pt", RECO_jet[0].Pt(), 1., 0., 1000., 1000);
  FillHist("RECO_SubLeadingJet_Pt", RECO_jet[1].Pt(), 1., 0., 1000., 1000);
  FillHist("RECO_OnShellW_Mass", RECO_onshellW.M(), 1., 0., 500., 500);
  FillHist("RECO_HeavyNeutrino_Mass", RECO_HN.M(), 1., 0., 2500., 2500);

  FillHist("RECO_DiJet_DeltaR", RECO_jet[0].DeltaR(RECO_jet[1]), 1., 0., 5., 50);
*/
  return;  

}


void HNOSDiLepton::GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index ){

  for( int i = it+1 ; i < truthColl.size(); i ++ ){
    if( truthColl.at(i).IndexMother() == it && truthColl.at(i).PdgId() == truthColl.at(it).PdgId() ){
      index.push_back(i);
    }
  }
  return;

}


std::vector<snu::KElectron> HNOSDiLepton::ShiftElectronEnergy(std::vector<snu::KElectron> beforeshift, TString el_ID, bool applyshift){

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
