// $Id: WWAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQWWAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "WWAnalyzer.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (WWAnalyzer);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
WWAnalyzer::WWAnalyzer() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("WWAnalyzer");
  
  Message("In WWAnalyzer constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void WWAnalyzer::InitialiseAnalysis() throw( LQError ) {
  
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


  return;
}


void WWAnalyzer::ExecuteEvents()throw( LQError ){

  DoTruthMCStudy();

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
  // ========== Trigger cut ====================
/*  triggerlist_m.push_back("HLT_IsoMu24_v");
  triggerlist_m.push_back("HLT_IsoTkMu24_v");*/
  std::vector<TString> triggerlist_m_All, triggerlist_e_All;
  std::vector<TString> triggerlist_mm_All, triggerlist_mm_H;
  std::vector<TString> triggerlist_em_All, triggerlist_em_BtoG, triggerlist_em_H, triggerlist_em_mu8el23, triggerlist_em_mu23el8;
  std::vector<TString> triggerlist_ee_All;
  triggerlist_m_All.clear(); triggerlist_e_All.clear();
  triggerlist_mm_All.clear(); triggerlist_mm_H.clear();
  triggerlist_em_All.clear(); triggerlist_em_BtoG.clear(); triggerlist_em_H.clear(); triggerlist_em_mu8el23.clear(); triggerlist_em_mu23el8.clear();
  triggerlist_ee_All.clear();

  triggerlist_m_All.push_back("HLT_IsoMu24_v");
  triggerlist_m_All.push_back("HLT_IsoTkMu24_v");
  triggerlist_e_All.push_back("HLT_Ele27_WPTight_Gsf");
  triggerlist_mm_All.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  triggerlist_mm_All.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  triggerlist_mm_All.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_mm_All.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  triggerlist_em_All.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_em_All.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_em_All.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_em_All.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_ee_All.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_mm_H.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_mm_H.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  triggerlist_em_BtoG.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_em_BtoG.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_em_H.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_em_H.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_em_mu8el23.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_em_mu8el23.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_em_mu23el8.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_em_mu23el8.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");

  if(k_channel.Contains("SingleMuon"))		if(!PassTriggerOR(triggerlist_m_All))  return;
  if(k_channel.Contains("SingleElectron"))	if(!PassTriggerOR(triggerlist_e_All))  return;
  if(k_channel.Contains("DoubleMuon"))		if(!PassTriggerOR(triggerlist_mm_All)) return;
  if(k_channel.Contains("MuonEG"))		if(!PassTriggerOR(triggerlist_em_All)) return;
  if(k_channel.Contains("DoubleEG"))		if(!PassTriggerOR(triggerlist_ee_All)) return;
  // ================================================================================


  // ========== Get Objects (muon, electron, jet) ====================
  std::vector<KLepton> leptons_before_pt_order; leptons_before_pt_order.clear();

  std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_POG_LOOSE", false);
  std::vector<snu::KMuon> muons;  muons.clear();
  std::vector<snu::KElectron> electronVetoColl = GetElectrons(false, false, "ELECTRON_POG_LOOSE");
  std::vector<snu::KElectron> electrons;  electrons.clear();

  for(unsigned int i=0; i<muonVetoColl.size(); i++){
    if(PassID(muonVetoColl.at(i), "MUON_POG_TIGHT")){
      muons.push_back(muonVetoColl.at(i));
      leptons_before_pt_order.push_back(muonVetoColl.at(i));
    }
  }
  for(unsigned int i=0; i<electronVetoColl.size(); i++){
    if(PassID(electronVetoColl.at(i), "ELECTRON_POG_TIGHT_CHARGE")){
      electrons.push_back(electronVetoColl.at(i));
      leptons_before_pt_order.push_back(electronVetoColl.at(i));
    }
  }

  if(muonVetoColl.size() != muons.size()) return;
  if(electronVetoColl.size() != electrons.size()) return;
  TString string_lepton = "NULL";

  if(muons.size() == 1 && electrons.size() == 0) string_lepton = "Mu1El0";
  if(muons.size() == 0 && electrons.size() == 1) string_lepton = "Mu0El1";

  if(muons.size() == 2 && electrons.size() == 0) string_lepton = "Mu2El0";
  if(muons.size() == 1 && electrons.size() == 1) string_lepton = "Mu1El1";
  if(muons.size() == 0 && electrons.size() == 2) string_lepton = "Mu0El2";

  if(muons.size() == 3 && electrons.size() == 0) string_lepton = "Mu3El0";
  if(muons.size() == 2 && electrons.size() == 1) string_lepton = "Mu2El1";
  if(muons.size() == 1 && electrons.size() == 2) string_lepton = "Mu1El2";
  if(muons.size() == 0 && electrons.size() == 3) string_lepton = "Mu0El3";

  if(muons.size() == 4 && electrons.size() == 0) string_lepton = "Mu4El0";
  if(muons.size() == 2 && electrons.size() == 2) string_lepton = "Mu2El2";
  if(muons.size() == 0 && electrons.size() == 4) string_lepton = "Mu0El4";

  if(string_lepton == "NULL") return;
  
  TString string_common = "NULL";
  if(muons.size() + electrons.size() == 1) string_common = "1L";
  if(muons.size() + electrons.size() == 2) string_common = "2L";
  if(muons.size() + electrons.size() == 3) string_common = "3L";
  if(muons.size() + electrons.size() == 4) string_common = "4L";

  std::vector<KLepton> leptons = SortByPtOrder( leptons_before_pt_order );
  if(leptons.at(leptons.size()-1).Pt() < 10) return;
  for(int i=0; i<leptons.size(); i++){
    for(int j=i+1; j<leptons.size(); j++){
      if((leptons.at(i) + leptons.at(j)).M() < 10) return;
    }
  }

  CorrectMuonMomentum(muons);
  CorrectedMETRochester(muons);
  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  snu::KParticle MET;
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);

  bool bool_OS_charge = false;
  bool bool_lepton_pt = false;
  TString string_trigger = "NULL";

  if(string_lepton == "Mu1El0"){
	if(PassTriggerOR(triggerlist_m_All)){
		if(muons.at(0).Pt() > 27) bool_lepton_pt = true;
		bool_OS_charge = true;
		string_trigger = "M";
	}
	else return;
  }
  if(string_lepton == "Mu0El1"){
	if(PassTriggerOR(triggerlist_e_All)){
		if(electrons.at(0).Pt() > 30) bool_lepton_pt = true;
		bool_OS_charge = true;
		string_trigger = "E";
	}
	else return;
  }

  if(string_lepton == "Mu2El0" 
  || string_lepton == "Mu3El0"
  || string_lepton == "Mu2El1"
  || string_lepton == "Mu4El0"
  || string_lepton == "Mu2El2"){
	if(PassTriggerOR(triggerlist_mm_All)){
                int muon_sum_charge = 0, electron_sum_charge = 0;
		for(int i=0; i<muons.size(); i++){
		  muon_sum_charge += muons.at(i).Charge();
                }
                for(int i=0; i<electrons.size(); i++){
                  electron_sum_charge += electrons.at(i).Charge();
                }
		if(!( abs(muon_sum_charge) > 1 || abs(electron_sum_charge) > 1 )) bool_OS_charge = true;

		if(muons.at(0).Pt() > 20 && muons.at(1).Pt() > 10) bool_lepton_pt = true;

                string_trigger = "MM";
	}
	else return;
  }  
  if(string_lepton == "Mu1El1"){
        if(PassTriggerOR(triggerlist_em_All)){
                if( muons.at(0).Charge() + electrons.at(0).Charge() == 0 ) bool_OS_charge = true;

		if(PassTriggerOR(triggerlist_em_mu8el23)){
		  if(muons.at(0).Pt() > 10 && electrons.at(0).Pt() > 25) bool_lepton_pt = true;
		}
		if(PassTriggerOR(triggerlist_em_mu23el8)){
                  if(!bool_lepton_pt){
                    if(muons.at(0).Pt() > 25 && electrons.at(0).Pt() > 10) bool_lepton_pt = true;
                  }
                }

		string_trigger = "EM";
        }
	else return;
  }
  if(string_lepton == "Mu0El2"
  || string_lepton == "Mu1El2"
  || string_lepton == "Mu0El3"
  || string_lepton == "Mu0El4"){
        if(PassTriggerOR(triggerlist_ee_All)){
                int muon_sum_charge = 0, electron_sum_charge = 0;
                for(int i=0; i<muons.size(); i++){
                  muon_sum_charge += muons.at(i).Charge();
                }
                for(int i=0; i<electrons.size(); i++){
                  electron_sum_charge += electrons.at(i).Charge();
                }
                if(!( abs(muon_sum_charge) > 1 || abs(electron_sum_charge) > 1 )) bool_OS_charge = true;

		if(electrons.at(0).Pt() > 25 && electrons.at(1).Pt() > 15) bool_lepton_pt = true;

		string_trigger = "EE";
        }
        else return;
  }

  if(!bool_OS_charge) return;
  if(!bool_lepton_pt) return;

  std::vector<snu::KJet> jetsPreColl = GetJets("JET_NOLEPTONVETO", 20, 2.4);
  std::vector<snu::KJet> jets, bjets; jets.clear(); bjets.clear();
  std::vector<snu::KJet> bjets_loose, bjets_medium, bjets_tight;
  bjets_loose.clear(); bjets_medium.clear(); bjets_tight.clear();
  for(int i=0; i<jetsPreColl.size(); i++){
    bool bool_away_from_lepton = true;
    snu::KJet this_jet=jetsPreColl.at(i);
    if( this_jet.IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium) ) bjets.push_back(this_jet);
    if( this_jet.IsBTagged(snu::KJet::CSVv2, snu::KJet::Loose) ) bjets_loose.push_back(this_jet);
    if( this_jet.IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium) ) bjets_medium.push_back(this_jet);
    if( this_jet.IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight) ) bjets_tight.push_back(this_jet);

    for(int j=0; j<leptons.size(); j++){
      if( this_jet.DeltaR( leptons.at(j) ) < 0.4 ){
        bool_away_from_lepton = false;
	break;
      }

    }
    if( bool_away_from_lepton ) jets.push_back( this_jet );
  }
  for(int i=0; i<bjets_loose.size(); i++)
    FillHist("LOOSE", bjets_loose.at(i).Eta(), 1., -3., 3., 60);
  for(int i=0; i<bjets_medium.size(); i++)
    FillHist("MEDIUM", bjets_medium.at(i).Eta(), 1., -3., 3., 60);
  for(int i=0; i<bjets_tight.size(); i++)
    FillHist("TIGHT", bjets_tight.at(i).Eta(), 1., -3., 3., 60);
  // ================================================================================

  // Get event weight
  double this_weight = 1.;
  if(!isData){
    this_weight = weight;
    this_weight *= MCweight;  
    this_weight *= mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muons, 0);
    this_weight *= mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muons, 0);
    this_weight *= mcdata_correction->MuonTrackingEffScaleFactor(muons);
    this_weight *= mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT_CHARGE", electrons);
    this_weight *= mcdata_correction->ElectronRecoScaleFactor(electrons);
    this_weight *= mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    this_weight *= GetKFactor();
  }
 

  if(string_trigger == "MM"){
    if(isData){
      if(GetDataPeriod() == 7) if(!PassTriggerOR(triggerlist_mm_H)) return;
    }
    if(!isData){
      double temp_weight = 0.;
      if(PassTriggerOR(triggerlist_mm_All)) temp_weight += 27257.617;
      if(PassTriggerOR(triggerlist_mm_H)) temp_weight += 8605.69;
      this_weight *= temp_weight;
    }
  }
  if(string_trigger == "EM"){
    if(isData){
      if(GetDataPeriod() < 7)  if(!PassTriggerOR(triggerlist_em_BtoG)) return;
      if(GetDataPeriod() == 7) if(!PassTriggerOR(triggerlist_em_H))  return;
    }
    if(!isData){
      double temp_weight = 0.;
      if(PassTriggerOR(triggerlist_em_BtoG)) temp_weight += 27257.617;
      if(PassTriggerOR(triggerlist_em_H))  temp_weight += 8605.69;
      this_weight *= temp_weight;
    }
  }
  if(string_trigger == "EE" || string_trigger == "M" || string_trigger == "E"){
    if(!isData){
      this_weight *= 35863.307;
    }
  }

  std::vector<snu::KMuon> pre_muonAntiIsoColl = GetMuons("MUON_POG_ANTIISO", true);
  std::vector<snu::KMuon> muonAntiColl; muonAntiColl.clear();
  for(int i=0; i<pre_muonAntiIsoColl.size(); i++){
    snu::KMuon this_muon = pre_muonAntiIsoColl.at(i);
    if(this_muon.RelIso04() > 0.25) muonAntiColl.push_back(this_muon);
  }
  // ================================================================================


  double projectedMET = GetProjectedMET(MET, leptons);
  int lepton_sum_charge = 0;
  double ST = 0.;
  int Wjet_index[2] = {0, 1};
  snu::KParticle Wleptonic[2];
  snu::KParticle WWdiboson;
  snu::KParticle MET_pz_recovered[2];
  bool imaginary_solution = false;

  for(int i=0; i<leptons.size(); i++){
    lepton_sum_charge += leptons.at(i).Charge();
    ST += leptons.at(i).Pt();
  }
  for(int i=0; i<jets.size(); i++){
    ST += jets.at(i).Pt();
  }
  ST += MET.Pt();

  TString this_string = "NULL";

  if(string_common == "2L"){
    std::vector<snu::KParticle> MT2assistedMET; MT2assistedMET.clear();
    double MT2 = GetMT2( leptons, MET, MT2assistedMET );

    bool img[2] = {false, false};
    std::vector<snu::KParticle> MT2assistedMETreco_p; MT2assistedMETreco_p.clear();
    std::vector<snu::KParticle> MT2assistedMETreco_m; MT2assistedMETreco_m.clear();

    GetNuPzFromWMass( leptons.at(0), MT2assistedMET.at(0), MT2assistedMETreco_p, img[0] );
    GetNuPzFromWMass( leptons.at(1), MT2assistedMET.at(1), MT2assistedMETreco_m, img[1] );

    this_string = string_lepton+"_preselection_";
    DrawHistColl(this_string, leptons, jets, bjets, MET, this_weight);
    DrawMAOSDifference(this_string, MT2assistedMETreco_p.at(0), MT2assistedMETreco_m.at(0), this_weight);
    FillHist(this_string+"ll_mass", (leptons.at(0)+leptons.at(1)).M(), this_weight, 0., 200., 200);
    FillHist(this_string+"ll_pt", (leptons.at(0)+leptons.at(1)).Pt(), this_weight, 0., 200., 200);
    FillHist(this_string+"ll_deltar", leptons.at(0).DeltaR(leptons.at(1)), this_weight, 0., 5., 50);
    FillHist(this_string+"antiisomuon_number", muonAntiColl.size(), this_weight, 0., 10., 10);
    FillHist(this_string+"projectedmet", projectedMET, this_weight, 0., 200., 200);
    FillHist(this_string+"chargesum", lepton_sum_charge, this_weight, -2., 3., 5);
    FillHist(this_string+"transversemass2", MT2, this_weight, 0., 500., 500);
    FillHist(this_string+"mt2assistedmet_metp_pt", MT2assistedMET.at(0).Pt(), this_weight, 0., 500., 500);
    FillHist(this_string+"mt2assistedmet_metm_pt", MT2assistedMET.at(1).Pt(), this_weight, 0., 500., 500);
    FillHist(this_string+"mt2assistedmet_metp_pzsol", MT2assistedMETreco_p.at(0).Pz(), this_weight, -200., 200., 400);
    FillHist(this_string+"mt2assistedmet_metm_pzsol", MT2assistedMETreco_m.at(0).Pz(), this_weight, -200., 200., 400);

    if(bjets.size() == 0){

      this_string = string_lepton+"_bjets0_";
      DrawHistColl(this_string, leptons, jets, bjets, MET, this_weight);
      DrawMAOSDifference(this_string, MT2assistedMETreco_p.at(0), MT2assistedMETreco_m.at(0), this_weight);
      FillHist(this_string+"ll_mass", (leptons.at(0)+leptons.at(1)).M(), this_weight, 0., 200., 200);
      FillHist(this_string+"ll_pt", (leptons.at(0)+leptons.at(1)).Pt(), this_weight, 0., 200., 200);
      FillHist(this_string+"ll_deltar", leptons.at(0).DeltaR(leptons.at(1)), this_weight, 0., 5., 50);
      FillHist(this_string+"antiisomuon_number", muonAntiColl.size(), this_weight, 0., 10., 10);
      FillHist(this_string+"projectedmet", projectedMET, this_weight, 0., 200., 200);
      FillHist(this_string+"chargesum", lepton_sum_charge, this_weight, -2., 3., 5);
      FillHist(this_string+"transversemass2", MT2, this_weight, 0., 500., 500);
      FillHist(this_string+"mt2assistedmet_metp_pt", MT2assistedMET.at(0).Pt(), this_weight, 0., 500., 500);
      FillHist(this_string+"mt2assistedmet_metm_pt", MT2assistedMET.at(1).Pt(), this_weight, 0., 500., 500);
      FillHist(this_string+"mt2assistedmet_metp_pzsol", MT2assistedMETreco_p.at(0).Pz(), this_weight, -200., 200., 400);
      FillHist(this_string+"mt2assistedmet_metm_pzsol", MT2assistedMETreco_m.at(0).Pz(), this_weight, -200., 200., 400);

      if(jets.size() < 2){

        this_string = string_lepton+"_jetssm2_";
        DrawHistColl(this_string, leptons, jets, bjets, MET, this_weight);
        DrawMAOSDifference(this_string, MT2assistedMETreco_p.at(0), MT2assistedMETreco_m.at(0), this_weight);
        FillHist(this_string+"ll_mass", (leptons.at(0)+leptons.at(1)).M(), this_weight, 0., 200., 200);
        FillHist(this_string+"ll_pt", (leptons.at(0)+leptons.at(1)).Pt(), this_weight, 0., 200., 200);
        FillHist(this_string+"ll_deltar", leptons.at(0).DeltaR(leptons.at(1)), this_weight, 0., 5., 50);
        FillHist(this_string+"antiisomuon_number", muonAntiColl.size(), this_weight, 0., 10., 10);
        FillHist(this_string+"projectedmet", projectedMET, this_weight, 0., 200., 200);
        FillHist(this_string+"chargesum", lepton_sum_charge, this_weight, -2., 3., 5);
        FillHist(this_string+"transversemass2", MT2, this_weight, 0., 500., 500);
        FillHist(this_string+"assistedmet_plus", MT2assistedMET.at(0).Pt(), this_weight, 0., 500., 500);
        FillHist(this_string+"assistedmet_minus", MT2assistedMET.at(1).Pt(), this_weight, 0., 500., 500);
        FillHist(this_string+"mt2assistedmet_metp_pt", MT2assistedMET.at(0).Pt(), this_weight, 0., 500., 500);
        FillHist(this_string+"mt2assistedmet_metm_pt", MT2assistedMET.at(1).Pt(), this_weight, 0., 500., 500);
        FillHist(this_string+"mt2assistedmet_metp_pzsol", MT2assistedMETreco_p.at(0).Pz(), this_weight, -200., 200., 400);
        FillHist(this_string+"mt2assistedmet_metm_pzsol", MT2assistedMETreco_m.at(0).Pz(), this_weight, -200., 200., 400);

      }
    }
  }


  return;
}// End of execute event loop

void WWAnalyzer::DrawHistColl( TString this_string,
			       std::vector<KLepton> leptons,
			       std::vector<snu::KJet> jets, std::vector<snu::KJet> bjets,
			       snu::KParticle MET,
			       double this_weight ){

  TString string_order[4] = {"1st_", "2nd_", "3rd_", "4th_"};
  TString this_obj_string = "";
  this_obj_string = this_string+"lepton_";
  for(int i=0; i<leptons.size(); i++){
    FillHist(this_obj_string+string_order[i]+"pt", leptons.at(i).Pt(), this_weight, 0., 400., 400);
    FillHist(this_obj_string+string_order[i]+"eta", leptons.at(i).Eta(), this_weight, -5., 5., 200);
    FillHist(this_obj_string+string_order[i]+"dxy", leptons.at(i).dXY(), this_weight, -1., 1., 200);
    FillHist(this_obj_string+string_order[i]+"dz", leptons.at(i).dZ(), this_weight, -1., 1., 200);
    FillHist(this_obj_string+string_order[i]+"reliso", leptons.at(i).RelIso(), this_weight, 0., 1., 200);
    FillHist(this_obj_string+string_order[i]+"flavour", leptons.at(i).LeptonFlavour(), this_weight, 0., 3., 3);
    FillHist(this_obj_string+"pt", leptons.at(i).Pt(), this_weight, 0., 400., 400);
    FillHist(this_obj_string+"eta", leptons.at(i).Eta(), this_weight, -5., 5., 200);
    FillHist(this_obj_string+"dxy", leptons.at(i).dXY(), this_weight, -1., 1., 200);
    FillHist(this_obj_string+"dz", leptons.at(i).dZ(), this_weight, -1., 1., 200);
    FillHist(this_obj_string+"reliso", leptons.at(i).RelIso(), this_weight, 0., 1., 200);
    FillHist(this_obj_string+"flavour", leptons.at(i).LeptonFlavour(), this_weight, 0., 3., 3);
  }

  this_obj_string = this_string+"jet_";
  for(int i=0; i<jets.size(); i++){
    if(i<4){
      FillHist(this_obj_string+string_order[i]+"pt", jets.at(i).Pt(), this_weight, 0., 400., 400);
      FillHist(this_obj_string+string_order[i]+"eta", jets.at(i).Eta(), this_weight, -5., 5., 200);
    }
    FillHist(this_obj_string+"pt", jets.at(i).Pt(), this_weight, 0., 400., 400);
    FillHist(this_obj_string+"eta", jets.at(i).Eta(), this_weight, -5., 5., 200);
  }

  this_obj_string = this_string+"bjet_";
  for(int i=0; i<bjets.size(); i++){
    if(i<4){
      FillHist(this_obj_string+string_order[i]+"pt", bjets.at(i).Pt(), this_weight, 0., 400., 400);
      FillHist(this_obj_string+string_order[i]+"eta", bjets.at(i).Eta(), this_weight, -5., 5., 200);
    }
    FillHist(this_obj_string+"pt", bjets.at(i).Pt(), this_weight, 0., 400., 400);
    FillHist(this_obj_string+"eta", bjets.at(i).Eta(), this_weight, -5., 5., 200);
  }

  FillHist(this_string+"event_number", 0., this_weight, 0., 1., 1);
  FillHist(this_string+"lepton_number", leptons.size(), this_weight, 0., 5., 5);
  FillHist(this_string+"jet_number", jets.size(), this_weight, 0., 5., 5);
  FillHist(this_string+"bjet_number", bjets.size(), this_weight, 0., 5., 5);
  FillHist(this_string+"met", MET.Pt(), this_weight, 0., 400., 400);

}

std::vector<KLepton> WWAnalyzer::SortByPtOrder( std::vector<KLepton> leptons ){
  int const N=leptons.size();
  int order[N];
  for(int i=0; i<N; i++){
    order[i] = i;
  }

  std::vector<KLepton> leptons_temp;
  bool can_stop = true;
  for(int i=0; i<leptons.size()-1; i++){
    KLepton this_lepton = leptons.at(i);
    KLepton next_lepton = leptons.at(i+1);

    if(this_lepton.Pt()<next_lepton.Pt()){
      order[i] = i+1;
      order[i+1] = i;
      can_stop = false;
    }
  }
  for(int i=0; i<leptons.size(); i++){
    leptons_temp.push_back(leptons.at(order[i]));
  }
  if(!can_stop) SortByPtOrder( leptons_temp );
  else return leptons_temp;
}

void WWAnalyzer::DrawMAOSDifference(TString this_string, snu::KParticle MT2assistedMETreco_p, snu::KParticle MT2assistedMETreco_m, double this_weight){

  if( !(k_sample_name.Contains("WWtoEMu"))) return;

  double diff_pt_p = (MT2assistedMETreco_p.Pt() - gen_nup.Pt())/gen_nup.Pt();
  double diff_pt_m = (MT2assistedMETreco_m.Pt() - gen_num.Pt())/gen_num.Pt();

  cout<< MT2assistedMETreco_p.Pt() << "\t" << gen_nup.Pt() <<endl;

  FillHist(this_string+"diff_pt_p", diff_pt_p, this_weight, -10., 10., 200);
  FillHist(this_string+"diff_pt_m", diff_pt_m, this_weight, -10., 10., 200);


}



void WWAnalyzer::DoTruthMCStudy( void ){

  //TruthPrintOut();
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

  int max = truthColl.size();
  std::vector<int> wp_vector, wm_vector; wp_vector.clear(); wm_vector.clear();
  for( int i=2; i<max ; i++){
    if( (truthColl.at(i).PdgId() == -24) &&
	((truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 21) || abs(truthColl.at(truthColl.at(i).IndexMother()).PdgId()) < 6) ){
		wp_vector.push_back(i);
		GENFindDecayIndex( truthColl, i, wp_vector );
		break;
    }
  }
  for( int i=2; i<max ; i++){
    if( (truthColl.at(i).PdgId() == 24) &&
        ((truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 21) || abs(truthColl.at(truthColl.at(i).IndexMother()).PdgId()) < 6) ){
		wm_vector.push_back(i);
                GENFindDecayIndex( truthColl, i, wm_vector );
		break;
    }
  }
  if(wp_vector.size() * wm_vector.size() == 0) return;
  int wp_index = wp_vector.back(), wm_index = wm_vector.back();

  std::vector<int> wpdecay_lep_vector, wmdecay_lep_vector, wpdecay_nu_vector, wmdecay_nu_vector;
  wpdecay_lep_vector.clear(); wmdecay_lep_vector.clear(); wpdecay_nu_vector.clear(); wmdecay_nu_vector.clear();

  for( int i=2; i<max ; i++){
    if( (abs(truthColl.at(i).PdgId()) == 11 || abs(truthColl.at(i).PdgId()) == 13) &&
        (truthColl.at(i).IndexMother() == wp_index) ){
                wpdecay_lep_vector.push_back(i);
                GENFindDecayIndex( truthColl, i, wpdecay_lep_vector );
                break;
    }
  }

  for( int i=2; i<max ; i++){
    if( (abs(truthColl.at(i).PdgId()) == 12 || abs(truthColl.at(i).PdgId()) == 14) &&
        (truthColl.at(i).IndexMother() == wp_index) ){
                wpdecay_nu_vector.push_back(i);
                GENFindDecayIndex( truthColl, i, wpdecay_nu_vector );
                break;
    }
  }  

  for( int i=2; i<max ; i++){
    if( (abs(truthColl.at(i).PdgId()) == 11 || abs(truthColl.at(i).PdgId()) == 13) &&
        (truthColl.at(i).IndexMother() == wm_index) ){
                wmdecay_lep_vector.push_back(i);
                GENFindDecayIndex( truthColl, i, wmdecay_lep_vector );
                break;
    }
  }

  for( int i=2; i<max ; i++){
    if( (abs(truthColl.at(i).PdgId()) == 12 || abs(truthColl.at(i).PdgId()) == 14) &&
        (truthColl.at(i).IndexMother() == wm_index) ){
                wmdecay_nu_vector.push_back(i);
                GENFindDecayIndex( truthColl, i, wmdecay_nu_vector );
                break;
    }
  }
  if(wpdecay_lep_vector.size() * wmdecay_lep_vector.size() * wpdecay_nu_vector.size() * wmdecay_nu_vector.size() == 0) return;
  int lepp_index = wpdecay_lep_vector.at(0), lepm_index = wmdecay_lep_vector.at(0);
  int nup_index = wpdecay_nu_vector.at(0), num_index = wmdecay_nu_vector.at(0);

  snu::KParticle gen_wp, gen_wm, gen_lepp, gen_lepm;//, gen_nup, gen_num;
  std::vector<snu::KParticle> gen_leptons;
  gen_lepp = truthColl.at( lepp_index );
  gen_nup = truthColl.at( nup_index );
  gen_lepm = truthColl.at( lepm_index );
  gen_num = truthColl.at( num_index );
  gen_wp = gen_lepp+gen_nup;
  gen_wm = gen_lepm+gen_num;
  gen_leptons.push_back( gen_lepp );
  gen_leptons.push_back( gen_lepm );

  std::vector<snu::KParticle> MT2_assisted_nu; MT2_assisted_nu.clear();
  double MT2 = GetMT2( gen_leptons, (gen_nup+gen_num), MT2_assisted_nu );

  FillHist("MCTRUTH_MT2", MT2, weight, 0., 200., 200);
  FillHist("MCTRUTH_gennup_pt", gen_nup.Pt(), weight, 0., 500., 500);
  FillHist("MCTRUTH_gennum_pt", gen_num.Pt(), weight, 0., 500., 500);
  FillHist("MCTRUTH_mt2assistednup_pt", MT2_assisted_nu.at(0).Pt(), weight, 0., 500., 500);
  FillHist("MCTRUTH_mt2assistednum_pt", MT2_assisted_nu.at(1).Pt(), weight, 0., 500., 500);
  FillHist("MCTRUTH_dev_nup",(MT2_assisted_nu.at(0).Pt() - gen_nup.Pt())/gen_nup.Pt(), weight, -50., 50., 200);
  FillHist("MCTRUTH_dev_num",(MT2_assisted_nu.at(1).Pt() - gen_num.Pt())/gen_num.Pt(), weight, -50., 50., 200);

  FillHist("MCTRUTH_PZSOLUTION_pzptratio",fabs(MT2_assisted_nu.at(0).Pz()/MT2_assisted_nu.at(0).Pt()), weight, 0., 10., 1000);
  FillHist("MCTRUTH_PZSOLUTION_pzptratio",fabs(MT2_assisted_nu.at(1).Pz()/MT2_assisted_nu.at(1).Pt()), weight, 0., 10., 1000);
  FillHist("MCTRUTH_PZSOLUTION_dev_pzptratio",fabs(fabs(MT2_assisted_nu.at(0).Pz()/MT2_assisted_nu.at(0).Pt())-fabs(MT2_assisted_nu.at(1).Pz()/MT2_assisted_nu.at(1).Pt())), weight, 0., 10., 1000);

  FillHist("MCTRUTH_lepton_pt", gen_lepp.Pt(), weight, 0., 500., 500);
  FillHist("MCTRUTH_lepton_pt", gen_lepm.Pt(), weight, 0., 500., 500);
  FillHist("MCTRUTH_leptonp_pt", gen_lepp.Pt(), weight, 0., 500., 500);
  FillHist("MCTRUTH_leptonm_pt", gen_lepm.Pt(), weight, 0., 500., 500);
  FillHist("MCTRUTH_lepton_eta", gen_lepp.Eta(), weight, -4., 4., 160);
  FillHist("MCTRUTH_lepton_eta", gen_lepm.Eta(), weight, -4., 4., 160);
  FillHist("MCTRUTH_leptonp_eta", gen_lepp.Eta(), weight, -4., 4., 160);
  FillHist("MCTRUTH_leptonm_eta", gen_lepm.Eta(), weight, -4., 4., 160);

  FillHist("MCTRUTH_Wp_pt", (gen_lepp+gen_nup).Pt(), weight, 0., 500., 500);
  FillHist("MCTRUTH_Wm_pt", (gen_lepm+gen_num).Pt(), weight, 0., 500., 500);
  FillHist("MCTRUTH_Wp_pz", (gen_lepp+gen_nup).Pz(), weight, -200., 200., 400);
  FillHist("MCTRUTH_Wm_pz", (gen_lepm+gen_num).Pz(), weight, -200., 200., 400);
  FillHist("MCTRUTH_Wp_eta", (gen_lepp+gen_nup).Eta(), weight, -4., 4., 160);
  FillHist("MCTRUTH_Wm_eta", (gen_lepm+gen_num).Eta(), weight, -4., 4., 160);

  FillHist("MCTRUTH_WpWm_pt", (gen_lepp+gen_nup+gen_lepm+gen_num).Pt(), weight, 0., 1000., 1000);
  FillHist("MCTRUTH_WpWm_pz", (gen_lepp+gen_nup+gen_lepm+gen_num).Pz(), weight, -500., 500., 1000);
  FillHist("MCTRUTH_WpWm_mass", (gen_lepp+gen_nup+gen_lepm+gen_num).M(), weight, 0., 3000., 3000);

  DrawAsymmetryHist("MCTRUTH_", (gen_lepp+gen_nup), (gen_lepm+gen_num), weight );

/* 
  //FillHist("2dhist", (gen_lep+gen_nu+gen_q1+gen_q2).M(), (reco_lep+MET+reco_q1+reco_q2).M(), 1., massarray, 8, massarray, 8);
*/
  return;
}

void WWAnalyzer::GENFindDecayIndex( std::vector<snu::KTruth> truthColl, int it, std::vector<int>& index ){

  for( int i=it+1; i<truthColl.size(); i++ ){
    if( truthColl.at(i).IndexMother() == index.back() && truthColl.at(i).PdgId() == truthColl.at(it).PdgId() ){
      index.push_back(i);
    }
  }
  return;

}

void WWAnalyzer::DrawAsymmetryHist( TString this_string, snu::KParticle wp, snu::KParticle wm, double this_weight ){

  snu::KParticle wpwm = wp + wm;
  double Pplus[2] = {0.,}, Pminus[2] = {0.,};

  Pplus[0] = ( wp.E() + wm.Pz() ) / TMath::Sqrt(2.);
  Pplus[1] = ( wp.E() - wm.Pz() ) / TMath::Sqrt(2.);
  Pminus[0] = ( wm.E() - wp.Pz() ) / TMath::Sqrt(2.);
  Pminus[1] = ( wm.E() - wp.Pz() ) / TMath::Sqrt(2.);

  double CosCSframe = 2.* (Pplus[0]*Pminus[1] - Pminus[0]*Pplus[1]) / (wpwm.M() * TMath::Sqrt( wpwm.M()*wpwm.M() + wpwm.Pt()*wpwm.Pt() ) );
  CosCSframe = CosCSframe * fabs(wpwm.Pz())/wpwm.Pz();

  int bin_size = 30;
  for(int i=0; i<bin_size; i++){
    if(i*(3000/bin_size) < wpwm.M())
    if(wpwm.M() < (i+1)*(3000/bin_size)){
	TString string_mass[2];
	string_mass[0].Form("%d",(i*(3000/bin_size)));        string_mass[1].Form("%d",((i+1)*(3000/bin_size)));
	TString string_massrange = "massrange"+string_mass[0]+"to"+string_mass[1];
	FillHist(this_string+"AFB_mean_"+string_massrange, fabs(CosCSframe)/CosCSframe, this_weight, -2., 3., 5);
	FillHist(this_string+"AFB_diboson_mass_"+string_massrange, wpwm.M(), this_weight, 0., 3000., 3000);
	FillHist(this_string+"AFB_event_number_"+string_massrange, 0., this_weight, 0., 1., 1);
	break;
    }
  }

  FillHist(this_string+"AFB_mean_massrangeGLOBAL", fabs(CosCSframe)/CosCSframe, this_weight, -2., 3., 5);
  FillHist(this_string+"AFB_diboson_mass_massrangeGLOBAL", wpwm.M(), this_weight, 0., 3000., 3000);
  FillHist(this_string+"AFB_event_number_massrangeGLOBAL", 0., this_weight, 0., 1., 1);

}


double WWAnalyzer::GetProjectedMET( snu::KParticle MET, std::vector<KLepton> leptons ){

  if(leptons.size() == 0) return 0.;

  double mindphi = (3.141592)/2.;
  for(int i=0; i<leptons.size(); i++){
    if( mindphi > fabs(leptons.at(i).DeltaPhi(MET)) ) mindphi = fabs(leptons.at(i).DeltaPhi(MET));
  }
  return (MET.Pt() * TMath::Sin(mindphi));


}

void WWAnalyzer::GetDijetMassClosest( std::vector<snu::KJet> jets, double target_mass, int& m, int& n){

  if(jets.size() < 2) return;

  double target_candidate_mass = (jets.at(0) + jets.at(1)).M();
  double target_temp_mass = 0.;;
  for(int i=0; i<jets.size(); i++){
    for(int j=i+1; j<jets.size(); j++){
      target_temp_mass = (jets.at(i) + jets.at(j)).M();
      if( fabs(target_temp_mass - target_mass) < fabs(target_candidate_mass - target_mass) ){
        m = i;
        n = j;
        target_candidate_mass = target_temp_mass;
      }
    }
  }

  return;
}

void WWAnalyzer::GetNuPzFromWMass( KLepton lepton, snu::KParticle MET, std::vector<snu::KParticle>& reconstructedMET, bool& img ){

  double m = 0., n = 0.;
  double target_mass = 80.4;
  double d = target_mass*target_mass - lepton.M()*lepton.M() + 2.* (lepton.Px()*MET.Px() + lepton.Py()*MET.Py());
  double a = lepton.E()*lepton.E() - lepton.Pz()*lepton.Pz();
  double b = d * lepton.Pz();
  double c = lepton.E()*lepton.E()*MET.E()*MET.E() - d*d/4.;
  if(b*b-4*a*c<0){
    img = true;
    m = b / (2.*a);
    n = b / (2.*a);
  }
  else{
    img = false;
    m =( b + TMath::Sqrt( b*b - 4.*a*c) ) / (2.*a);
    n =( b - TMath::Sqrt( b*b - 4.*a*c) ) / (2.*a);
  }

  if( fabs(n) < fabs(m) ){
    double temp = n;
    n = m;
    m = temp;
  }

  snu::KParticle temp[2];
  temp[0].SetPxPyPzE(MET.Px(), MET.Py(), m, TMath::Sqrt(MET.Px()*MET.Px() + MET.Py()*MET.Py() + m*m));
  temp[1].SetPxPyPzE(MET.Px(), MET.Py(), n, TMath::Sqrt(MET.Px()*MET.Px() + MET.Py()*MET.Py() + n*n));

  reconstructedMET.push_back(temp[0]);
  reconstructedMET.push_back(temp[1]);

}

double WWAnalyzer::GetMT( KLepton this_lepton, double this_MET ){

  double MT = 0.;
  MT = TMath::Sqrt( this_lepton.M()*this_lepton.M() + 0. + 2 * ( this_MET * (this_lepton.E() - this_lepton.Pt()) ));

  return MT;
}

double WWAnalyzer::GetMT( snu::KParticle this_lepton, double this_MET ){

  double MT = 0.;
  MT = TMath::Sqrt( this_lepton.M()*this_lepton.M() + 0. + 2 * ( this_MET * (this_lepton.E() - this_lepton.Pt()) ));

  return MT;
}

double WWAnalyzer::GetMT2( std::vector<KLepton> leptons, snu::KParticle MET, std::vector<snu::KParticle>& MT2_assisted_nu ){

  double MT2 = 99999999.;
  double MT2_cand = 0.;
  double METp_cand = 0.;
  int it_break = 0; 
 
  for(int it_MET = 0.; it_MET < MET.Pt(); it_MET++){

    MT2_cand = TMath::Max( GetMT(leptons.at(0), it_MET), GetMT(leptons.at(1), MET.Pt() - it_MET) );
    if(MT2 > MT2_cand){
      MT2 = MT2_cand;
      METp_cand = it_MET;
    }
    else{
      it_break++;
      if(it_break > 4) break;
    }
  }

  snu::KParticle MT2_assisted_nup;
  MT2_assisted_nup.SetPtEtaPhiM( METp_cand, 0., MET.Phi(), 0. );
  MT2_assisted_nu.push_back( MT2_assisted_nup );
  MT2_assisted_nu.push_back( (MET - MT2_assisted_nup) );

  return MT2;

}

double WWAnalyzer::GetMT2( std::vector<snu::KParticle> leptons, snu::KParticle MET, std::vector<snu::KParticle>& MT2_assisted_nu ){

  double MT2 = 99999999.;
  double MT2_cand = 0.;
  double METp_cand = 0.;
  int it_break = 0;

  for(int it_MET = 0.; it_MET < MET.Pt(); it_MET++){

    MT2_cand = TMath::Max( GetMT(leptons.at(0), it_MET), GetMT(leptons.at(1), MET.Pt() - it_MET) );
    if(MT2 > MT2_cand){
      MT2 = MT2_cand;
      METp_cand = it_MET;
    }
    else{
      it_break++;
      if(it_break > 4) break;
    }
  }

  snu::KParticle MT2_assisted_nup;
  MT2_assisted_nup.SetPtEtaPhiM( METp_cand, 0., MET.Phi(), 0. );
  MT2_assisted_nu.push_back( MT2_assisted_nup );
  MT2_assisted_nu.push_back( (MET - MT2_assisted_nup) );

  return MT2;

}


void WWAnalyzer::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void WWAnalyzer::BeginCycle() throw( LQError ){
  
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

WWAnalyzer::~WWAnalyzer() {
  
  Message("In WWAnalyzer Destructor" , INFO);
  
}


void WWAnalyzer::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void WWAnalyzer::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this WWAnalyzerCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void WWAnalyzer::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



