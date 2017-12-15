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

  if(std::find(k_flags.begin(), k_flags.end(), "truth") !=k_flags.end()) DoTruthMCStudy();

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
  std::vector<TString> triggerlist_mm_All, triggerlist_mm_H;
  std::vector<TString> triggerlist_em_All, triggerlist_em_BtoG, triggerlist_em_H, triggerlist_em_mu8el23, triggerlist_em_mu23el8;
  std::vector<TString> triggerlist_ee_All;
  triggerlist_mm_All.clear(); triggerlist_mm_H.clear();
  triggerlist_em_All.clear(); triggerlist_em_BtoG.clear(); triggerlist_em_H.clear(); triggerlist_em_mu8el23.clear(); triggerlist_em_mu23el8.clear();
  triggerlist_ee_All.clear();

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

  if(k_channel.Contains("DoubleMuon"))	if(!PassTriggerOR(triggerlist_mm_All)) return;
  if(k_channel.Contains("MuonEG"))	if(!PassTriggerOR(triggerlist_em_All)) return;
  if(k_channel.Contains("DoubleEG"))	if(!PassTriggerOR(triggerlist_ee_All)) return;
  // ================================================================================


  // ========== Get Objects (muon, electron, jet) ====================
  std::vector<KLepton> leptons_before_pt_order; leptons_before_pt_order.clear();

  std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_POG_LOOSE", true);
  std::vector<snu::KMuon> muons;  muons.clear();
  std::vector<snu::KElectron> electronVetoColl = GetElectrons(true, true, "ELECTRON_POG_LOOSE");
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
		if(!( abs(muon_sum_charge) > 2 || abs(electron_sum_charge) > 2 )) bool_OS_charge = true;

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
                if(!( abs(muon_sum_charge) > 2 || abs(electron_sum_charge) > 2 )) bool_OS_charge = true;

		if(electrons.at(0).Pt() > 25 && electrons.at(1).Pt() > 15) bool_lepton_pt = true;

		string_trigger = "EE";
        }
        else return;
  }

  if(!bool_OS_charge) return;
  if(!bool_lepton_pt) return;

/*  for(int i=0; i<leptons.size(); i++){
    cout<<leptons.at(i).Pt()<<endl;
  }*/

  std::vector<snu::KJet> jetsPreColl = GetJets("JET_NOLEPTONVETO", 30, 5.0);
  std::vector<snu::KJet> jets, bjets; jets.clear(); bjets.clear();
  
  for(int i=0; i<jetsPreColl.size(); i++){
    bool bool_away_from_lepton = true;
    snu::KJet this_jet=jetsPreColl.at(i);
    if( this_jet.IsBTagged(snu::KJet::CSVv2, snu::KJet::Loose) ) bjets.push_back(this_jet);
    for(int j=0; j<leptons.size(); j++){
      if( this_jet.DeltaR( leptons.at(j) ) < 0.4 ){
        bool_away_from_lepton = false;
	break;
      }

    }
    if( bool_away_from_lepton ) jets.push_back( this_jet );
  }

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
  if(string_trigger == "EE"){
    if(!isData){
      this_weight *= 35863.307;
    }
  }
  // ================================================================================


  TString this_string = "NULL";
  if(string_common == "2L"){

    this_string = string_lepton+"_preselection_";
    DrawHistColl(this_string, leptons, jets, bjets, MET, this_weight);
    FillHist(this_string+"ll_mass", (leptons.at(0)+leptons.at(1)).M(), this_weight, 0., 200., 200);
    FillHist(this_string+"ll_pt", (leptons.at(0)+leptons.at(1)).Pt(), this_weight, 0., 200., 200);
    FillHist(this_string+"ll_deltar", leptons.at(0).DeltaR(leptons.at(1)), this_weight, 0., 5., 50);
    if(bjets.size() == 0){
      this_string = string_lepton+"_bjets0_";
      DrawHistColl(this_string, leptons, jets, bjets, MET, this_weight);
      FillHist(this_string+"ll_mass", (leptons.at(0)+leptons.at(1)).M(), this_weight, 0., 200., 200);
      FillHist(this_string+"ll_pt", (leptons.at(0)+leptons.at(1)).Pt(), this_weight, 0., 200., 200);
      FillHist(this_string+"ll_deltar", leptons.at(0).DeltaR(leptons.at(1)), this_weight, 0., 5., 50);
      if(jets.size() < 2){
        this_string = string_lepton+"_jetssm2_";
        DrawHistColl(this_string, leptons, jets, bjets, MET, this_weight);
        FillHist(this_string+"ll_mass", (leptons.at(0)+leptons.at(1)).M(), this_weight, 0., 200., 200);
        FillHist(this_string+"ll_pt", (leptons.at(0)+leptons.at(1)).Pt(), this_weight, 0., 200., 200);
        FillHist(this_string+"ll_deltar", leptons.at(0).DeltaR(leptons.at(1)), this_weight, 0., 5., 50);
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


  FillHist(this_string+"number", leptons.size(), this_weight, 0., 5., 5);
  FillHist(this_string+"jet_number", jets.size(), this_weight, 0., 5., 5);
  FillHist(this_string+"bjet_number", bjets.size(), this_weight, 0., 5., 5);
  FillHist(this_string+"met", MET.Pt(), this_weight, 0., 400., 400);
  FillHist(this_string+"eventnumber", 0., this_weight, 0., 1., 1);

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


void WWAnalyzer::DoTruthMCStudy( void ){
  //TruthPrintOut();
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

  int max = truthColl.size();

  std::vector<int> wp_index, wm_index;
  for( int i=2 ; i<10 ; i++){
    if( ((truthColl.at(i).PdgId()) == 24) ){
      wp_index.push_back(i);
    }
    if( ((truthColl.at(i).PdgId()) == -24) ){
      wm_index.push_back(i);
    }
  }
//  cout<<wp_index.size()<<" "<<wm_index.size()<<endl;

  if(wp_index.size() * wm_index.size() ==0){  TruthPrintOut(); FillHist("?", 0., 1., 0., 1., 1); return; }

  snu::KParticle Wp, Wm, WW;
  Wp.SetPxPyPzE(truthColl.at(wp_index.at(0)).Px(), truthColl.at(wp_index.at(0)).Py(), truthColl.at(wp_index.at(0)).Pz(), truthColl.at(wp_index.at(0)).E());
  Wm.SetPxPyPzE(truthColl.at(wm_index.at(0)).Px(), truthColl.at(wm_index.at(0)).Py(), truthColl.at(wm_index.at(0)).Pz(), truthColl.at(wm_index.at(0)).E());
  WW = Wp + Wm;

  double P_Wp[2] = {0.,}, P_Wm[2] = {0.,};
  P_Wp[0] = (Wp.E() + Wp.Pz())/TMath::Sqrt(2.);
  P_Wp[1] = (Wp.E() - Wp.Pz())/TMath::Sqrt(2.);
  P_Wm[0] = (Wm.E() + Wm.Pz())/TMath::Sqrt(2.);
  P_Wm[1] = (Wm.E() - Wm.Pz())/TMath::Sqrt(2.);

  double CosCSangle = (WW.Pz()/fabs(WW.Pz())) * 2 * (P_Wp[0]*P_Wm[1] - P_Wp[1]*P_Wm[0]) / (WW.M() * TMath::Sqrt( WW.M()*WW.M() + WW.Pt()*WW.Pt() ));

  FillHist("CosCSangle", CosCSangle, 1., -1., 1., 50);
  FillHist("A_FB_massrange_Global", CosCSangle/fabs(CosCSangle), 1., -2., 3., 5);
  if(WW.M() > 400.){
    FillHist("A_FB_massrange_400toInf", CosCSangle/fabs(CosCSangle), 1., -2., 3., 5);
  }
  else if(WW.M() > 300.){
    FillHist("A_FB_massrange_300to400", CosCSangle/fabs(CosCSangle), 1., -2., 3., 5);
  }
  else if(WW.M() > 250.){
    FillHist("A_FB_massrange_250to300", CosCSangle/fabs(CosCSangle), 1., -2., 3., 5);
  }
  else if(WW.M() > 200.){
    FillHist("A_FB_massrange_200to250", CosCSangle/fabs(CosCSangle), 1., -2., 3., 5);
  }
  else{
    FillHist("A_FB_massrange_0to200", CosCSangle/fabs(CosCSangle), 1., -2., 3., 5);
  }


  FillHist("mass_WW", WW.M(), 1., 0., 500., 500);

  return;
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



