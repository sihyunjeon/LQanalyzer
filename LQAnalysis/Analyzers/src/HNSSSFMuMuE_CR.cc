// $Id: HNSSSFMuMuE_CR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNSSSFMuMuE_CR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNSSSFMuMuE_CR.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNSSSFMuMuE_CR);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNSSSFMuMuE_CR::HNSSSFMuMuE_CR() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNSSSFMuMuE_CR");
  
  Message("In HNSSSFMuMuE_CR constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNSSSFMuMuE_CR::InitialiseAnalysis() throw( LQError ) {
  
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

  string lqdir = getenv("LQANALYZER_DIR");

  return;
}


void HNSSSFMuMuE_CR::ExecuteEvents()throw( LQError ){
  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  FillCutFlow("NoCut", weight);

  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);

  if(!PassMETFilter()) return;     /// Initial event cuts :
  FillCutFlow("EventCut", weight);

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex
  float pileup_reweight=(1.0);
  if(!k_isdata){ pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);}

  std::vector<TString> triggerlist_mumu, triggerlist_muel;
  triggerlist_mumu.clear();
  triggerlist_muel.clear();

  bool pass_mumu_trig = false;
  bool pass_elel_trig = false;
  bool pass_muel_trig = false;

  triggerlist_mumu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_mumu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  if(PassTriggerOR(triggerlist_mumu)) pass_mumu_trig = true;

  if(PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) pass_elel_trig = true;

  triggerlist_muel.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_muel.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_muel.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if(PassTriggerOR(triggerlist_muel)) pass_muel_trig = true;

  if( !pass_mumu_trig && !pass_elel_trig) return;

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_HN", 30., 2.4);

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TRI_TIGHT",false);

  std::vector<snu::KElectron> electronLooseColl = GetElectrons(true,false,"ELECTRON16_POG_FAKELOOSE_CC_d0");
  std::vector<snu::KElectron> electronTightColl = GetElectrons(true,false,"ELECTRON16_FR_POG_TIGHT_CC");

  CorrectMuonMomentum(muonLooseColl);
  float muon_trkeff = mcdata_correction->MuonTrackingEffScaleFactor(muonLooseColl);
  float electron_idsf = mcdata_correction->ElectronScaleFactor("ELECTRON16_FR_POG_TIGHT_CC", electronLooseColl);
  float electron_reco = mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);

  if(!isData){
    weight *= pileup_reweight;
    weight *= electron_idsf;
    weight *= electron_reco;
    weight *= muon_trkeff;
  }
 
  int period_index = 0;
  period_index = GetPeriodIndex();

  int nbjet = NBJet(jetLooseColl, snu::KJet::CSVv2, snu::KJet::Medium, period_index);

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  snu::KParticle MET;
  METPt = CorrectedMETRochester(muonLooseColl, METPt, METPhi, true);
  METPhi = CorrectedMETRochester(muonLooseColl, METPt, METPhi, false);
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);

  double Z_mass = 91.1876;
  bool is_mue = ( ((muonLooseColl.size() == 1) && (electronLooseColl.size() == 1))
                 && ((muonTightColl.size() == 1) && (electronTightColl.size() == 1)) );
  bool is_mumue = ( ((muonLooseColl.size() == 2) && (electronLooseColl.size() == 1))
		 && ((muonTightColl.size() == 2) && (electronTightColl.size() == 1)) );
  bool is_mumuee = ( ((muonLooseColl.size() == 2) && (electronLooseColl.size() == 2))
		 && ((muonTightColl.size() == 2) && (electronTightColl.size() == 2)) );
  bool is_mumumumu = ( ((muonLooseColl.size() == 4) && (electronLooseColl.size() == 0))
		 && ((muonTightColl.size() == 4) && (electronTightColl.size() == 0)) );

  double weight_err = -999.;

  snu::KParticle lep[4];
  if((is_mue)){
    lep[0] = muonLooseColl.at(0);
    lep[1] = electronLooseColl.at(0);
    lep[2].SetPxPyPzE(999., 999., 999., 9999999.);
    lep[3].SetPxPyPzE(999., 999., 999., 9999999.);
  }

  if((is_mumue)){
    lep[0] = muonLooseColl.at(0);
    lep[1] = muonLooseColl.at(1);
    lep[2] = electronLooseColl.at(0);
    lep[3].SetPxPyPzE(999., 999., 999., 9999999.);
  }

  if((is_mumuee)){
    lep[0] = muonLooseColl.at(0);
    lep[1] = muonLooseColl.at(1);
    lep[2] = electronLooseColl.at(0);
    lep[3] = electronLooseColl.at(1);
  }
  
  if((is_mumumumu)){
    lep[0] = muonLooseColl.at(0);
    lep[1] = muonLooseColl.at(1);
    lep[2] = muonLooseColl.at(2);
    lep[3] = muonLooseColl.at(3);
  }
  
  if( is_mue && pass_muel_trig ){
    float weight_trigger = WeightByTrigger(triggerlist_mumu, TargetLumi);
    float this_weight, this_weight_err;
    this_weight = weight;
    if(!isData) this_weight*=weight_trigger;
    this_weight_err = 0.;

    bool is_preselection=false;
    bool is_DY=false;
    bool is_ttbar=false;

    bool p_lepton_OSlep = ((lep[0].Charge() != lep[1].Charge()));

    bool p_lepton_pt = ((lep[0].Pt() > 15) && (lep[1].Pt() > 25));

    bool p_MET_30 = (MET.Pt() > 30);
    bool p_MET_20 = (MET.Pt() > 20);

    bool p_bjet_N_0 = ( nbjet == 0);

    if( p_lepton_OSlep && p_lepton_pt ) is_preselection = true;

    if( is_preselection ){
      TString suffix = "preselection_mue";
      FillCLHist(sssf_mumue, suffix, eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, suffix+"_up", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight+this_weight_err);
      FillCLHist(sssf_mumue, suffix+"_down", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight-this_weight_err);
      FillUpDownHist("number_of_events_"+suffix, 0., this_weight, this_weight_err, 0., 1., 1);
      FillUpDownHist("PFMET_"+suffix, MET.Pt(), this_weight, this_weight_err, 0., 500., 500);
      FillUpDownHist("NJets_"+suffix, jetLooseColl.size(), this_weight, this_weight_err, 0., 5., 5);
      FillUpDownHist("dilepton_mass_"+suffix, (lep[0]+lep[1]).M(), this_weight, this_weight_err, 0., 250., 250);
    }

  }

  if( is_mumue && pass_mumu_trig ){

    float weight_trigger = WeightByTrigger(triggerlist_mumu, TargetLumi);
    float this_weight, this_weight_err;
    this_weight = weight;
    if(!isData) this_weight*=weight_trigger;
    this_weight_err = 0.;

    bool is_preselection=false;
    bool is_WZ=false;
    bool is_Zjet=false;

    bool p_lepton_OSlep = ((lep[0].Charge()) != (lep[1].Charge()));

    snu::KParticle Z_candidate;
    Z_candidate = lep[0] + lep[1];
    bool p_Z_candidate_mass = ((fabs(Z_candidate.M() - Z_mass) < 10));

    bool p_W_transverse_mass = (MT(lep[2],MET) > 30);

    bool p_Zlepton_pt = (((lep[0].Pt() > 20) && (lep[1].Pt() > 10)));
    bool p_Wlepton_pt = (lep[2].Pt() > 20);

    bool p_MET_30 = (MET.Pt() > 30);
    bool p_MET_20 = (MET.Pt() > 20);

    bool p_trilepton_mass = ( ((lep[0]+lep[1]+lep[2]).M() > 100) );

    bool p_bjet_N_0 = ( nbjet == 0 );

    if( Z_candidate.M() < 4 ) goto stop;

    if( p_lepton_OSlep && p_Zlepton_pt												) is_preselection= true;
    if( p_lepton_OSlep && p_Zlepton_pt && p_Wlepton_pt && p_Z_candidate_mass                        && p_MET_30 && p_trilepton_mass && p_bjet_N_0) is_WZ= true;
    if( p_lepton_OSlep && p_Zlepton_pt && p_Z_candidate_mass &&!p_W_transverse_mass &&!p_MET_20 && p_trilepton_mass && p_bjet_N_0) is_Zjet = true;

    if( is_preselection ){
      TString suffix = "preselection_mumue";
      FillCLHist(sssf_mumue, suffix, eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, suffix+"_up", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight+this_weight_err);
      FillCLHist(sssf_mumue, suffix+"_down", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight-this_weight_err);
      FillUpDownHist("number_of_events_"+suffix, 0., this_weight, this_weight_err, 0., 1., 1);
      FillUpDownHist("PFMET_"+suffix, MET.Pt(), this_weight, this_weight_err, 0., 500., 500);
      FillUpDownHist("NJets_"+suffix, jetLooseColl.size(), this_weight, this_weight_err, 0., 5., 5);
      FillUpDownHist("Z_candidate_mass_"+suffix, Z_candidate.M(), this_weight, this_weight_err, 0., 200., 200);
      FillUpDownHist("W_transverse_mass_"+suffix, MT(lep[2],MET), this_weight, this_weight_err, 0., 150., 150);
      FillUpDownHist("trilepton_mass_"+suffix, (lep[0]+lep[1]+lep[2]).M(), this_weight, this_weight_err, 0., 250., 250);
    }

    if( is_WZ ){
      TString suffix = "WZ_mumue";
      FillCLHist(sssf_mumue, suffix, eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, suffix+"_up", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight+this_weight_err);
      FillCLHist(sssf_mumue, suffix+"_down", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight-this_weight_err);
      FillUpDownHist("number_of_events_"+suffix, 0., this_weight, this_weight_err, 0., 1., 1);
      FillUpDownHist("PFMET_"+suffix, MET.Pt(), this_weight, this_weight_err, 0., 500., 500);
      FillUpDownHist("NJets_"+suffix, jetLooseColl.size(), this_weight, this_weight_err, 0., 5., 5);
      FillUpDownHist("Z_candidate_mass_"+suffix, Z_candidate.M(), this_weight, this_weight_err, 70., 110., 40); 
      FillUpDownHist("W_transverse_mass_"+suffix, MT(lep[2],MET), this_weight, this_weight_err, 0., 150., 150);
      FillUpDownHist("trilepton_mass_"+suffix, (lep[0]+lep[1]+lep[2]).M(), this_weight, this_weight_err, 0., 250., 250); 
    }

    if( is_Zjet ){
      TString suffix = "Zjet_mumue";
      FillCLHist(sssf_mumue, suffix, eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, suffix+"_up", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight+this_weight_err);
      FillCLHist(sssf_mumue, suffix+"_down", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight-this_weight_err);
      FillUpDownHist("number_of_events_"+suffix, 0., this_weight, this_weight_err, 0., 1., 1);
      FillUpDownHist("PFMET_"+suffix, MET.Pt(), this_weight, this_weight_err, 0., 500., 500);
      FillUpDownHist("NJets_"+suffix, jetLooseColl.size(), this_weight, this_weight_err, 0., 5., 5);
      FillUpDownHist("Z_candidate_mass_"+suffix, Z_candidate.M(), this_weight, this_weight_err, 70., 110., 40);
      FillUpDownHist("W_transverse_mass_"+suffix, MT(lep[2],MET), this_weight, this_weight_err, 0., 150., 150);
      FillUpDownHist("trilepton_mass_"+suffix, (lep[0]+lep[1]+lep[2]).M(), this_weight, this_weight_err, 0., 250., 250);
    }

  }



  if( is_mumuee && pass_mumu_trig ){

    float weight_trigger = WeightByTrigger(triggerlist_mumu, TargetLumi);
    float this_weight, this_weight_err;
    this_weight = weight;
    if(!isData) this_weight*=weight_trigger;
    this_weight_err = 0.;

    bool is_preselection=false;
    bool is_ZZ=false;

    bool p_lepton_OSlep = ((lep[0].Charge() != lep[1].Charge())
		        && (lep[2].Charge() != lep[3].Charge()));

    snu::KParticle tetralep = lep[0] + lep[1] + lep[2] + lep[3];

    snu::KParticle Z_candidate[2];
    Z_candidate[0] = lep[0] + lep[1];
    Z_candidate[1] = lep[2] + lep[3];
    bool p_Z_candidate_mass = ((fabs(Z_candidate[0].M() - Z_mass) < 10)
			    && (fabs(Z_candidate[1].M() - Z_mass) < 10));

    bool p_Zlepton_pt = ((lep[0].Pt() > 20 && lep[1].Pt() > 10)
		      && (lep[2].Pt() > 20 && lep[3].Pt() > 10));

    if( Z_candidate[0].M() < 4 || Z_candidate[1].M() < 4 ) goto stop; 

    if( p_lepton_OSlep && p_Zlepton_pt  		     ) is_preselection= true;
    if( p_lepton_OSlep && p_Zlepton_pt && p_Z_candidate_mass ) is_ZZ= true;

    if( is_preselection ){
      TString suffix = "preselection_mumuee";
      FillUpDownHist("h_leadingMuonPt_"+suffix, lep[0].Pt(), this_weight, this_weight_err, 0., 400., 400);
      FillUpDownHist("h_leadingMuonEta_"+suffix, lep[0].Eta(), this_weight, this_weight_err, -3., 3., 30);
      FillUpDownHist("h_secondMuonPt_" +suffix, lep[1].Pt(), this_weight, this_weight_err, 0., 400., 400);
      FillUpDownHist("h_secondMuonEta_" +suffix, lep[1].Eta(), this_weight, this_weight_err, -3., 3., 30);
      FillUpDownHist("h_leadingElectronPt_"+suffix, lep[2].Pt(), this_weight, this_weight_err, 0., 400., 400);
      FillUpDownHist("h_leadingElectronEta_"+suffix, lep[2].Eta(), this_weight, this_weight_err, -3., 3., 30);
      FillUpDownHist("h_secondElectronPt_" +suffix, lep[3].Pt(), this_weight, this_weight_err, 0., 400., 400);
      FillUpDownHist("h_secondElectronEta_" +suffix, lep[3].Eta(), this_weight, this_weight_err, -3., 3., 30);
      FillUpDownHist("number_of_events_"+suffix, 0., this_weight, this_weight_err, 0., 1., 1);
      FillUpDownHist("PFMET_"+suffix, MET.Pt(), this_weight, this_weight_err, 0., 500., 500);
      FillUpDownHist("NJets_"+suffix, jetLooseColl.size(), this_weight, this_weight_err, 0., 5., 5);
      FillUpDownHist("Z_mumu_candidate_mass_"+suffix, Z_candidate[0].M(), this_weight, this_weight_err, 0., 200., 200);
      FillUpDownHist("Z_elel_candidate_mass_"+suffix, Z_candidate[1].M(), this_weight, this_weight_err, 0., 200., 200);
      FillUpDownHist("tetralepton_mass_"+suffix, tetralep.M(), this_weight, this_weight_err, 0., 500., 500);
    }
    if( is_ZZ ){
      TString suffix = "ZZ_mumuee";
      FillUpDownHist("h_leadingMuonPt_"+suffix, lep[0].Pt(), this_weight, this_weight_err, 0., 400., 400);
      FillUpDownHist("h_leadingMuonEta_"+suffix, lep[0].Eta(), this_weight, this_weight_err, -3., 3., 30);
      FillUpDownHist("h_secondMuonPt_" +suffix, lep[1].Pt(), this_weight, this_weight_err, 0., 400., 400);
      FillUpDownHist("h_secondMuonEta_" +suffix, lep[1].Eta(), this_weight, this_weight_err, -3., 3., 30);
      FillUpDownHist("h_leadingElectronPt_"+suffix, lep[2].Pt(), this_weight, this_weight_err, 0., 400., 400);
      FillUpDownHist("h_leadingElectronEta_"+suffix, lep[2].Eta(), this_weight, this_weight_err, -3., 3., 30);
      FillUpDownHist("h_secondElectronPt_" +suffix, lep[3].Pt(), this_weight, this_weight_err, 0., 400., 400);
      FillUpDownHist("h_secondElectronEta_" +suffix, lep[3].Eta(), this_weight, this_weight_err, -3., 3., 30);
      FillUpDownHist("number_of_events_"+suffix, 0., this_weight, this_weight_err, 0., 1., 1);
      FillUpDownHist("PFMET_"+suffix, MET.Pt(), this_weight, this_weight_err, 0., 500., 500);
      FillUpDownHist("NJets_"+suffix, jetLooseColl.size(), this_weight, this_weight_err, 0., 5., 5);
      FillUpDownHist("Z_mumu_candidate_mass_"+suffix, Z_candidate[0].M(), this_weight, this_weight_err, 70., 110., 40);
      FillUpDownHist("Z_elel_candidate_mass_"+suffix, Z_candidate[1].M(), this_weight, this_weight_err, 70., 110., 40);
      FillUpDownHist("tetralepton_mass_"+suffix, tetralep.M(), this_weight, this_weight_err, 0., 500., 500);
    }
  }


  stop:

  return;

}// End of execute event loop
  


void HNSSSFMuMuE_CR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNSSSFMuMuE_CR::BeginCycle() throw( LQError ){
  
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

HNSSSFMuMuE_CR::~HNSSSFMuMuE_CR() {
  
  Message("In HNSSSFMuMuE_CR Destructor" , INFO);
  
}


void HNSSSFMuMuE_CR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNSSSFMuMuE_CR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNSSSFMuMuE_CRCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNSSSFMuMuE_CR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
 
  out_muons.clear();
  out_electrons.clear();
}


int HNSSSFMuMuE_CR::GetPeriodIndex(void){
  if( isData ){
    if( k_sample_name.Contains("B") ||  k_sample_name.Contains("C") || k_sample_name.Contains("D") || k_sample_name.Contains("E") || k_sample_name.Contains("F") ){
      return 1;
    }
    else return 7;
  }
  else return 0;
}

