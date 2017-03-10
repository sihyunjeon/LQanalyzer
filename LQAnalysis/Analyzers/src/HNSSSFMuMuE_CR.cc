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

  std::vector<TString> triggerlist_emu, triggerlist_mumu;
  triggerlist_emu.clear();
  triggerlist_mumu.clear();

  bool pass_emu_trig = false;
  bool pass_mumu_trig = false;

  triggerlist_emu.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_emu.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_emu.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if(PassTriggerOR(triggerlist_emu)) pass_emu_trig = true;

  triggerlist_mumu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_mumu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  if(PassTriggerOR(triggerlist_mumu)) pass_mumu_trig = true;

  if(!pass_emu_trig && !pass_mumu_trig) return;

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_HN", 30., 2.4);

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TRI_TIGHT",false);
  std::vector<snu::KMuon> muonKeepfakeLooseColl = GetMuons("MUON_HN_TRI_LOOSE",true);
  std::vector<snu::KMuon> muonKeepfakeTightColl = GetMuons("MUON_HN_TRI_TIGHT",true);
  std::vector<snu::KMuon> muonNonpromptColl;
  muonNonpromptColl.clear();

  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false,false,"ELECTRON16_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,false,"ELECTRON16_HN_TIGHT");
  std::vector<snu::KElectron> electronKeepfakeColl = GetElectrons(true,true,"ELECTRON16_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronNonpromptColl;
  electronNonpromptColl.clear();

  bool DoMCClosure = true;//std::find(k_flags.begin(), k_flags.end(), "MCClosure") != k_flags.end();

  if( DoMCClosure ){
    std::vector<snu::KMuon> muonPromptColl = GetHNTriMuonsByLooseRelIso(0.4, false);
    if( !(muonPromptColl.size() == 2) ) muonNonpromptColl = GetHNTriMuonsByLooseRelIso(0.4, true);
  }
  if( DoMCClosure ){
/*    std::vector<snu::KElectron> electronPromptColl = GetHNElectronsByLooseRelIso(0.5, false);
    cout << electronPromptColl.size() <<endl;
    if( !(electronPromptColl.size() == 2) ) electronNonpromptColl = GetHNElectronsByLooseRelIso(0.5, true);*/
  }

  CorrectMuonMomentum(muonLooseColl);
  float muon_trkeff = mcdata_correction->MuonTrackingEffScaleFactor(muonLooseColl);
  float electron_idsf = mcdata_correction->ElectronScaleFactor("ELECTRON16_HN_FAKELOOSE", electronLooseColl);
  float electron_reco = mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);

  double ev_weight = weight;
  if(!isData){
    weight *= pileup_reweight;
    weight *= electron_idsf;
    weight *= electron_reco;
    weight *= muon_trkeff;
    //ev_weight = w * trigger_sf * id_iso_sf *  pu_reweight*trigger_ps;
  }
 
  int period_index = 0;
  period_index = GetPeriodIndex();

  int nbjet = NBJet(jetLooseColl, snu::KJet::CSVv2, snu::KJet::Medium, period_index);

  HT = 0.;
  for( int i=0; i<jetLooseColl.size(); i++){
    HT += jetLooseColl.at(i).Pt();
  }

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  snu::KParticle MET;
  METPt = CorrectedMETRochester(muonLooseColl, METPt, METPhi, true);
  METPhi = CorrectedMETRochester(muonLooseColl, METPt, METPhi, false);
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);

  ST = eventbase->GetEvent().SumET();

  looseLT = 0.;
  for( int i=0; i<muonLooseColl.size(); i++){
    looseLT += muonLooseColl.at(i).Pt();
  }
  for( int i=0; i<electronLooseColl.size(); i++){
    looseLT += electronLooseColl.at(i).Pt();
  }
  tightLT = 0.;
  for( int i=0; i<muonTightColl.size(); i++){
    tightLT += muonTightColl.at(i).Pt();
  }
  for( int i=0; i<electronTightColl.size(); i++){
    tightLT += electronTightColl.at(i).Pt();
  }


  float STarray [] = {0., 50., 100., 150., 200., 250., 300., 500., 700., 1000.};
  float looseLTarray [] = {0., 50., 100., 150., 200., 250., 300., 500., 700., 1000.};
  float tightLTarray [] = {0., 50., 100., 150., 200., 250., 300., 500., 700., 1000.};

  double Z_mass = 91.1876;

  bool is_mumue = ( ((muonLooseColl.size() == 2) && (electronLooseColl.size() == 1)) && ((muonTightColl.size() == 2) && (electronTightColl.size() == 1)) );
  bool is_emu = ( ((muonLooseColl.size() == 1) && (electronLooseColl.size() == 1)) && ((muonTightColl.size() == 1) && (electronTightColl.size() == 1)) );
  
  // =============================================================================
  // == MuMuE Selection===========================================================
  // ====================================

  if( is_mumue && pass_mumu_trig ){

    float weight_trigger = WeightByTrigger(triggerlist_mumu, TargetLumi);
    float this_weight;
    if( isData ) this_weight = weight;
    if( !isData ) this_weight = weight * weight_trigger;

    SF[0] = muonLooseColl.at(0);
    SF[1] = muonLooseColl.at(1);
    OF = electronLooseColl.at(0);

    bool is_mumu_trig_WZ_mumue=false;
    bool is_mumu_trig_Zjet_mumue=false;

    // == OS muon pair cut
    bool p_lepton_OSSF = ((SF[0].Charge()) != (SF[1].Charge()));
    // == Z candidate mass cut
    snu::KParticle Z_candidate;
    Z_candidate = SF[0] + SF[1];
    bool p_Z_candidate_mass = ((fabs(Z_candidate.M() - Z_mass) < 10));
    // == W transverse mass cut
    bool p_W_transverse_mass = (MT(OF,MET) > 30);
    // == lepton Pt cut
    bool p_lepton_pt = (((SF[0].Pt() > 20) && (SF[1].Pt() > 10)) && (OF.Pt() > 10));
    // == MET cut
    bool p_MET_30 = (MET.Pt() > 30);
    bool p_MET_20 = (MET.Pt() > 20);
    // == trilepton mass cut
    bool p_trilepton_mass = ( ((SF[0]+SF[1]+OF).M() > 100) );
    // == jet cut
    bool p_jet_N_2 = ( (jetLooseColl.size() > 1) );
    // == bjet cut
    bool p_bjet_N_0 = ( nbjet == 0 );

    if( p_lepton_OSSF && p_Z_candidate_mass                        && p_lepton_pt && p_MET_30 && p_trilepton_mass ) is_mumu_trig_WZ_mumue = true;
    if( p_lepton_OSSF && p_Z_candidate_mass &&!p_W_transverse_mass && p_lepton_pt &&!p_MET_20 && p_trilepton_mass ) is_mumu_trig_Zjet_mumue = true;

    if( p_lepton_OSSF && p_lepton_pt ){
      DrawHistograms("mumu_trig_mumue", SF, OF, MET, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, "mumu_trig_mumue", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillHist("NBJets_mumu_trig_mumue", nbjet, this_weight, 0., 5., 5);
    }

    // ========== WZ selection =====================================
    if( is_mumu_trig_WZ_mumue ){
      DrawHistograms("mumu_trig_WZ_mumue", SF, OF, MET, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, "mumu_trig_WZ_mumue", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillHist("NBJets_mumu_trig_WZ_mumue", nbjet, this_weight, 0., 5., 5);
    }

    // ========== Z+lepton selection =====================================
    if( is_mumu_trig_Zjet_mumue ){
      DrawHistograms("mumu_trig_Zjet_mumue", SF, OF, MET, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, "mumu_trig_Zjet_mumue", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillHist("NBJets_mumu_trig_Zjet_mumue", nbjet, this_weight, 0., 5., 5);
    }

  }



  if( is_mumue && pass_emu_trig ){

    float weight_trigger = WeightByTrigger(triggerlist_emu, TargetLumi);
    float this_weight;
    if( isData ) this_weight = weight;
    if( !isData ) this_weight = weight * weight_trigger;

    SF[0] = muonLooseColl.at(0);
    SF[1] = muonLooseColl.at(1);
    OF = electronLooseColl.at(0);

    bool is_emu_trig_ttbarjet_mumue=false;

    // == OS muon pair cut
    bool p_lepton_OSOF = ((SF[0].Charge()) != (OF.Charge()));
    bool p_lepton_SSSF = ((SF[0].Charge()) == (SF[1].Charge()));
    // == lepton Pt cut
    bool p_lepton_pt = (((SF[0].Pt() > 20) && (SF[1].Pt() > 10) && (OF.Pt() > 20)));
    // == MET cut
    bool p_MET_20 = (MET.Pt() > 20);
    // == jet cut
    bool p_jet_N_2 = ( (jetLooseColl.size() > 1) );
    // == bjet cut
    bool p_bjet_N_1 = ( nbjet > 0 );
    // == dilepton mass cut
    bool p_dilepton_mass = true;//( ((SF[0]+OF).M()) > 20 ) ;
      
    if( (p_lepton_OSOF && p_lepton_SSSF) && p_lepton_pt && p_jet_N_2 && p_bjet_N_1 && p_MET_20 ) is_emu_trig_ttbarjet_mumue = true;

    if( (p_lepton_OSOF && p_lepton_SSSF) && p_lepton_pt ){
      FillHist("dilepton_mass_emu_trig_mumue", ((SF[0]+OF).M()), this_weight, 0., 500., 100);
      DrawHistograms("emu_trig_mumue", SF, OF, MET, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, "emu_trig_mumue", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillHist("NBJets_emu_trig_mumue", nbjet, this_weight, 0., 5., 5);
      FillHist("looseLT_emu_trig_mumue", looseLT, this_weight, 0., 1000., 100);
      FillHist("tightLT_emu_trig_mumue", tightLT, this_weight, 0., 1000., 100);
      FillHist("ST_emu_trig_mumue", ST, this_weight, 0., 3000., 150);
      FillHist("HT_emu_trig_mumue", HT, this_weight, 0., 1000., 100);
    }

    if( is_emu_trig_ttbarjet_mumue ){
      FillHist("dilepton_mass_emu_trig_ttbarjet_mumue", ((SF[0]+OF).M()), this_weight, 0., 500., 100);
      DrawHistograms("emu_trig_ttbarjet_mumue", SF, OF, MET, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, "emu_trig_ttbarjet_mumue", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillHist("NBJets_emu_trig_ttbarjet_mumue", nbjet, this_weight, 0., 5., 5);
      FillHist("looseLT_emu_trig_ttbarjet_mumue", looseLT, this_weight, 0., 1000., 100);
      FillHist("tightLT_emu_trig_ttbarjet_mumue", tightLT, this_weight, 0., 1000., 100);
      FillHist("ST_emu_trig_ttbarjet_mumue", ST, this_weight, 0., 3000., 150);
      FillHist("HT_emu_trig_ttbarjet_mumue", HT, this_weight, 0., 1000., 100);
    }

    if( (p_lepton_OSOF && p_lepton_SSSF) && p_lepton_pt && (ST > 300) ){
      FillHist("dilepton_mass_emu_trig_ttbarjet_mumue_ST", ((SF[0]+OF).M()), this_weight, 0., 500., 100);
      DrawHistograms("emu_trig_ttbarjet_mumue_ST", SF, OF, MET, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, "emu_trig_ttbarjet_mumue_ST", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillHist("NBJets_emu_trig_ttbarjet_mumue_ST", nbjet, this_weight, 0., 5., 5);
      FillHist("looseLT_emu_trig_ttbarjet_mumue_ST", looseLT, this_weight, 0., 1000., 100);
      FillHist("tightLT_emu_trig_ttbarjet_mumue_ST", tightLT, this_weight, 0., 1000., 100);
      FillHist("ST_emu_trig_ttbarjet_mumue_ST", ST, this_weight, 0., 3000., 150);
      FillHist("HT_emu_trig_ttbarjet_mumue_ST", HT, this_weight, 0., 1000., 100);

    }
  }


          
  // =============================================================================
  // == EMu Selection============================================================
  // ====================================

  if( is_emu && pass_emu_trig ){

    float weight_trigger = WeightByTrigger(triggerlist_emu, TargetLumi);
    float this_weight;
    if( isData ) this_weight = weight;
    if( !isData ) this_weight = weight * weight_trigger;

    SF[0] = muonLooseColl.at(0);
    OF = electronLooseColl.at(0);
    SF[1].SetPxPyPzE(0,0,0,0);

    snu::KMuon mu = muonLooseColl.at(0);
    snu::KElectron el = electronLooseColl.at(0);
    bool tautau =false;
    if(!isData && ( mu.MCFromTau() && el.MCFromTau()) ) tautau=true;

    bool is_emu_trig_ttbarjet_emu=false;

    // == OS muon pair cut
    bool p_lepton_OSOF = ((SF[0].Charge()) != (OF.Charge()));
    // == lepton Pt cut
    bool p_lepton_pt = (((SF[0].Pt() > 35) && (OF.Pt() > 35)));
    // == MET cut
    bool p_MET_30 = (MET.Pt() > 40);
    // == jet cut
    bool p_jet_N_2 = ( (jetLooseColl.size() > 1) );
    // == bjet cut
    bool p_bjet_N_2 = ( nbjet > 1 );
    // == dilepton mass cut
    bool p_dilepton_mass = true;//( ((SF[0]+OF).M()) > 20 ) ;

    if( p_lepton_OSOF && p_lepton_pt && p_jet_N_2 && p_bjet_N_2 && p_MET_30 ) is_emu_trig_ttbarjet_emu = true;

    if( p_lepton_OSOF && p_lepton_pt ){
      FillHist("dilepton_mass_emu_trig_emu", ((SF[0]+OF).M()), this_weight, 0., 500., 100);
      if(tautau){
	FillHist("dilepton_mass_tautau_emu_trig_emu", ((SF[0]+OF).M()), this_weight, 0., 500., 100);
        FillHist("PFMET_tautau_emu_trig_emu", MET.Pt(), this_weight, 0., 500., 100);
      } 
      DrawHistograms("emu_trig_emu", SF, OF, MET, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, "emu_trig_emu", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillHist("NBJets_emu_trig_emu", nbjet, this_weight, 0., 5., 5);
      FillHist("looseLT_emu_trig_emu", looseLT, this_weight, 0., 1000., 100);
      FillHist("tightLT_emu_trig_emu", tightLT, this_weight, 0., 1000., 100);
      FillHist("ST_emu_trig_emu", ST, this_weight, 0., 2000., 100);
      FillHist("HT_emu_trig_emu", HT, this_weight, 0., 1000., 100);

    }
   
    if( is_emu_trig_ttbarjet_emu ){
      FillHist("dilepton_mass_emu_trig_ttbarjet_emu", ((SF[0]+OF).M()), this_weight, 0., 500., 100);
      if(tautau){
        FillHist("dilepton_mass_tautau_emu_trig_ttbarjet_emu", ((SF[0]+OF).M()), this_weight, 0., 500., 100);
        FillHist("PFMET_tautau_emu_trig_ttbarjet_emu", MET.Pt(), this_weight, 0., 500., 100);
      }
      DrawHistograms("emu_trig_ttbarjet_emu", SF, OF, MET, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, "emu_trig_ttbarjet_emu", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillHist("NBJets_emu_trig_ttbarjet_emu", nbjet, this_weight, 0., 5., 5);
      FillHist("looseLT_emu_trig_ttbarjet_emu", looseLT, this_weight, 0., 1000., 100);
      FillHist("tightLT_emu_trig_ttbarjet_emu", tightLT, this_weight, 0., 1000., 100);
      FillHist("ST_emu_trig_ttbarjet_emu", ST, this_weight, 0., 2000., 100);
      FillHist("HT_emu_trig_ttbarjet_emu", HT, this_weight, 0., 1000., 100);
    }

    if( p_lepton_OSOF && p_lepton_pt && (ST > 300) ){
      FillHist("dilepton_mass_emu_trig_ttbarjet_emu_ST", ((SF[0]+OF).M()), this_weight, 0., 500., 100);
      if(tautau){
        FillHist("dilepton_mass_tautau_emu_trig_ttbarjet_emu_ST", ((SF[0]+OF).M()), this_weight, 0., 500., 100);
        FillHist("PFMET_tautau_emu_trig_ttbarjet_emu_ST", MET.Pt(), this_weight, 0., 500., 100);
      }
      DrawHistograms("emu_trig_ttbarjet_emu_ST", SF, OF, MET, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, "emu_trig_ttbarjet_emu_ST", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillHist("NBJets_emu_trig_ttbarjet_emu_ST", nbjet, this_weight, 0., 5., 5);
      FillHist("looseLT_emu_trig_ttbarjet_emu_ST", looseLT, this_weight, 0., 1000., 100);
      FillHist("tightLT_emu_trig_ttbarjet_emu_ST", tightLT, this_weight, 0., 1000., 100);
      FillHist("ST_emu_trig_ttbarjet_emu_ST", ST, this_weight, 0., 2000., 100);
      FillHist("HT_emu_trig_ttbarjet_emu_ST", HT, this_weight, 0., 1000., 100);
    }

  }


 



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


void HNSSSFMuMuE_CR::DrawHistograms(TString suffix, snu::KParticle SF[], snu::KParticle OF, snu::KParticle MET,  std::vector<snu::KJet> jetLooseColl, double weight){

  if( suffix == "2mu" || suffix == "2mu_ttW" ){
    FillHist("number_of_events_"+suffix, 0., weight, 0., 1., 1);

    FillHist("PFMET_"+suffix, MET.Pt(), weight, 0., 500., 500);
    FillHist("NJets_"+suffix, jetLooseColl.size(), weight, 0., 5., 5);

    FillHist("W_transverse_mass_"+suffix, MT(OF,MET), weight, 0., 500., 500);
    FillHist("Z_candidate_mass_"+suffix, (SF[0]+SF[1]).M(), weight, 0., 500., 500);

    FillHist("SFleadingLeptonPt_"+suffix, SF[0].Pt(), weight, 0., 500., 500);
    FillHist("SFleadingLeptonEta_"+suffix, SF[0].Eta(), weight, -3., 3., 60);
    FillHist("SFsecondLeptonPt_"+suffix, SF[1].Pt(), weight, 0., 500., 500);
    FillHist("SFsecondLeptonEta_"+suffix, SF[1].Eta(), weight, -3., 3., 60);

  } 

  else{
    FillHist("number_of_events_"+suffix, 0., weight, 0., 1., 1);

    FillHist("PFMET_"+suffix, MET.Pt(), weight, 0., 500., 500);
    FillHist("NJets_"+suffix, jetLooseColl.size(), weight, 0., 5., 5);

    FillHist("W_transverse_mass_"+suffix, MT(OF,MET), weight, 0., 500., 500);
    FillHist("Z_candidate_mass_"+suffix, (SF[0]+SF[1]).M(), weight, 0., 500., 500);

    FillHist("SFleadingLeptonPt_"+suffix, SF[0].Pt(), weight, 0., 500., 500);
    FillHist("SFleadingLeptonEta_"+suffix, SF[0].Eta(), weight, -3., 3., 60);
    FillHist("SFsecondLeptonPt_"+suffix, SF[1].Pt(), weight, 0., 500., 500);
    FillHist("SFsecondLeptonEta_"+suffix, SF[1].Eta(), weight, -3., 3., 60);

    FillHist("OFLeptonPt_"+suffix, OF.Pt(), weight, 0., 500., 500);
    FillHist("OFLeptonEta_"+suffix, OF.Eta(), weight, -3., 3., 60);

  }

  return;

}


void HNSSSFMuMuE_CR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
 
  SF[0].SetPxPyPzE(0,0,0,0); SF[1].SetPxPyPzE(0,0,0,0); OF.SetPxPyPzE(0,0,0,0);
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
