// $Id: HNSSSFMuMuE_CR_FR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNSSSFMuMuE_CR_FR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNSSSFMuMuE_CR_FR.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"

//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNSSSFMuMuE_CR_FR);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNSSSFMuMuE_CR_FR::HNSSSFMuMuE_CR_FR() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNSSSFMuMuE_CR_FR");
  
  Message("In HNSSSFMuMuE_CR_FR constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNSSSFMuMuE_CR_FR::InitialiseAnalysis() throw( LQError ) {
  
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


void HNSSSFMuMuE_CR_FR::ExecuteEvents()throw( LQError ){
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  
  FillCutFlow("NoCut", weight);

  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);

  if(!PassMETFilter()) return;     /// Initial event cuts :
  FillCutFlow("EventCut", weight);

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex

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

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_HN", 30);

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TRI_TIGHT",false);
  std::vector<snu::KMuon> muonKeepfakeLooseColl = GetMuons("MUON_HN_TRI_LOOSE",true);
  std::vector<snu::KMuon> muonKeepfakeTightColl = GetMuons("MUON_HN_TRI_TIGHT",true);
  std::vector<snu::KMuon> muonNonpromptColl;
  muonNonpromptColl.clear();

  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,false,"ELECTRON_HN_TIGHT");
  std::vector<snu::KElectron> electronKeepfakeColl = GetElectrons(true,true,"ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronNonpromptColl;
  electronNonpromptColl.clear();

  bool DoMCClosure = false;//std::find(k_flags.begin(), k_flags.end(), "MCClosure") != k_flags.end();

  if( DoMCClosure ){
    std::vector<snu::KMuon> muonPromptColl = GetHNTriMuonsByLooseRelIso(0.4, false);
    if( !(muonPromptColl.size() == 2) ) muonNonpromptColl = GetHNTriMuonsByLooseRelIso(0.4, true);
  }
  if( DoMCClosure ){
    std::vector<snu::KElectron> electronPromptColl = GetHNElectronsByLooseRelIso(0.5, false);
    cout << electronPromptColl.size() <<endl;
    if( !(electronPromptColl.size() == 2) ) electronNonpromptColl = GetHNElectronsByLooseRelIso(0.5, true);
    cout<<electronNonpromptColl.size() << endl;
    cout<< "                             " << endl;
  }

  //CorrectMuonMomentum(muonLooseColl);
  float muon_trkeff = mcdata_correction->MuonTrackingEffScaleFactor(muonLooseColl);
  float electron_idsf = mcdata_correction->ElectronScaleFactor("ELECTRON_HN_FAKELOOSE", electronLooseColl);
  float electron_reco = mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);

  int period_index = 0;
  period_index = GetPeriodIndex();

  int nbjet = NBJet(jetLooseColl, snu::KJet::CSVv2, snu::KJet::Medium, period_index);
  double HT = 0.;
  for( int i=0; i<jetLooseColl.size(); i++){
    HT += jetLooseColl.at(i).Pt();
  }

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  snu::KParticle MET;
  METPt = CorrectedMETRochester(muonLooseColl, METPt, METPhi, true);
  METPhi = CorrectedMETRochester(muonLooseColl, METPt, METPhi, false);
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);

  double Z_mass = 91.1876;

  bool is_mumue = ( ((muonLooseColl.size() == 2) && (electronLooseColl.size() == 1)) && !((muonTightColl.size() == 2) && (electronTightColl.size() == 1)) );
  bool is_emu = ( ((muonLooseColl.size() == 1) && (electronLooseColl.size() == 1)) && !((muonTightColl.size() == 1) && (electronTightColl.size() == 1)) );

  if((is_mumue)){
    weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_TIGHT", 1);
    weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_TIGHT", 1);
  }
  else if((is_emu)){
    weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 1, electronLooseColl, "ELECTRON_HN_TIGHT", 1);
    weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 1, electronLooseColl, "ELECTRON_HN_TIGHT", 1);
  }
  // =============================================================================
  // == MuMuE Selection===========================================================
  // ====================================

  if( is_mumue && pass_mumu_trig ){

    SF[0] = muonLooseColl.at(0);
    SF[1] = muonLooseColl.at(1);
    OF = electronLooseColl.at(0);

    bool is_mumu_trig_WZ=false; bool is_mumu_trig_Zjet=false;

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

    if( p_lepton_OSSF && p_Z_candidate_mass                        && p_lepton_pt && p_MET_30 && p_trilepton_mass ) CR_WZ_mumue = true;
    if( p_lepton_OSSF && p_Z_candidate_mass &&!p_W_transverse_mass && p_lepton_pt &&!p_MET_20 && p_trilepton_mass ) CR_Zjet_mumue = true;

    // ========== WZ selection =====================================
    if( is_mumu_trig_WZ ){
      DrawHistograms("mumu_trig_WZ", SF, OF, MET, jetLooseColl, weight);
      FillCLHist(sssf_mumue, "mumu_trig_WZ", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, weight);
      FillUpDownHist("NBJets_mumu_trig_WZ", nbjet, weight, weight_err, 0., 5., 5);

    }

    // ========== Z+lepton selection =====================================
    if( is_mumu_trig_Zjet ){
      DrawHistograms("mumu_trig_Zjet", SF, OF, MET, jetLooseColl, weight);
      FillCLHist(sssf_mumue, "mumu_trig_Zjet", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, weight);
      FillUpDownHist("NBJets_mumu_trig_Zjet", nbjet, weight, weight_err, 0., 5., 5);
    }

  }



  if( is_mumue && pass_emu_trig ){

    SF[0] = muonLooseColl.at(0);
    SF[1] = muonLooseColl.at(1);
    OF = electronLooseColl.at(0);

    bool is_emu_trig_ttbarjet=false;

    // == OS muon pair cut
    bool p_lepton_OSOF = ((SF[0].Charge()) != (OF.Charge()));
    bool p_lepton_SSSF = ((SF[0].Charge()) == (SF[1].Charge()));
    // == lepton Pt cut
    bool p_lepton_pt = (((SF[0].Pt() > 15) && (SF[1].Pt() > 15) && (OF.Pt() > 25)));
    // == MET cut
    bool p_MET_30 = (MET.Pt() > 30);
    // == jet cut
    bool p_jet_N_2 = ( (jetLooseColl.size() > 1) );
    // == bjet cut
    bool p_bjet_N_1 = ( nbjet > 0 );
      
    if( (p_lepton_OSOF && p_lepton_SSSF) && p_lepton_pt && p_jet_N_2 && p_bjet_N_1 ) CR_ttbarjet_mumue = true;

    // ========== ttbar+lep selection =====================================
    if( is_emu_trig_ttbarjet ){
      DrawHistograms("emu_trig_ttbarjet", SF, OF, MET, jetLooseColl, weight);
      FillCLHist(sssf_mumue, "emu_trig_ttbarjet", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, weight);
      FillUpDownHist("NBJets_emu_trig_ttbarjet", nbjet, weight, weight_err, 0., 5., 5);
    }

  }

          
  // =============================================================================
  // == EMu Selection============================================================
  // ====================================

  if( is_emu && pass_emu_trig ){

    SF[0] = muonLooseColl.at(0);
    OF = electronLooseColl.at(0);
    SF[1].SetPxPyPzE(0,0,0,0);

    bool is_emu_trig_ttbar=false;

    // == OS muon pair cut
    bool p_lepton_OSOF = ((SF[0].Charge()) != (OF.Charge()));
    // == lepton Pt cut
    bool p_lepton_pt = (((SF[0].Pt() > 35) && (OF.Pt() > 35)));
    // == MET cut
    bool p_MET_40 = (MET.Pt() > 40);
    // == jet cut
    bool p_jet_N_2 = ( (jetLooseColl.size() > 1) );
    // == bjet cut
    bool p_bjet_N_2 = ( nbjet > 1 );
    // == dilepton mass cut
    bool p_dilepton_mass = true;//( ((SF[0]+OF).M()) > 20 ) ;

    if( p_lepton_OSOF && p_lepton_pt && p_jet_N_2 && p_bjet_N_2 && p_MET_40 ) CR_ttbar_emu = true;

    if( p_lepton_OSOF && p_lepton_pt ){
      FillHist("dilepton_mass_emu", ((SF[0]+OF).M()), weight, 0., 500., 100);
      DrawHistograms("emu", SF, OF, MET, jetLooseColl, weight);
      FillCLHist(sssf_mumue, "emu", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, weight);
      FillHist("NBJets_emu", nbjet, weight, 0., 5., 5);
    }
   
    if( CR_ttbar_emu ){
      FillHist("dilepton_mass_emu_ttbar", ((SF[0]+OF).M()), weight, 0., 500., 100);
      DrawHistograms("emu_ttbar", SF, OF, MET, jetLooseColl, weight);
      FillCLHist(sssf_mumue, "emu_ttbar", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, weight);
      FillHist("NBJets_emu", nbjet, weight, 0., 5., 5);
    }

  }


 



  return;

}// End of execute event loop
  


void HNSSSFMuMuE_CR_FR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNSSSFMuMuE_CR_FR::BeginCycle() throw( LQError ){
  
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

HNSSSFMuMuE_CR_FR::~HNSSSFMuMuE_CR_FR() {
  
  Message("In HNSSSFMuMuE_CR_FR Destructor" , INFO);
  
}


void HNSSSFMuMuE_CR_FR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNSSSFMuMuE_CR_FR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNSSSFMuMuE_CR_FRCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNSSSFMuMuE_CR_FR::DrawHistograms(TString suffix, snu::KParticle SF[], snu::KParticle OF, snu::KParticle MET,  std::vector<snu::KJet> jetLooseColl, double weight){

  if( suffix == "2mu" || suffix == "2mu_ttW" ){
    FillUpDownHist("number_of_events_"+suffix, 0., weight, weight_err, 0., 1., 1);

    FillUpDownHist("PFMET_"+suffix, MET.Pt(), weight, weight_err, 0., 500., 500);
    FillUpDownHist("NJets_"+suffix, jetLooseColl.size(), weight, weight_err, 0., 5., 5);

    FillUpDownHist("W_transverse_mass_"+suffix, MT(OF,MET), weight, weight_err, 0., 500., 500);
    FillUpDownHist("Z_candidate_mass_"+suffix, (SF[0]+SF[1]).M(), weight, weight_err, 0., 500., 500);

    FillUpDownHist("SFleadingLeptonPt_"+suffix, SF[0].Pt(), weight, weight_err, 0., 500., 500);
    FillUpDownHist("SFleadingLeptonEta_"+suffix, SF[0].Eta(), weight, weight_err, -3., 3., 60);
    FillUpDownHist("SFsecondLeptonPt_"+suffix, SF[1].Pt(), weight, weight_err, 0., 500., 500);
    FillUpDownHist("SFsecondLeptonEta_"+suffix, SF[1].Eta(), weight, weight_err, -3., 3., 60);

  }

  else{
    FillUpDownHist("number_of_events_"+suffix, 0., weight, weight_err, 0., 1., 1);

    FillUpDownHist("PFMET_"+suffix, MET.Pt(), weight, weight_err, 0., 500., 500);
    FillUpDownHist("NJets_"+suffix, jetLooseColl.size(), weight, weight_err, 0., 5., 5);

    FillUpDownHist("W_transverse_mass_"+suffix, MT(OF,MET), weight, weight_err, 0., 500., 500);
    FillUpDownHist("Z_candidate_mass_"+suffix, (SF[0]+SF[1]).M(), weight, weight_err, 0., 500., 500);

    FillUpDownHist("SFleadingLeptonPt_"+suffix, SF[0].Pt(), weight, weight_err, 0., 500., 500);
    FillUpDownHist("SFleadingLeptonEta_"+suffix, SF[0].Eta(), weight, weight_err, -3., 3., 60);
    FillUpDownHist("SFsecondLeptonPt_"+suffix, SF[1].Pt(), weight, weight_err, 0., 500., 500);
    FillUpDownHist("SFsecondLeptonEta_"+suffix, SF[1].Eta(), weight, weight_err, -3., 3., 60);

    FillUpDownHist("OFLeptonPt_"+suffix, OF.Pt(), weight, weight_err, 0., 500., 500);
    FillUpDownHist("OFLeptonEta_"+suffix, OF.Eta(), weight, weight_err, -3., 3., 60);

  }

  return;

}


void HNSSSFMuMuE_CR_FR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
 
  SF[0].SetPxPyPzE(0,0,0,0); SF[1].SetPxPyPzE(0,0,0,0); OF.SetPxPyPzE(0,0,0,0);
  CR_WZ_mumue=false; CR_Zjet_mumue=false; CR_ttbar_emu=false; CR_ttbarjet_mumue=false;
  out_muons.clear();
  out_electrons.clear();
}


int HNSSSFMuMuE_CR_FR::GetPeriodIndex(void){
  if( isData ){
    if( k_sample_name.Contains("B") ||  k_sample_name.Contains("C") || k_sample_name.Contains("D") || k_sample_name.Contains("E") || k_sample_name.Contains("F") ){
      return 1;
    }
    else return 7;
  }
  else return 0;
}
