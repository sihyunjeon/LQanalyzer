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

  if ( k_sample_name.Contains("DYJets") || k_sample_name.Contains("WJets") || k_sample_name.Contains("TTJets") ){
    MCClosure();
    return;
  }



  std::vector<TString> triggerlist_mumu;
  triggerlist_mumu.clear();

  bool pass_mumu_trig = false;

  triggerlist_mumu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_mumu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  if(PassTriggerOR(triggerlist_mumu)) pass_mumu_trig = true;
  if( !pass_mumu_trig ) return;

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_HN", 30., 2.4);

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TRI_TIGHT",false);

  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false,false,"ELECTRON_HN_LOWDXY_FAKELOOSE");
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,false,"ELECTRON_HN_LOWDXY_TIGHT");

  int period_index = 0;
  period_index = GetPeriodIndex();

  int nbjet = NBJet(jetLooseColl, snu::KJet::CSVv2, snu::KJet::Medium, period_index);

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  snu::KParticle MET;
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);

  double Z_mass = 91.1876;

  bool is_mumue = ( ((muonLooseColl.size() == 2) && (electronLooseColl.size() == 1)) );

  double weight_err = -999.;
  if((is_mumue)){
    weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", 1);
    weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", 1);

    lep[0] = muonLooseColl.at(0);
    lep[1] = muonLooseColl.at(1);
    lep[2] = electronLooseColl.at(0);
  }
  
  // =============================================================================
  // == MuMuE Selection===========================================================
  // ====================================

  if( is_mumue && pass_mumu_trig ){

    float this_weight, this_weight_err;
    this_weight = weight;
    this_weight_err = this_weight_err;

    lep[0] = muonLooseColl.at(0);
    lep[1] = muonLooseColl.at(1);
    lep[2] = electronLooseColl.at(0);

    bool is_preselection=false;
    bool is_WZ=false;
    bool is_Zjet=false;

    bool p_lepton_OSlep = ((lep[0].Charge()) != (lep[1].Charge()));

    snu::KParticle Z_candidate;
    Z_candidate = lep[0] + lep[1];
    bool p_Z_candidate_mass = ((fabs(Z_candidate.M() - Z_mass) < 10));

    bool p_W_transverse_mass = (MT(lep[2],MET) > 30);

    bool p_lepton_pt = (((lep[0].Pt() > 20) && (lep[1].Pt() > 10)) && (lep[2].Pt() > 10));

    bool p_MET_30 = (MET.Pt() > 30);
    bool p_MET_20 = (MET.Pt() > 20);

    bool p_trilepton_mass = ( ((lep[0]+lep[1]+lep[2]).M() > 100) );

    bool p_bjet_N_0 = ( nbjet == 0 );

    if( p_lepton_OSlep && p_lepton_pt                                                                                           ) is_preselection= true;
    if( p_lepton_OSlep && p_Z_candidate_mass                        && p_lepton_pt && p_MET_30 && p_trilepton_mass && p_bjet_N_0) is_WZ= true;
    if( p_lepton_OSlep && p_Z_candidate_mass &&!p_W_transverse_mass && p_lepton_pt &&!p_MET_20 && p_trilepton_mass && p_bjet_N_0) is_Zjet = true;

    if( is_preselection ){
      TString suffix = "preselection_mumue";
      FillCLHist(sssf_mumue, suffix, eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight);
      FillCLHist(sssf_mumue, suffix+"_up", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight+this_weight_err);
      FillCLHist(sssf_mumue, suffix+"_down", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, this_weight-this_weight_err);
      FillUpDownHist("number_of_events_"+suffix, 0., this_weight, this_weight_err, 0., 1., 1);
      FillUpDownHist("PFMET_"+suffix, MET.Pt(), this_weight, this_weight_err, 0., 500., 500);
      FillUpDownHist("NJets_"+suffix, jetLooseColl.size(), this_weight, this_weight_err, 0., 5., 5);
      FillUpDownHist("Z_candidate_mass_"+suffix, Z_candidate.M(), this_weight, this_weight_err, 0., 120., 120);
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
      FillUpDownHist("Z_candidate_mass_"+suffix, Z_candidate.M(), this_weight, this_weight_err, 0., 120., 120); 
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
      FillUpDownHist("Z_candidate_mass_"+suffix, Z_candidate.M(), this_weight, this_weight_err, 0., 120., 120);
      FillUpDownHist("W_transverse_mass_"+suffix, MT(lep[2],MET), this_weight, this_weight_err, 0., 150., 150);
      FillUpDownHist("trilepton_mass_"+suffix, (lep[0]+lep[1]+lep[2]).M(), this_weight, this_weight_err, 0., 250., 250);
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


void HNSSSFMuMuE_CR_FR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
 
  lep[0].SetPxPyPzE(0,0,0,0); lep[1].SetPxPyPzE(0,0,0,0); lep[2].SetPxPyPzE(0,0,0,0);
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


void HNSSSFMuMuE_CR_FR::MCClosure(void){

  double this_weight = 1;

  this_weight = 1 *MCweight;

  bool pass_mm_trig = false;
  if(PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")) pass_mm_trig = true;
  if( !pass_mm_trig ) return;

  std::vector<snu::KMuon> muonPRColl = GetHNTriMuonsByLooseRelIso(0.4, false);
  std::vector<snu::KMuon> muonBaseColl = GetHNTriMuonsByLooseRelIso(0.4, true);
  std::vector<snu::KMuon> muonPredColl;
  std::vector<snu::KMuon> muonObsColl;
  muonPredColl.clear();
  muonObsColl.clear();

  std::vector<snu::KElectron> electronLooseColl = GetElectrons(true,true,"ELECTRON_HN_LOWDXY_FAKELOOSE");
  if( electronLooseColl.size() > 0 ) return;

  if( (muonPRColl.size() == 2) ){
    if( (muonPRColl.at(0).Charge() == muonPRColl.at(1).Charge()) ) FillHist("PR_SS", 0., 1., 0., 1., 1);
  }

  if( (muonPRColl.size() < 2) ){
    if( (muonBaseColl.size() == 2) ){
      if( (muonBaseColl.at(0).Charge() == muonBaseColl.at(1).Charge()) ){
        if( (muonBaseColl.at(0).Pt() > 20) && (muonBaseColl.at(1).Pt() > 10) ){

          for(int i=0; i<muonBaseColl.size(); i++){
            if(eventbase->GetMuonSel()->MuonPass(muonBaseColl.at(i), "MUON_HN_TRI_TIGHT")){
              muonObsColl.push_back(muonBaseColl.at(i));
            }
          }

          if( muonObsColl.size() == 2 ){

            FillHist("n_Observed_electrons", 0., this_weight, 0., 1., 1);

            FillHist("Pt_Observed_electrons", muonObsColl.at(0).Pt(), this_weight, 0., 200., 20);
            FillHist("Pt_Observed_electrons", muonObsColl.at(1).Pt(), this_weight, 0., 200., 20);

            FillHist("Eta_Observed_electrons", muonObsColl.at(0).Eta(), this_weight, -3., 3., 10);
            FillHist("Eta_Observed_electrons", muonObsColl.at(1).Eta(), this_weight, -3., 3., 10);

            FillHist("Pt_Observed_leading_electron", muonObsColl.at(0).Pt(), this_weight, 0., 200., 20);
            FillHist("Pt_Observed_subleading_electron", muonObsColl.at(1).Pt(), this_weight, 0., 200., 20);

            FillHist("Eta_Observed_leading_electron", muonObsColl.at(0).Eta(), this_weight, -3., 3., 10);
            FillHist("Eta_Observed_subleading_electron", muonObsColl.at(1).Eta(), this_weight, -3., 3., 10);

          }
	}
      }
    }
  }

  if( (muonBaseColl.size() == 2) ){
    if( (muonBaseColl.at(0).Charge() == muonBaseColl.at(1).Charge()) ){
      if( (muonBaseColl.at(0).Pt() > 20) && (muonBaseColl.at(1).Pt() > 10) ){

        for(int i=0; i<muonBaseColl.size(); i++){
          muonPredColl.push_back(muonBaseColl.at(i));
        }
	

        m_datadriven_bkg->GetFakeObj()->SetUseQCDFake(true);
        this_weight = this_weight*( m_datadriven_bkg->Get_DataDrivenWeight(false, muonPredColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", 0));

        FillHist("n_Predicted_electrons", 0., this_weight, 0., 1., 1);

        FillHist("Pt_Predicted_electrons", muonPredColl.at(0).Pt(), this_weight, 0., 200., 20);
        FillHist("Pt_Predicted_electrons", muonPredColl.at(1).Pt(), this_weight, 0., 200., 20);

        FillHist("Eta_Predicted_electrons", muonPredColl.at(0).Eta(), this_weight, -3., 3., 10);
        FillHist("Eta_Predicted_electrons", muonPredColl.at(1).Eta(), this_weight, -3., 3., 10);

        FillHist("Pt_Predicted_leading_electron", muonPredColl.at(0).Pt(), this_weight, 0., 200., 20);
        FillHist("Pt_Predicted_subleading_electron", muonPredColl.at(1).Pt(), this_weight, 0., 200., 20);

        FillHist("Eta_Predicted_leading_electron", muonPredColl.at(0).Eta(), this_weight, -3., 3., 10);
        FillHist("Eta_Predicted_subleading_electron", muonPredColl.at(1).Eta(), this_weight, -3., 3., 10);

      }
    }
  }
 
  

  return;
}

/*
void HNSSSFMuMuE_CR_FR::MCClosure(void){

  double this_weight = 1;
  
  this_weight*=MCweight;

  bool pass_ee_trig = false;
  if(PassTrigger("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v")) pass_ee_trig = true;
  if( !pass_ee_trig ) return;

  std::vector<snu::KElectron> electronPRColl = GetHNElectronsByLooseRelIso(0.5, false);
  std::vector<snu::KElectron> electronBaseColl = GetHNElectronsByLooseRelIso(0.5, true);
  std::vector<snu::KElectron> electronPredColl;
  std::vector<snu::KElectron> electronObsColl;
  electronPredColl.clear();
  electronObsColl.clear();

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",true);
  if( muonLooseColl.size() > 0 ) return;

  FillHist("N_Prompt_El", electronPRColl.size(), weight, 0., 10., 10);

  if( (electronPRColl.size() < 2) ){

    FillHist("N_Prompt_El_cut2", electronPRColl.size(), weight, 0., 10., 10);
    FillHist("N_Base_El", electronBaseColl.size(), weight, 0., 10., 10);

    if( (electronBaseColl.size() != 2) ) return;

    FillHist("N_Base_El_cut2", electronBaseColl.size(), weight, 0., 10., 10);

    bool is_SS = false;
    if( (electronBaseColl.at(0).Charge() != electronBaseColl.at(1).Charge()) ) return;//is_SS = true;

    if( (electronBaseColl.at(0).Pt() < 20) || (electronBaseColl.at(1).Pt() < 15) ) return;

    for(int i=0; i<electronBaseColl.size(); i++){
      if(eventbase->GetElectronSel()->ElectronPass(electronBaseColl.at(i), "ELECTRON_HN_LOWDXY_TIGHT")){
        electronObsColl.push_back(electronBaseColl.at(i));
      }
      electronPredColl.push_back(electronBaseColl.at(i));
    }

    if( electronObsColl.size() == 2 ){

      FillHist("n_Observed_electrons", 0., this_weight, 0., 1., 1);

      FillHist("Pt_Observed_electrons", electronObsColl.at(0).Pt(), this_weight, 0., 200., 20);
      FillHist("Pt_Observed_electrons", electronObsColl.at(1).Pt(), this_weight, 0., 200., 20);

      FillHist("Eta_Observed_electrons", electronObsColl.at(0).Eta(), this_weight, -3., 3., 10);
      FillHist("Eta_Observed_electrons", electronObsColl.at(1).Eta(), this_weight, -3., 3., 10);

      FillHist("Pt_Observed_leading_electron", electronObsColl.at(0).Pt(), this_weight, 0., 200., 20);
      FillHist("Pt_Observed_subleading_electron", electronObsColl.at(1).Pt(), this_weight, 0., 200., 20);
    
      FillHist("Eta_Observed_leading_electron", electronObsColl.at(0).Eta(), this_weight, -3., 3., 10);
      FillHist("Eta_Observed_subleading_electron", electronObsColl.at(1).Eta(), this_weight, -3., 3., 10);

    }

    if( electronPredColl.size() == 2 ){
      if( electronObsColl.size() == 2 ) return;
 
      m_datadriven_bkg->GetFakeObj()->SetUseQCDFake(true); 
      this_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 0, electronPredColl, "ELECTRON_HN_LOWDXY_TIGHT", 2);

      FillHist("n_Predicted_electrons", 0., this_weight, 0., 1., 1);

      FillHist("Pt_Predicted_electrons", electronPredColl.at(0).Pt(), this_weight, 0., 200., 20);
      FillHist("Pt_Predicted_electrons", electronPredColl.at(1).Pt(), this_weight, 0., 200., 20);

      FillHist("Eta_Predicted_electrons", electronPredColl.at(0).Eta(), this_weight, -3., 3., 10);
      FillHist("Eta_Predicted_electrons", electronPredColl.at(1).Eta(), this_weight, -3., 3., 10);

      FillHist("Pt_Predicted_leading_electron", electronPredColl.at(0).Pt(), this_weight, 0., 200., 20);
      FillHist("Pt_Predicted_subleading_electron", electronPredColl.at(1).Pt(), this_weight, 0., 200., 20);

      FillHist("Eta_Predicted_leading_electron", electronPredColl.at(0).Eta(), this_weight, -3., 3., 10);
      FillHist("Eta_Predicted_subleading_electron", electronPredColl.at(1).Eta(), this_weight, -3., 3., 10);

    }
  }
}*/
