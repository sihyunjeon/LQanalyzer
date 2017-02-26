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

  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton
  
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                                              
  float pileup_reweight=(1.0);
  if (!k_isdata) {   pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());}
    
  TString mumu_trigger="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
  vector<TString> trignames;
  trignames.push_back(mumu_trigger);

  float weight_trigger = WeightByTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", TargetLumi);

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_NOCUT");
  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
  int nbjet = NBJet(GetJets("JET_HN"));

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TRI_TIGHT",false);

  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,false,"ELECTRON_HN_TIGHT");

  bool trig_pass=PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  if(!trig_pass) return;
  CorrectMuonMomentum(muonLooseColl);
   
  double ev_weight = weight;
  if(!isData){
    weight *= weight_trigger;
    weight *= pileup_reweight;
    //ev_weight = w * trigger_sf * id_iso_sf *  pu_reweight*trigger_ps;
  }


  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();

  snu::KParticle MET;
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);

  bool is_WGtoLNuEE=false, is_WGtoLNuMM=false;
  if(k_sample_name.Contains("WGtoLNuEE")) is_WGtoLNuEE=true;
  if(k_sample_name.Contains("WGtoLNuMM")) is_WGtoLNuMM=true;

  if( muonLooseColl.size() == 0 || electronLooseColl.size() == 0 ) return;

  // Control Region divided in 3 classes
  // =============================================================================
  // MuMuE
  // =========================================
  // == CLASS 1
  // ========== WZ selection SMP-16-002
  // =========================================
  // == CLASS 2
  // ========== Z(mumu) + nonprompt lepton(electron)
  //
  //
  // EEMu
  // =========================================
  // == CLASS 3
  // ========== Z(elel) + nonprompt lepton(muon)


  
  // =============================================================================
  // == MuMuE Selection===========================================================
  // ====================================

  if( (muonLooseColl.size() == 2 && muonTightColl.size() == 2) && (electronLooseColl.size() == 1 && electronTightColl.size() == 1) && !is_WGtoLNuEE ){
    
    SF[0] = muonLooseColl.at(0);
    SF[1] = muonLooseColl.at(1);
    OF = electronLooseColl.at(0);

    if( SF[0].Charge() != SF[1].Charge() ){    

      snu::KParticle Z_candidate;
      Z_candidate = SF[0] + SF[1];

      bool p_Z_candidate_mass = ((Z_candidate.M() > 81) && (Z_candidate.M() < 101));
      bool p_lepton_pt = (((SF[0].Pt() > 20) && (SF[1].Pt() > 10)) && (OF.Pt() > 20));
      bool p_MET_pt = (MET.Pt() > 30);
      bool p_dilep_mass = ( ((SF[0]+SF[1]).M() > 4) && ((SF[0]+OF).M() > 4) && ((SF[1]+OF).M() > 4) );
      bool p_trilep_mass = ( ((SF[0]+SF[1]+OF).M() > 100) );
      bool p_bjet = ( nbjet == 0 );

      if(OF.Pt() < 20) return;
      DrawHistograms("2mu1e", SF, OF, MET, jetTightColl, weight);
      FillCLHist(sssf_mumue, "2mu1e", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);

      // ========== WZ selection =====================================
      if( p_Z_candidate_mass && p_lepton_pt && p_MET_pt && p_dilep_mass && p_trilep_mass && p_bjet ){
        DrawHistograms("2mu1e_WZ", SF, OF, MET, jetTightColl, weight);
        FillCLHist(sssf_mumue, "2mu1e_WZ", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);
      }

      // ========== Z+lepton selection =====================================
      if( p_Z_candidate_mass && p_lepton_pt && !p_MET_pt && p_dilep_mass && p_trilep_mass && p_bjet ){
        DrawHistograms("2mu1e_Zjet", SF, OF, MET, jetTightColl, weight);
        FillCLHist(sssf_mumue, "2mu1e_Zjet", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);
      }
    }
  }

/*

  // =============================================================================
  // == EEMu Selection============================================================
  // ====================================

  if( (muonLooseColl.size() == 1 && muonTightColl.size() == 1) && (electronLooseColl.size() == 2 && electronTightColl.size() == 2) && !is_WGtoLNuMM){


    SF[0] = electronLooseColl.at(0);
    SF[1] = electronLooseColl.at(1);
    OF = muonLooseColl.at(0);

    if( SF[0].Charge() != SF[1].Charge() ){

      snu::KParticle Z_candidate;
      Z_candidate = SF[0] + SF[1];

      bool p_Z_candidate_mass = ((Z_candidate.M() > 76) && (Z_candidate.M() < 106));
      bool p_lepton_pt = (((SF[0].Pt() > 20) && (SF[1].Pt() > 10)) && (OF.Pt() > 20));
      bool p_MET_pt = (MET.Pt() > 30);
      bool p_dilep_mass = ( ((SF[0]+SF[1]).M() > 4) && ((SF[0]+OF).M() > 4) && ((SF[1]+OF).M() > 4) );
      bool p_trilep_mass = ( ((SF[0]+SF[1]+OF).M() > 100) );
      bool p_bjet = ( nbjet == 0 );

      DrawHistograms("1mu2e", SF, OF, MET, jetTightColl, weight);

      // ========== Z+lepton selection =====================================
      if( p_Z_candidate_mass && p_lepton_pt && !p_MET_pt && p_dilep_mass && p_trilep_mass && p_bjet ){
        DrawHistograms("1mu2e_Zjet", SF, OF, MET, jetTightColl, weight);
      }
    }
  }
*/

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

  SF[0].SetPxPyPzE(0,0,0,0); SF[1].SetPxPyPzE(0,0,0,0); OF.SetPxPyPzE(0,0,0,0); 
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


void HNSSSFMuMuE_CR::DrawHistograms(TString suffix, snu::KParticle SF[], snu::KParticle OF, snu::KParticle MET,  std::vector<snu::KJet> jetTightColl, double weight){

  FillHist("number_of_events_"+suffix, 0., weight, 0., 1., 1);

  FillHist("PFMET_"+suffix, MET.Pt(), weight, 0., 500., 500);
  FillHist("NJets_"+suffix, jetTightColl.size(), weight, 0., 5., 5);

  FillHist("W_transverse_mass_"+suffix, MT(OF,MET), weight, 0., 500., 500);
  FillHist("Z_candidate_mass_"+suffix, (SF[0]+SF[1]).M(), weight, 0., 500., 500);

  FillHist("SFleadingLeptonPt_"+suffix, SF[0].Pt(), weight, 0., 500., 500);
  FillHist("SFleadingLeptonEta_"+suffix, SF[0].Eta(), weight, -3., 3., 60);
  FillHist("SFsecondLeptonPt_"+suffix, SF[1].Pt(), weight, 0., 500., 500);
  FillHist("SFsecondLeptonEta_"+suffix, SF[1].Eta(), weight, -3., 3., 60);

  FillHist("OFLeptonPt_"+suffix, OF.Pt(), weight, 0., 500., 500);
  FillHist("OFLeptonEta_"+suffix, OF.Eta(), weight, -3., 3., 60);

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



