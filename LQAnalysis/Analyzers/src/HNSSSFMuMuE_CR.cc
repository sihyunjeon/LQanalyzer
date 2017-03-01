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
  TString mumu_trigger="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
  vector<TString> trignames;
  trignames.push_back(mumu_trigger);

  float weight_trigger = WeightByTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", TargetLumi);

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_HN", 15.);
  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN", 20.);

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TRI_TIGHT",false);

  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,false,"ELECTRON_HN_TIGHT");

  bool trig_pass=PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  if(!trig_pass) return;

  CorrectMuonMomentum(muonLooseColl);
  float electron_idsf = mcdata_correction->ElectronScaleFactor("ELECTRON_HN_FAKELOOSE", electronLooseColl);
  float electron_reco = mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);

  double ev_weight = weight;
  if(!isData){
    weight *= weight_trigger;
    weight *= pileup_reweight;
    weight *= electron_idsf;
    weight *= electron_reco;
    //ev_weight = w * trigger_sf * id_iso_sf *  pu_reweight*trigger_ps;
  }
 
  int period_index = 0;
  period_index = GetPeriodIndex();

  int nbjet = NBJet(jetLooseColl, snu::KJet::CSVv2, snu::KJet::Loose, period_index);
  double HT = 0.;
  for( int i=0; i<jetLooseColl.size(); i++){
    HT += jetLooseColl.at(i).Pt();
  }

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  snu::KParticle MET;
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);

  double Z_mass = 91.1876;

  bool is_mumue = ( ((muonLooseColl.size() == 2) && (electronLooseColl.size() == 1)) && ((muonTightColl.size() == 2) && (electronTightColl.size() == 1)) );
  
  // =============================================================================
  // == MuMuE Selection===========================================================
  // ====================================

  if( is_mumue ){

    SF[0] = muonLooseColl.at(0);
    SF[1] = muonLooseColl.at(1);
    OF = electronLooseColl.at(0);

    CR_WZ_mumue=false; CR_Zjet_mumue=false; CR_ttW_mumue=false;

    // == OS muon pair cut
    bool p_lepton_OSSF = ((SF[0].Charge()) != (SF[1].Charge()));
    // == Z candidate mass cut
    snu::KParticle Z_candidate;
    Z_candidate = SF[0] + SF[1];
    bool p_Z_candidate_mass = ((fabs(Z_candidate.M() - Z_mass) < 10));
    // == lepton Pt cut
    bool p_lepton_pt = (((SF[0].Pt() > 20) && (SF[1].Pt() > 10)) && (OF.Pt() > 10));
    // == MET cut
    bool p_MET_30 = (MET.Pt() > 30);
    bool p_MET_20 = (MET.Pt() > 20);
    // == trilepton mass cut
    bool p_trilepton_mass = ( ((SF[0]+SF[1]+OF).M() > 100) );
    // == jet cut
    bool p_jet_N_2 = ( (jetTightColl.size() > 1) );
    // == bjet cut
    bool p_bjet_N_0 = ( nbjet == 0 );

    if(  p_lepton_OSSF  &&  p_Z_candidate_mass  &&  p_lepton_pt  &&  p_MET_30  &&  p_trilepton_mass                                ) CR_WZ_mumue = true;
    if(  p_lepton_OSSF  &&  p_Z_candidate_mass  &&  p_lepton_pt  && !p_MET_20  &&  p_trilepton_mass                                ) CR_Zjet_mumue = true;
    if( !p_lepton_OSSF                          &&  p_lepton_pt  &&  p_MET_30                        &&  p_jet_N_2  && !p_bjet_N_0 ) CR_ttW_mumue = true;

    if( p_lepton_OSSF && p_lepton_pt ){
      DrawHistograms("2mu1e", SF, OF, MET, jetLooseColl, weight);
      FillCLHist(sssf_mumue, "2mu1e", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, weight);
    }

    // ========== WZ selection =====================================
    if( CR_WZ_mumue ){
      DrawHistograms("2mu1e_WZ", SF, OF, MET, jetLooseColl, weight);
      FillCLHist(sssf_mumue, "2mu1e_WZ", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, weight);
    }

    // ========== Z+lepton selection =====================================
    if( CR_Zjet_mumue ){
      DrawHistograms("2mu1e_Zjet", SF, OF, MET, jetLooseColl, weight);
      FillCLHist(sssf_mumue, "2mu1e_Zjet", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, weight);
    }

    // ========== ttW selection =====================================
    if( CR_ttW_mumue ){
      DrawHistograms("2mu1e_ttW", SF, OF, MET, jetLooseColl, weight);
      FillCLHist(sssf_mumue, "2mu1e_ttW", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetLooseColl, weight);
    }

    FillHist("weight", weight, weight, 0., 2., 2000);
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

  int period_index = 0;
  period_index = GetPeriodIndex();
  int nbjet = NBJet(jetLooseColl, snu::KJet::CSVv2, snu::KJet::Loose, period_index);
  FillHist("NBJets_"+suffix, nbjet, weight, 0., 5., 5);

  return;

}


void HNSSSFMuMuE_CR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
 
  SF[0].SetPxPyPzE(0,0,0,0); SF[1].SetPxPyPzE(0,0,0,0); OF.SetPxPyPzE(0,0,0,0);
  CR_WZ_mumue=false; CR_Zjet_mumue=false; CR_ttW_mumue=false;
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
