// $Id: FakeRateCalculator_Mu.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFakeRateCalculator_Mu Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "FakeRateCalculator_Mu.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FakeRateCalculator_Mu);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
FakeRateCalculator_Mu::FakeRateCalculator_Mu() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("FakeRateCalculator_Mu");
  
  Message("In FakeRateCalculator_Mu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void FakeRateCalculator_Mu::InitialiseAnalysis() throw( LQError ) {
  
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


void FakeRateCalculator_Mu::ExecuteEvents()throw( LQError ){

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
   if (!k_isdata) {   pileup_reweight = TempPileupWeight();}
 
  TString trig_mu8="HLT_Mu8_v";
  TString trig_mu17="HLT_Mu17_v";
  TString trig_mu50="HLT_Mu50_v";

  vector<TString> trig_list;
  trig_list.push_back(trig_mu8);
  trig_list.push_back(trig_mu17);
  trig_list.push_back(trig_mu50);

  std::vector<snu::KJet> jetColl =   GetJets("JET_HN");
  int nbjet = NBJet(GetJets("JET_HN"));

  std::vector<snu::KMuon> muonLooseColl =GetMuons("MUON_HN_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl =GetMuons("MUON_HN_TIGHT",false); 

  bool trigpass_mu8 = PassTrigger("HLT_Mu8_v");
  bool trigpass_mu17 = PassTrigger("HLT_Mu17_v");
  bool trigpass_mu50 = PassTrigger("HLT_Mu50_v");

  if( !(trigpass_mu8) && !(trigpass_mu17) && !(trigpass_mu50) ) return;

  CorrectMuonMomentum(muonLooseColl);
   
  double ev_weight = weight;
  if(!isData){
    //weight = weight * trigger_sf;
    //ev_weight = w * trigger_sf * id_iso_sf *  pu_reweight*trigger_ps;
  }

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.PFMET();

  double HT=0.;
  for(int i=0; i<jetColl.size(); i++){
    HT+=jetColl.at(i).Pt();
  }

  if(muonLooseColl.size() == 0) return;

  // Single Muon Trigger used for study
  if( trigpass_mu8 || trigpass_mu17 ){

    weight = weight*GetPrescale(muonLooseColl, trigpass_mu8, trigpass_mu17);

    FillHist("SingleMuonTrigger_NEvents_TrigPass", 0., weight, 0., 1., 1);
    FillHist("SingleMuonTrigger_NMuons_TrigPass", muonLooseColl.size(), weight, 0., 10., 10);
    FillHist("SingleMuonTrigger_NJets_TrigPass", jetColl.size(), weight, 0., 10., 10);

    // Use Dijet enriched events
    if( !(jetColl.size() == 0) ){

      for(int i=0; i<muonLooseColl.size(); i++){
        snu::KMuon muon;
        muon=muonLooseColl.at(i);

	StudyMuon("withJets", muon, weight); 

        FillHist("SingleMuonTrigger_NEvents_withJets", 0., weight, 0., 1., 1);
        FillHist("SingleMuonTrigger_NMuons_withJets", muonLooseColl.size(), weight, 0., 10., 10);
        FillHist("SingleMuonTrigger_NJets_withJets", jetColl.size(), weight, 0., 10., 10);
      }

      for(int i=0; i<jetColl.size(); i++){
	snu::KParticle jet;
	jet=jetColl.at(i);

	FillHist("SingleMuonTrigger_JetPt_withJets", jet.Pt(), weight, 0., 500., 500);
	FillHist("SingleMuonTrigger_JetEta_withJets", jet.Eta(), weight, -3., 3., 60);
	FillHist("SingleMuonTrigger_HT_withJets", HT, weight, 0., 500., 500);
      }

      if( muonLooseColl.size() == 1 ){
        snu::KMuon muon;
        muon=muonLooseColl.at(0);

	StudyMuon("oneMuonwithJets", muon, weight);

        FillHist("SingleMuonTrigger_NEvents_oneMuonwithJets", 0., weight, 0., 1., 1);
        FillHist("SingleMuonTrigger_NMuons_oneMuonwithJets", muonLooseColl.size(), weight, 0., 10., 10);
        FillHist("SingleMuonTrigger_NJets_oneMuonwithJets", jetColl.size(), weight, 0., 10., 10);

	float dR=-999., dPhi=-999., Pt_mu_over_jet=-999.;
        for( int i=0; i<jetColl.size(); i++){
	  snu::KParticle jet;
	  jet=jetColl.at(i);

   	  dR = jet.DeltaR(muon);
	  dPhi = jet.DeltaPhi(muon);
	  Pt_mu_over_jet = ( muon.Pt() / jet.Pt() );

	  if( dR > 1. ){
	    if( dPhi > 2.5 ){

	      StudyMuon("F0", muon, weight);
	      FillHist("SingleMuonTrigger_NMuons_F0", muonLooseColl.size(), weight, 0., 10., 10);
	      FillHist("SingleMuonTrigger_NJets_F0", jetColl.size(), weight, 0., 10., 10);

	      if( muonTightColl.size() == 1 ){

		StudyMuon("F", muon, weight);
                FillHist("SingleMuonTrigger_NMuons_F", muonLooseColl.size(), weight, 0., 10., 10);
                FillHist("SingleMuonTrigger_NJets_F", jetColl.size(), weight, 0., 10., 10);

	      }

	      goto fill_F0;
	    }
	  }
	}

	fill_F0:;
      } 
    }
  }

}// End of execute event loop
  


void FakeRateCalculator_Mu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void FakeRateCalculator_Mu::BeginCycle() throw( LQError ){
  
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

FakeRateCalculator_Mu::~FakeRateCalculator_Mu() {
  
  Message("In FakeRateCalculator_Mu Destructor" , INFO);
  
}


void FakeRateCalculator_Mu::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void FakeRateCalculator_Mu::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FakeRateCalculator_MuCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FakeRateCalculator_Mu::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


double FakeRateCalculator_Mu::GetPrescale(std::vector<snu::KMuon> muonColl, bool passlow, bool passhigh){

  double prescale_trigger = 0.;

  if( muonColl.at(0).Pt() >= 20. ){
    if( passhigh ){
      prescale_trigger = WeightByTrigger("HLT_Mu17_v", TargetLumi);
    }
    else prescale_trigger = 0.;
  }
  else{
    if( passlow ){
      prescale_trigger = WeightByTrigger("HLT_Mu8_v", TargetLumi);
    }
    else prescale_trigger = 0.;
  }

  if(isData && prescale_trigger == 0) return prescale_trigger;
  else if(isData && !(prescale_trigger == 0)) return 1.;
  else if(!(isData) && (prescale_trigger == 0)) return prescale_trigger;
  else if(!(isData) && !(prescale_trigger == 0)) return prescale_trigger;
  else return -999999.;

}


void FakeRateCalculator_Mu::StudyMuon(TString suffix, snu::KMuon muon, double weight){

//  FillHist("SingleMuonTrigger_NEvents_"+suffix, 0., weight, 0., 1., 1);
  FillHist("SingleMuonTrigger_MuonPt_"+suffix, muon.Pt(), weight, 0., 500., 500);
  FillHist("SingleMuonTrigger_MuonEta_"+suffix, muon.Eta(), weight, -3., 3., 60);
  FillHist("SingleMuonTrigger_MuondXY_"+suffix, muon.dXY(), weight, -0.5, 0.5, 100);
  FillHist("SingleMuonTrigger_MuondZ_"+suffix, muon.dZ(), weight, -1., 1., 200);
  FillHist("SingleMuonTrigger_MuonRelIso04_"+suffix, muon.RelIso04(), weight, 0., 1., 100);
//  FillHist("SingleMuonTrigger_MET_"+suffix, MET, weight, 0., 500., 500);
//  FillHist("SingleMuonTrigger_HT_"+suffix, HT, weight, 0., 500., 500);

  return;

}
