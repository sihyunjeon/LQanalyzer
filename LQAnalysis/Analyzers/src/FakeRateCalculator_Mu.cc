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
   if (!k_isdata) {   pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());}
 
  TString trig_mu8="HLT_Mu8_v";
  TString trig_mu17="HLT_Mu17_v";
  TString trig_mu50="HLT_Mu50_v";

  vector<TString> trig_list;
  trig_list.push_back(trig_mu8);
  trig_list.push_back(trig_mu17);
  trig_list.push_back(trig_mu50);

  std::vector<snu::KJet> jetColl =   GetJets("JET_HN");
  int nbjet = NBJet(GetJets("JET_HN"));

  std::vector<snu::KMuon> muonTriLooseColl =GetMuons("MUON_HN_TRI_LOOSE",false);
  std::vector<snu::KMuon> muonTriTightColl =GetMuons("MUON_HN_TRI_TIGHT",false); 

  std::vector<snu::KMuon> muonTriHighdXYLooseColl =GetMuons("MUON_HN_TRI_HIGHDXY_LOOSE",false);
  std::vector<snu::KMuon> muonTriHighdXYTightColl =GetMuons("MUON_HN_TRI_HIGHDXY_TIGHT",false); 

  bool trigpass_mu8 = PassTrigger("HLT_Mu8_v");
  bool trigpass_mu17 = PassTrigger("HLT_Mu17_v");
  bool trigpass_mu50 = PassTrigger("HLT_Mu50_v");

  if( !(trigpass_mu8) && !(trigpass_mu17) ) return;

//  CorrectMuonMomentum(muonTriLooseColl);
   
  double ev_weight = weight;
  if(!isData){
    weight *= pileup_reweight;
    //ev_weight = w * trigger_sf * id_iso_sf *  pu_reweight*trigger_ps;
  }

  snu::KEvent Evt = eventbase->GetEvent();
  double MET = Evt.PFMET();

  double HT=0.;
  for(int i=0; i<jetColl.size(); i++){
    HT+=jetColl.at(i).Pt();
  }

  // Single Muon Trigger used for study
  if( trigpass_mu8 || trigpass_mu17 ){

    //////////////////////////////////////////////////////////
    ///////////////////////// Use Dijet enriched event
    if( muonTriLooseColl.size() == 1 ){

      double dijet_weight = weight*GetPrescale(muonTriLooseColl, trigpass_mu8, trigpass_mu17);

      if( !(jetColl.size() == 0) ){

        snu::KMuon muon;
        muon=muonTriLooseColl.at(0);

	float dR=-999., dPhi=-999., Pt_mu_over_jet=-999.;
        for( int i=0; i<jetColl.size(); i++ ){
	  snu::KParticle jet;
	  jet=jetColl.at(i);

   	  dR = jet.DeltaR(muon);
	  dPhi = jet.DeltaPhi(muon);
	  Pt_mu_over_jet = ( muon.Pt() / jet.Pt() );

	  if( dR > 1. ){
	    if( dPhi > 2.5 ){
	      if( Pt_mu_over_jet < 1. ){

	        StudyDijetMuon("F0", muon, dijet_weight);

	        if( muonTriTightColl.size() == 1 ){

	          StudyDijetMuon("F", muon, dijet_weight);

	        }

	      goto fill_F0;

	      }
	    }
	  }
	}
	fill_F0:;
      }
    }//end of dijet


    /////////////////////////////////////////////////////////////////
    ///////////////////////// Use Loose to Tight method
    if( muonTriHighdXYLooseColl.size() == 1 ){

      double highdxy_weight = weight*GetPrescale(muonTriHighdXYLooseColl, trigpass_mu8, trigpass_mu17);

      snu::KMuon muon;
      muon=muonTriHighdXYLooseColl.at(0);

      StudyHighdXYMuon("F0", muon, highdxy_weight);
      
      if( muonTriHighdXYTightColl.size() == 1 ){

	StudyHighdXYMuon("F", muon, highdxy_weight);

      }
    }//end of highdxy
  }//end of trigger fires

  return;

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


void FakeRateCalculator_Mu::StudyDijetMuon(TString suffix, snu::KMuon muon, double weight){

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

  FillHist("SingleMuonTrigger_Dijet_MuonPt_"+suffix, muon.Pt(), weight, 0., 500., 500);
  FillHist("SingleMuonTrigger_Dijet_MuonEta_"+suffix, muon.Eta(), weight, -3., 3., 60);
  FillHist("SingleMuonTrigger_Dijet_MuondXY_"+suffix, muon.dXY(), weight, -0.5, 0.5, 100);
  FillHist("SingleMuonTrigger_Dijet_MuondXYSig_"+suffix, muon.dXYSig(), weight, -10., 10., 200);
  FillHist("SingleMuonTrigger_Dijet_MuondZ_"+suffix, muon.dZ(), weight, -1., 1., 200);
  FillHist("SingleMuonTrigger_Dijet_MuonRelIso04_"+suffix, muon.RelIso04(), weight, 0., 1., 100);
  FillHist("SingleMuonTrigger_Dijet_NEvents_"+suffix, muon.Pt(), fabs(muon.Eta()), weight, ptarray, 9, etaarray, 4);

  return;

}

void FakeRateCalculator_Mu::StudyHighdXYMuon(TString suffix, snu::KMuon muon, double weight){

  float etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  float ptarray [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

  FillHist("SingleMuonTrigger_HighdXY_MuonPt_"+suffix, muon.Pt(), weight, 0., 500., 500);
  FillHist("SingleMuonTrigger_HighdXY_MuonEta_"+suffix, muon.Eta(), weight, -3., 3., 60);
  FillHist("SingleMuonTrigger_HighdXY_MuondXY_"+suffix, muon.dXY(), weight, -0.5, 0.5, 100);
  FillHist("SingleMuonTrigger_HighdXY_MuondXYSig_"+suffix, muon.dXYSig(), weight, -10., 10., 200);
  FillHist("SingleMuonTrigger_HighdXY_MuondZ_"+suffix, muon.dZ(), weight, -1., 1., 200);
  FillHist("SingleMuonTrigger_HighdXY_MuonRelIso04_"+suffix, muon.RelIso04(), weight, 0., 1., 100);
  FillHist("SingleMuonTrigger_HighdXY_NEvents_"+suffix, muon.Pt(), fabs(muon.Eta()), weight, ptarray, 9, etaarray, 4);

  return;

}
