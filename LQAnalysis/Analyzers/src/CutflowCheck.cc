// $Id: CutflowCheck.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQCutflowCheck Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "CutflowCheck.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (CutflowCheck);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
CutflowCheck::CutflowCheck() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("CutflowCheck");
  
  Message("In CutflowCheck constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void CutflowCheck::InitialiseAnalysis() throw( LQError ) {
  
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


void CutflowCheck::ExecuteEvents()throw( LQError ){

  FillHist("PROCESSED", 0., 1., 0., 1., 1);

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
   if(!PassMETFilter()) return;     /// Initial event cuts : 
   if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex  
   


   float pileup_reweight=(1.0);
   if (!k_isdata) {   pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());}
     
   
  TString dimuon_trigmuon_trig1="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
  TString dimuon_trigmuon_trig2="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
 
  vector<TString> trignames;
  trignames.push_back(dimuon_trigmuon_trig1);
  trignames.push_back(dimuon_trigmuon_trig2);

  if(!PassTriggerOR(trignames)) return;   
  float weight_trigger = WeightByTrigger(trignames, TargetLumi);

  if(!isData){
    weight*=weight_trigger;
    weight*=pileup_reweight;
  }
  FillHist("STEP0_event_cleaning_selection_and_trigger", 0., weight, 0., 1., 1);

  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TIGHT",false);
  if(muonTightColl.size()!=2) return;
  if((muonTightColl.at(0).Pt() < 20.) || (muonTightColl.at(1).Pt() < 10.)) return;
  FillHist("STEP1_exactly_two_hn_tight_muon", 0., weight, 0., 1., 1);

  std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_HN_VETO", false);
  for(int i=0; i<muonVetoColl.size(); i++){
    if((!PassID(muonVetoColl.at(i), "MUON_HN_TIGHT"))) return;
  }
  FillHist("STEP2_loose_muon_veto", 0., weight, 0., 1., 1);

  std::vector<snu::KElectron> electronVetoColl = GetElectrons(false, false, "ELECTRON_HN_VETO");
  if((electronVetoColl.size())!=0) return;
  FillHist("STEP3_electron_veto", 0., weight, 0., 1., 1);
  
  if((muonTightColl.at(0)+muonTightColl.at(1)).M() < 10) return;
  FillHist("STEP4_mass_bigger_than_10", 0., weight, 0., 1., 1);

  std::vector<snu::KJet> jetColl = GetJets("JET_HN");
  if(jetColl.size() <2) return;
  FillHist("STEP5_more_than_1_jet", 0., weight, 0., 1., 1);

  double MET = eventbase->GetEvent().MET();
  MET = CorrectedMETRochester(muonTightColl, true);
  if(MET > 50) return;
  FillHist("STEP6_MET_cut", 0., weight, 0., 1., 1);

  snu::KParticle W_candidate, W_selection=(jetColl.at(0)+jetColl.at(1));
  for(int i=0; i<jetColl.size(); i++){
    for(int j=i+1; j<jetColl.size(); j++){
      W_candidate = (jetColl.at(i)+jetColl.at(j));
      if(fabs(W_selection.M() - 80.4) > fabs(W_candidate.M() - 80.4)){
        W_selection = W_candidate;
      }
      FillHist("W_candidate", W_candidate.M(), 1., 0., 400, 200);
    }
  }
  if(W_selection.M() > 200) return;
  FillHist("W_selection", W_selection.M(), 1., 0., 400, 200);
  FillHist("STEP7_mass_of_jj_smaller_than_200", 0., weight, 0., 1., 1);

  int nbjet=0;
  for(int i=0; i<jetColl.size(); i++){
    if(IsBTagged(jetColl.at(i), snu::KJet::CSVv2, snu::KJet::Medium)){
      nbjet++;
      return;
    }
  }
  FillHist("[check]_n_bjet", nbjet, weight, 0., 5., 5);
  FillHist("STEP8_b_jet_veto_1", 0., weight, 0., 1., 1);

  int nbjet2=0;
  std::vector<snu::KJet> jetNoLeptonVetoColl = GetJets("JET_NOLEPTONVETO");
  for(int i=0; i<jetNoLeptonVetoColl.size(); i++){
    if(IsBTagged(jetNoLeptonVetoColl.at(i), snu::KJet::CSVv2, snu::KJet::Medium)){
      nbjet2++;
      return;
    }
  }
  FillHist("[check]_n_bjet2", nbjet2, weight, 0., 5., 5);
  FillHist("STEP9_b_jet_veto_2", 0., weight, 0., 1., 1);


   return;
}// End of execute event loop
  


void CutflowCheck::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void CutflowCheck::BeginCycle() throw( LQError ){
  
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

CutflowCheck::~CutflowCheck() {
  
  Message("In CutflowCheck Destructor" , INFO);
  
}


void CutflowCheck::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void CutflowCheck::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this CutflowCheckCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void CutflowCheck::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



