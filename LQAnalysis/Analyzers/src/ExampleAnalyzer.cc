// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQExampleAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ExampleAnalyzer.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ExampleAnalyzer);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ExampleAnalyzer::ExampleAnalyzer() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ExampleAnalyzer");
  
  Message("In ExampleAnalyzer constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void ExampleAnalyzer::InitialiseAnalysis() throw( LQError ) {
  
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


void ExampleAnalyzer::ExecuteEvents()throw( LQError ){

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
     
   
  TString dimuon_trigmuon_trig1="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
  TString dimuon_trigmuon_trig2="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v";
  TString dimuon_trigmuon_trig3="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v";
  TString dimuon_trigmuon_trig4="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
   // Now you should do an OR of 4 triggers 
 
  vector<TString> trignames;
  trignames.push_back(dimuon_trigmuon_trig1);
  trignames.push_back(dimuon_trigmuon_trig2);
  trignames.push_back(dimuon_trigmuon_trig3);
  trignames.push_back(dimuon_trigmuon_trig4);

  if(!PassTriggerOR(triggerlist)) return;   

  std::vector<snu::KMuon> muonL = GetMuons("MUON_HN_TRI_LOOSE",false); 
  std::vector<snu::KMuon> muonT = GetMuons("MUON_HN_TRI_TIGHT",false);

  std::vector<snu::KElectron> electronTMVA = GetElectrons(false,false, "ELECTRON_MVA_TIGHT");
  std::vector<snu::KElectron> electronTHNMVA = GetElectrons(false,false, "ELECTRON_HN_TIGHT");
  std::vector<snu::KElectron> electronTGENT = GetElectrons(false, false, "ELECTRON_GENT_TIGHT");

  std::vector<snu::KElectron> electronLMVA = GetElectrons(false,false, "ELECTRON_MVA_FAKELOOSE");
  std::vector<snu::KElectron> electronLHNMVA = GetElectrons(false,false, "ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronLGENT = GetElectrons(false, false, "ELECTRON_GENT_FAKELOOSE");

  if(muonL.size() != 2) return;

  if((muonL.at(0).Charge() != muonL.at(1).Charge())) return;
  if((muonL.at(0).Pt() < 20 || muonL.at(1).Pt() < 10)) return;

  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN", 30., 2.4);
  for(int j=0; j<jetTightColl.size(); j++){
    if(jetTightColl.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) return;
  } 

  if(k_running_nonprompt){
    if(electronLMVA.size() == 1){
      weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonL, "MUON_HN_TRI_TIGHT", 2, electronLMVA, "ELECTRON_MVA_TIGHT", 1, "ELECTRON_MVA_FAKELOOSE", "dijet_ajet40");
      FillHist("N_ELECTRON_TIGHT", 0., weight, 0., 2., 2);
    }
    if(electronLHNMVA.size() == 1){
      weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonL, "MUON_HN_TRI_TIGHT", 2, electronLHNMVA, "ELECTRON_HN_TIGHT", 1, "ELECTRON_HN_FAKELOOSE", "mva");
      FillHist("N_ELECTRON_TIGHT", 1., weight, 0., 2., 2);
    }
  }
  else{
    if(k_sample_name.Contains("HN")) weight =1.;

    if(electronTMVA.size() == 1 && muonT.size() == 2){
      FillHist("N_ELECTRON_TIGHT", 0., weight, 0., 2., 2);
    }
    if(electronTHNMVA.size() == 1 && muonT.size() == 2){
      FillHist("N_ELECTRON_TIGHT", 1., weight, 0., 2., 2);
    }
  }
//  if(electronGENT.size() == 1) FillHist("N_ELECTRON_TIGHT", 2., 1., 0., 3., 3);

   //   std::vector<snu::KElectron> electrons2 =  GetElectrons(BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0);

	    
   
   return;
}// End of execute event loop
  


void ExampleAnalyzer::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ExampleAnalyzer::BeginCycle() throw( LQError ){
  
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

ExampleAnalyzer::~ExampleAnalyzer() {
  
  Message("In ExampleAnalyzer Destructor" , INFO);
  
}


void ExampleAnalyzer::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ExampleAnalyzer::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ExampleAnalyzerCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void ExampleAnalyzer::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



