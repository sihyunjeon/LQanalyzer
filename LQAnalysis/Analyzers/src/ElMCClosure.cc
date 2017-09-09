// $Id: ElMCClosure.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQElMCClosure Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ElMCClosure.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ElMCClosure);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ElMCClosure::ElMCClosure() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ElMCClosure");
  
  Message("In ElMCClosure constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void ElMCClosure::InitialiseAnalysis() throw( LQError ) {
  
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


void ElMCClosure::ExecuteEvents()throw( LQError ){

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
     
   
  TString diel_trig="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
   // Now you should do an OR of 4 triggers 
 
  vector<TString> trignames;
  trignames.push_back(diel_trig);
  bool trig_pass=PassTriggerOR(trignames);
  if(!trig_pass) return;

  if(!isData){
    weight *=pileup_reweight;
  }

  std::vector<snu::KElectron> electronPromptLooseColl =  GetElectrons(false,false, "ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false,true, "ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,true, "ELECTRON_HN_TIGHT");
  if(electronPromptLooseColl.size() == 3) return;
  if(electronLooseColl.size() != 3) return;

  snu::KParticle el[3];
  for(int i=0;i<3;i++)
    el[i] = electronLooseColl.at(i);

  if(el[0].Pt() < 25) return;
  if(el[1].Pt() < 15) return;
  if(el[2].Pt() < 10) return;

  if(el[0].Charge() == el[1].Charge())
    if(el[1].Charge() == el[2].Charge()) return;

  std::vector<snu::KMuon> muonEmptyColl;
  muonEmptyColl.clear();

  double fake_weight     = m_datadriven_bkg->Get_DataDrivenWeight(false, muonEmptyColl, "MUON_HN_TRI_TIGHT", 0, electronLooseColl, "ELECTRON_HN_TIGHT", 3, "ELECTRON_HN_FAKELOOSE", "mva");
  FillHist("MCClosure_Predicted_NEvents", 0., fake_weight, 0., 1., 1);
  if(electronTightColl.size() == 3){
    FillHist("MCClosure_Observed_NEvents", 0., 1, 0., 1., 1);
  }
   
  return;
}// End of execute event loop
  


void ElMCClosure::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ElMCClosure::BeginCycle() throw( LQError ){
  
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

ElMCClosure::~ElMCClosure() {
  
  Message("In ElMCClosure Destructor" , INFO);
  
}


void ElMCClosure::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ElMCClosure::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ElMCClosureCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void ElMCClosure::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



