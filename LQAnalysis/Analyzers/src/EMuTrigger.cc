// $Id: EMuTrigger.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQEMuTrigger Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "EMuTrigger.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (EMuTrigger);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
EMuTrigger::EMuTrigger() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("EMuTrigger");
  
  Message("In EMuTrigger constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();


}


void EMuTrigger::InitialiseAnalysis() throw( LQError ) {
  
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


void EMuTrigger::ExecuteEvents()throw( LQError ){

  if(!PassMETFilter()) return;     /// Initial event cuts :
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex

  std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_HN_VETO", true);
  std::vector<snu::KMuon> muons = GetMuons("MUON_HN_TIGHT", true);

  std::vector<snu::KElectron> electronVetoColl = GetElectrons(false,true,"ELECTRON_HN_VETO");
  std::vector<snu::KElectron> electrons = GetElectrons(false,true,"ELECTRON_HN_TIGHTv4");

  if(muonVetoColl.size() != 1) return;
  if(electronVetoColl.size() != 1) return;

  if(muons.size() != 1) return;
  if(electrons.size() != 1) return;

  if(muons.at(0).Charge() == electrons.at(0).Charge()) return;

  if(!isData) weight*=MCweight;
  std::vector<TString> muontag; muontag.clear();
  muontag.push_back("HLT_IsoMu24_v");
  
  if(!PassTriggerOR(muontag)) return;
  if(muons.at(0).Pt() < 27.) return;
  if(fabs(muons.at(0).Eta()) > 2.1) return;

  std::vector<TString> electronprobe_leg1, electronprobe_leg2;
  electronprobe_leg1.clear(); electronprobe_leg2.clear();
  electronprobe_leg1.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  electronprobe_leg2.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");

  float ptarray[] = {10,15,20,25,30,50,90,150,500};
  float etaarray[] = {-2.5,-2.0,-1.566,-1.4442, -0.8, 0.0, 0.8, 1.4442, 1.566, 2.0, 2.5};

  FillHist("ALL_den_electron_pt", electrons.at(0).Pt(), weight, 0., 200., 200);
  FillHist("ALL_den_electron_eta", electrons.at(0).SCEta(), weight, -3., 3., 60);
  FillHist("ALL_den_electron_pteta", electrons.at(0).SCEta(), electrons.at(0).Pt(), weight, etaarray, 10, ptarray, 8);


  if(PassTriggerOR(electronprobe_leg1)){

    TString string = "LEG1";
    FillHist(string + "_num_electron_pt", electrons.at(0).Pt(), weight, 0., 200., 200);
    FillHist(string + "_num_electron_eta", electrons.at(0).SCEta(), weight, -3., 3., 60);
    FillHist(string + "_num_electron_pteta", electrons.at(0).SCEta(), electrons.at(0).Pt(), weight, etaarray, 10, ptarray, 8);

  } 

  if(PassTriggerOR(electronprobe_leg2)){

    TString string = "LEG2";
    FillHist(string + "_num_electron_pt", electrons.at(0).Pt(), weight, 0., 200., 200);
    FillHist(string + "_num_electron_eta", electrons.at(0).SCEta(), weight, -3., 3., 60);
    FillHist(string + "_num_electron_pteta", electrons.at(0).SCEta(), electrons.at(0).Pt(), weight, etaarray, 10, ptarray, 8);

  }
 

 
  return;
}// End of execute event loop
  


void EMuTrigger::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void EMuTrigger::BeginCycle() throw( LQError ){
  
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

EMuTrigger::~EMuTrigger() {
  
  Message("In EMuTrigger Destructor" , INFO);
  
}


void EMuTrigger::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void EMuTrigger::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this EMuTriggerCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void EMuTrigger::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


