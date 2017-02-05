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
  if (!k_isdata) {   pileup_reweight = TempPileupWeight();}
    
  std::vector<snu::KJet> jetLooseColl = GetJets("JET_NOCUT");
  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
  int nbjet = NBJet(GetJets("JET_HN"));

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TIGHT",false);

  std::vector<snu::KElectron> electronLooseColl = GetElectrons("ELECTRON_HN_FAKELOOSE", false);
  std::vector<snu::KElectron> electronTightColl = GetElectrons("ELECTRON_HN_TIGHT", false);

  bool trig_pass_mu9mu9e9 = PassTrigger("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v");
  bool trig_pass_mu17mu8 = PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  bool trig_pass_mu17tkmu8 = PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  bool electron_size = ((electronLooseColl.size() == 1) && (electronTightColl.size() == 1));
  bool muon_size = ((muonLooseColl.size() == 2) && (muonTightColl.size() == 2));

  cout << "e: " << electron_size << " mu: " << muon_size << "   " << trig_pass_mu9mu9e9 << trig_pass_mu17mu8 <<trig_pass_mu17tkmu8 << endl;

  if( electron_size && muon_size ){
    if(trig_pass_mu9mu9e9){
      if( (muonLooseColl.at(0).Pt() > 10) && (muonLooseColl.at(1).Pt() > 10) && (electronLooseColl.at(0).Pt() > 10) ){
        if( muonLooseColl.at(0).Charge() == muonLooseColl.at(1).Charge() ){
          if( muonLooseColl.at(0).Charge() == electronLooseColl.at(0).Charge() ){
            FillHist("DiMu9Ele9_SS-DiMuon_SS-SingleElectron", 0., 1., 0., 1., 1);
          }
          else if( muonLooseColl.at(0).Charge() != electronLooseColl.at(0).Charge() ){
            FillHist("DiMu9Ele9_SS-DiMuon_OS-SingleElectron", 0., 1., 0., 1., 1);
          }
        }
        else if( muonLooseColl.at(0).Charge() != muonLooseColl.at(1).Charge() ){
          FillHist("DiMu9Ele9_OS-DiMuon_SingleElectron", 0., 1., 0., 1., 1);
        }
      }
    }
    
    if(trig_pass_mu17mu8 || trig_pass_mu17tkmu8){
      if( (muonLooseColl.at(0).Pt() > 20) && (muonLooseColl.at(1).Pt() > 10) && (electronLooseColl.at(0).Pt() > 10) ){
        if( muonLooseColl.at(0).Charge() == muonLooseColl.at(1).Charge() ){
          if( muonLooseColl.at(0).Charge() == electronLooseColl.at(0).Charge() ){
            FillHist("Mu17Mu8_SS-DiMuon_SS-SingleElectron", 0., 1., 0., 1., 1);
          }
          else if( muonLooseColl.at(0).Charge() != electronLooseColl.at(0).Charge() ){
            FillHist("Mu17Mu8_SS-DiMuon_OS-SingleElectron", 0., 1., 0., 1., 1);
          }
        }
        else if( muonLooseColl.at(0).Charge() != muonLooseColl.at(1).Charge() ){
          FillHist("Mu17Mu8_OS-DiMuon_SingleElectron", 0., 1., 0., 1., 1);
        }
      }
    }
  }


  

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
