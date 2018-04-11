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

  TruthPrintOut();

  //if( !(k_sample_name.Contains("Tchannel"))) return;
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

//  TruthPrintOut();
//  return;

  int electron_ID = 11, muon_ID = 13, W_ID = 24, HN_ID = 9900012;
  int parton_ID = 1;

  //TruthPrintOut();

  int max = truthColl.size();
  vector<int> HN_index, onshell_W_index, lep1_index, lep2_index, quark1_index, quark2_index, vbf_quark_index;
  HN_index.clear(); onshell_W_index.clear(); lep1_index.clear(); lep2_index.clear(); quark1_index.clear(); quark2_index.clear(); vbf_quark_index.clear();

  // Look for Heavy Neutrino using PdgId
  for( int i = 2 ; i < max ; i++ ){
    if( abs(truthColl.at(i).PdgId()) == HN_ID ){
      HN_index.push_back(i);
      GENFindDecayIndex( truthColl, i, HN_index );
      break;
    }
  }
  if( HN_index.size() == 0 ){
    FillHist("GEN_HN_not_found", 0., 1., 0., 1., 1);
    return;
  }
  int HN_mother_index = truthColl.at( HN_index.at(0) ).IndexMother(); // Save HN mother index for lep1

  // Look for muon1 using PdgId and mother index (sister : HN)
  for( int i = 2 ; i < max ; i++ ){
    if( ((abs(truthColl.at(i).PdgId()) == muon_ID) || (abs(truthColl.at(i).PdgId()) == electron_ID)) && truthColl.at(i).IndexMother() == HN_mother_index ){
      lep1_index.push_back(i);
      GENFindDecayIndex( truthColl, i, lep1_index);
      break;
    }
  }
  if( lep1_index.size () == 0 ){
    FillHist("GEN_lep1_not_found", 0., 1., 0., 1., 1);
    return;
  }

  // Look for muon2 using PdgId and mother index (mother : HN)
  int HN_decay_index = 0;
  for( int index = 0 ; index < HN_index.size() ; index++ ){ // One of the HN index decays into muon
    for( int i = 2 ; i < max ; i++ ){
      if( ((abs(truthColl.at(i).PdgId()) == muon_ID) || (abs(truthColl.at(i).PdgId()) == electron_ID)) && truthColl.at(i).IndexMother() == HN_index.at(index) ){
        lep2_index.push_back(i);
        GENFindDecayIndex( truthColl, i, lep2_index);
        HN_decay_index = index;
        break;
      }
    }
  }

  snu::KTruth GEN_lep[2], GEN_quark[2], GEN_onshellW, GEN_HN, GEN_vbfquark;
  GEN_lep[0] = truthColl.at( lep1_index.at(0) );
  GEN_lep[1] = truthColl.at( lep2_index.at(0) );

  FillHist("pt1", GEN_lep[0].Pt(), 1., 0., 500., 500);
  FillHist("pt2", GEN_lep[1].Pt(), 1., 0., 500., 500);
  if(GEN_lep[0].PdgId() * GEN_lep[1].PdgId() > 0) FillHist("SSvsOS", 0., 1., 0., 2., 2);
  if(GEN_lep[0].PdgId() * GEN_lep[1].PdgId() < 0) FillHist("SSvsOS", 1., 1., 0., 2., 2);



 
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

void ExampleAnalyzer::GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index ){

  for( int i = it+1 ; i < truthColl.size(); i ++ ){
    if( truthColl.at(i).IndexMother() == it && truthColl.at(i).PdgId() == truthColl.at(it).PdgId() ){
      index.push_back(i);
    }
  }
  return;

}

