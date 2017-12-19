// $Id: SampleValidation.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQSampleValidation Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "SampleValidation.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (SampleValidation);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
SampleValidation::SampleValidation() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("SampleValidation");
  
  Message("In SampleValidation constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void SampleValidation::InitialiseAnalysis() throw( LQError ) {
  
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


void SampleValidation::ExecuteEvents()throw( LQError ){

  TruthPrintOut();
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

  int max = truthColl.size();

  std::vector<int> hn_index, lep1_index, lep2_index;
  for( int i=2 ; i<max ; i++){
    if( ((truthColl.at(i).PdgId()) == 9900012) ){
      hn_index.push_back(i);
      break;
    }
  }
  if(hn_index.size() == 0) return;
  for( int i=2; i<max ; i++){
    if( (fabs(truthColl.at(i).PdgId()) == 13 || fabs(truthColl.at(i).PdgId()) == 11) && truthColl.at(hn_index.at(0)).IndexMother()==truthColl.at(i).IndexMother()){
      lep1_index.push_back(i);
      break;
    }
  }
  for( int i=2; i<max ; i++){
    if( (fabs(truthColl.at(i).PdgId()) == 13|| fabs(truthColl.at(i).PdgId()) == 11)  && (truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 9900012)){
      lep2_index.push_back(i);
      break;
    }
  }
  if(lep1_index.size() * lep2_index.size() == 0) return;

  snu::KTruth lep[2];
  lep[0] = truthColl.at(lep1_index.at(0));
  lep[1] = truthColl.at(lep2_index.at(0));

  if(lep[0].Pt() < lep[1].Pt()){
    snu::KTruth temp;
    temp = lep[0];
    lep[0] = lep[1];
    lep[1] = temp;
  } 

  TString s1 = "";
  if(lep[0].PdgId() * lep[1].PdgId() > 0)  s1 = "SS";
  else s1 = "OS";
/*  if(s1 == "SS"){   TruthPrintOut();
cout<<lep[0].PdgId() <<" "<<lep[1].PdgId()<<endl;}*/

  FillHist("HeavyNeutrinoMass_"+s1, truthColl.at(hn_index.at(0)).M(), 1., 0., 1000., 1000);
  FillHist("LeadingPt_"+s1, lep[0].Pt(), 1., 0., 100., 20);
  FillHist("SubLeadingPt_"+s1, lep[1].Pt(), 1., 0., 100., 20);
  FillHist("LeadingEta_"+s1, lep[0].Eta(), 1., -4., 4., 80);
  FillHist("SubLeadingEta_"+s1, lep[1].Eta(), 1., -4., 4., 80);


  return;
}// End of execute event loop
  


void SampleValidation::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void SampleValidation::BeginCycle() throw( LQError ){
  
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

SampleValidation::~SampleValidation() {
  
  Message("In SampleValidation Destructor" , INFO);
  
}


void SampleValidation::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void SampleValidation::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this SampleValidationCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void SampleValidation::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
