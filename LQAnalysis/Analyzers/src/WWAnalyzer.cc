// $Id: WWAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQWWAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "WWAnalyzer.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (WWAnalyzer);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
WWAnalyzer::WWAnalyzer() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("WWAnalyzer");
  
  Message("In WWAnalyzer constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void WWAnalyzer::InitialiseAnalysis() throw( LQError ) {
  
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


void WWAnalyzer::ExecuteEvents()throw( LQError ){

//TruthPrintOut();
  

std::vector<snu::KTruth> truthColl;
eventbase->GetTruthSel()->Selection(truthColl);

int max = truthColl.size();


  std::vector<int> lep_index, neu_index , q1_index , q2_index; 


for( int i=3 ; i<max ; i++){
   if( (abs(truthColl.at(i).PdgId()) == 13 || abs(truthColl.at(i).PdgId()) == 11) && abs(truthColl.at( truthColl.at(i).IndexMother()).PdgId()) == 24){
       lep_index.push_back(i);
       break;
   }
  }


for( int i=3 ; i<max ; i++){
   if( (abs(truthColl.at(i).PdgId()) == 14 || abs(truthColl.at(i).PdgId()) == 12) && abs(truthColl.at( truthColl.at(i).IndexMother() ).PdgId()) == 24){
         neu_index.push_back(i); 
         break;
   }
  }

for( int i=3 ; i<max ; i++){
   if( abs(truthColl.at(i).PdgId()) < 6 && abs(truthColl.at( truthColl.at(i).IndexMother()).PdgId()) == 24){
       q1_index.push_back(i);
       break;
   }
  }


for( int i=3 ; i<max ; i++){
   if( abs(truthColl.at(i).PdgId()) < 6 && abs(truthColl.at( truthColl.at(i).IndexMother() ).PdgId()) == 24){   
       if( q1_index.at(0) != i){
         q2_index.push_back(i);
         break;
       }
   }
  }
//  cout << "ael_index :::::" << ael_index << endl;


  bool charge_is_minus = false;

  bool mc_found = false;

  if( lep_index.size() * neu_index.size() * q1_index.size() * q2_index.size() != 0) { 
     mc_found = true;
     }
 
  if(!mc_found) return;

  snu::KTruth lep, neu, q1, q2;
  lep = truthColl.at(lep_index.at(0));
  neu = truthColl.at(neu_index.at(0));
  q1 = truthColl.at(q1_index.at(0));
  q2 = truthColl.at(q2_index.at(0));

  if(lep.PdgId() > 0 ) charge_is_minus = true;
  FillHist("mass_of_lv", (lep+neu).M(), 1., 0., 200, 40); //histname variable weight=1 xmin xmax nbins
  FillHist("pt_of_q", q1.Pt(), 1., 0., 200, 40);
  FillHist("pt_of_q", q2.Pt(), 1., 0., 200, 40);
  FillHist("pt_of_q1", q1.Pt(), 1., 0., 200, 40);
  FillHist("pt_of_q2", q2.Pt(), 1., 0., 200, 40);
  FillHist("mass_of_q1q2", (q1+q2).M(), 1., 0., 200, 40);
  FillHist("mass_of_qqlv", (lep+neu+q1+q2).M(), 1., 0., 500, 50);
  FillHist("deltaR_of_w+_and_w-", (lep+neu).DeltaR(q1+q2), 1., 0., 5., 25);


 // pt of q1 q2 // FillHist("pt_of_q", q1.Pt(),  // FillHist("pt_of_q", q2.Pt(),
 // pt of q1
 // pt of q2
 // pt of neutrino  
 // pt of lepton
 // mass of lep neu
 // mass of q1 q2
 // mass of lep neu q1 q2
 // deltaR of w+ and w- FillHist( , (lep+neu).DeltaR(q1+q2), 1., 0., 5., 25);
 // 
  if(!charge_is_minus) FillHist("charge_of_w", 1., 1., -2., 3., 5); //!!!!!!!!
  else FillHist("charge_of_w", -1., 1., -2., 3., 5);


   
return;
}// End of execute event loop
  


void WWAnalyzer::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void WWAnalyzer::BeginCycle() throw( LQError ){
  
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

WWAnalyzer::~WWAnalyzer() {
  
  Message("In WWAnalyzer Destructor" , INFO);
  
}


void WWAnalyzer::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void WWAnalyzer::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this WWAnalyzerCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void WWAnalyzer::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



