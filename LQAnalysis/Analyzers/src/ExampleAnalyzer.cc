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

  if(!PassMETFilter()) return;     /// Initial event cuts : 

//  bool trig_pass = (PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v") || PassTrigger("HLT_Ele17_CaloIdL_GsfTrkIdVL_v"));
  bool trig_pass =PassTrigger("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v");
  if(!trig_pass ) return;

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                
  std::vector<snu::KElectron> electronTightColl, electronPromptColl;
  electronPromptColl.clear();
  bool dxy001 = false;
  if(std::find(k_flags.begin(), k_flags.end(), "p1") !=k_flags.end()) dxy001 = true;
  bool dxy002 = false;
  if(std::find(k_flags.begin(), k_flags.end(), "p2") !=k_flags.end()) dxy002 = true;

  if(dxy001 && !dxy002) electronTightColl =  GetElectrons(true,false,"ELECTRON_HN_LOWDXY_TIGHT_001");
  else if(!dxy001 && dxy002) electronTightColl =  GetElectrons(true,false,"ELECTRON_HN_LOWDXY_TIGHT_002");
  else {FillHist("[ERROR]flag_not_defined", 0., 1., 0., 1., 1); return;}

  if(electronTightColl.size() != 2) return;

  if(electronTightColl.at(0).Charge() != electronTightColl.at(1).Charge() ) return;

  cout<<"have two SS electrons................................."<<endl;

  snu::KParticle RECOel[2];
  RECOel[0] = electronTightColl.at(0);
  RECOel[1] = electronTightColl.at(1);

  if(RECOel[0].Pt() <100) return;
  if(RECOel[1].Pt() <20) return;

  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  vector<int> el_index, anti_el_index;
  el_index.clear(); anti_el_index.clear();

  for( int i = 2 ; i < truthColl.size() ; i++ ){
    if( (truthColl.at(i).PdgId()) == 11 ){
      el_index.push_back(i);
      GENFindDecayIndex( truthColl, i, el_index);
      break;
    }
  }
  for( int i = 2 ; i < truthColl.size() ; i++ ){
    if( (truthColl.at(i).PdgId()) == -11 ){
      anti_el_index.push_back(i);
      GENFindDecayIndex( truthColl, i, anti_el_index);
      break;
    }
  }

  snu::KTruth GENel[2];
  if(el_index.size() == 0 || anti_el_index.size() == 0) return;

  if( truthColl.at(el_index.back()).Pt() > truthColl.at(anti_el_index.back()).Pt() ){
    GENel[0] = truthColl.at(el_index.back());
    GENel[1] = truthColl.at(anti_el_index.back());
  }
  else{
    GENel[1] = truthColl.at(el_index.back());
    GENel[0] = truthColl.at(anti_el_index.back());
  }

  TruthPrintOut();

  cout<<"RECO============================================================"<<endl;
  cout<<" Chagre :: "<<RECOel[0].Charge() << "  " << RECOel[1].Charge()<<endl;
  cout<<"  leading Pt  "<<endl;
  cout<<"     "<<RECOel[0].Pt()<<endl;
  cout<<"  subleading Pt  "<<endl;
  cout<<"     "<<RECOel[1].Pt()<<endl;
  cout<<"  leading Eta  "<<endl;
  cout<<"     "<<RECOel[0].Eta()<<endl;
  cout<<"  subleading Eta  "<<endl;
  cout<<"     "<<RECOel[1].Eta()<<endl;
  cout<<"  leading Phi  "<<endl;
  cout<<"     "<<RECOel[0].Phi()<<endl;
  cout<<"  subleading Phi  "<<endl;
  cout<<"     "<<RECOel[1].Phi()<<endl;
  cout<<"================================================================="<<endl;


 /*

  bool is_CF = false;

  for(int i=0; i<electronPromptColl.size(); i++){

    is_CF = false;

    snu::KElectron this_lep;
    this_lep = electronPromptColl.at(i);

    if( (this_lep.MCIsCF()) ) is_CF = true;

    FillHist("PROMPT_PT", this_lep.Pt(), 1., 200., 1000., 800);

    if( is_CF ){

cout<<"################################################################################################"<<endl;
cout<<"event # : " << eventbase->GetEvent().EventNumber()<< endl;
cout<<"Pt : " << this_lep.Pt() << endl;
cout<<"Eta : "<< this_lep.Eta() << endl;

      FillHist("PROMPT_CF_PT", this_lep.Pt(), 1., 200., 1000., 800);
      FillHist("PROMPT_CF_ETA", this_lep.Eta(), 1., -3., 3., 60);

    }


  }
*/
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
