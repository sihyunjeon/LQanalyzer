// $Id: PRCalculatgnames.push_back(dor.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQPRCalculator Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "PRCalculator.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (PRCalculator);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
PRCalculator::PRCalculator() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("PRCalculator");
  
  Message("In PRCalculator constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void PRCalculator::InitialiseAnalysis() throw( LQError ) {
  
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


void PRCalculator::ExecuteEvents()throw( LQError ){

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
     
   
  if(!PassTrigger("HLT_Ele27_WPTight_Gsf_v")) return;   
  std::vector<snu::KElectron> electronTHNMVA = GetElectrons(false,false, "ELECTRON_HN_TIGHT");
  std::vector<snu::KElectron> electronLHNMVA = GetElectrons(false,false, "ELECTRON_HN_FAKELOOSE");

  if(!(electronLHNMVA.size() != 2)) return;
  if(!(electronTHNMVA.size() == 0)) return;

  if(electronLHNMVA.at(0).Pt() < 30) return;
  if(electronTHNMVA.at(0).Pt() < 30) return;

  double Z_mass = 91.1876;

  float ptarray[8] = {10., 20., 40., 60., 80., 100., 200., 500.};
  float etaarray[5] = {0.0, 0.9, 1.4442, 1.556, 2.5};

  snu::KParticle elTAG = electronTHNMVA.at(0);
  for(int i=0; i<electronLHNMVA.size(); i++){
      snu::KParticle elPROBE = electronLHNMVA.at(i);
      if(elPROBE.Charge()!=elTAG.Charge()){
        if(fabs((elPROBE+elTAG).M() - Z_mass) < 10){
          FillHist("den_Pt_eta_global", fabs(elPROBE.Eta()), elPROBE.Pt(), 1, etaarray, 4, ptarray, 7);
          FillHist("den_Z_mass_global", (elPROBE+elTAG).M(), 1, 60., 120., 30);
          if(PassID(electronLHNMVA.at(i), "ELECTRON_HN_TIGHT")){
            FillHist("num_Pt_eta_global", fabs(elPROBE.Eta()), elPROBE.Pt(), 1, etaarray, 4, ptarray, 7);
            FillHist("num_Z_mass_global", (elPROBE+elTAG).M(), 1, 60., 120., 30);
cout<< "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " <<endl;
          }
        }
      }
  }

  

//  if(electronGENT.size() == 1) FillHist("N_ELECTRON_TIGHT", 2., 1., 0., 3., 3);

   //   std::vector<snu::KElectron> electrons2 =  GetElectrons(BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0);

	    
   
   return;
}// End of execute event loop
  


void PRCalculator::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void PRCalculator::BeginCycle() throw( LQError ){
  
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

PRCalculator::~PRCalculator() {
  
  Message("In PRCalculator Destructor" , INFO);
  
}


void PRCalculator::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void PRCalculator::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this PRCalculatorCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void PRCalculator::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



