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
  
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void ExampleAnalyzer::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  

  return;

}


void ExampleAnalyzer::ExecuteEvents()throw( LQError ){

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                                              
  std::vector<TString> triggerlist_mm_All; triggerlist_mm_All.clear();

  triggerlist_mm_All.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  triggerlist_mm_All.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  triggerlist_mm_All.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_mm_All.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  if(!PassTriggerOR(triggerlist_mm_All)) return;

  std::vector<snu::KMuon> loosemuons, muons; loosemuons.clear(); muons.clear();
  loosemuons = GetMuons("MUON_POG_LOOSE", true);

  for(unsigned int i=0; i<loosemuons.size(); i++){
    if(PassID(loosemuons.at(i), "MUON_POG_TIGHT")){
      muons.push_back(loosemuons.at(i));
    }
  }

  if(muons.size() != 2) return;

  CorrectMuonMomentum(muons);
  if((muons.at(0).Pt() < 20) || (muons.at(1).Pt() < 10)) return;
  if(muons.at(0).Charge() == muons.at(1).Charge()) return;

  CorrectedMETRochester(muons);

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  snu::KParticle MET;
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);

  std::vector<snu::KElectron> vetoelectrons = GetElectrons(true, true, "ELECTRON_POG_TIGHT");

  double this_weight = 1.;
  if(!isData){
    this_weight = weight;
    this_weight *= MCweight;
    this_weight *= mcdata_correction->MuonScaleFactor("MUON_POG_TIGHT", muons, 0);
    this_weight *= mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muons, 0);
    this_weight *= mcdata_correction->MuonTrackingEffScaleFactor(muons);
    this_weight *= mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    this_weight *= GetKFactor();
  }

  DrawHistograms(this_weight, muons, MET, "preselection");
  if(MET.Pt() < 100) DrawHistograms(this_weight, muons, MET, "met50cut");
  if(vetoelectrons.size() == 0) DrawHistograms(this_weight, muons, MET, "noelectrons");
  if(fabs((muons.at(0) + muons.at(1)).M() - 91.1876) < 20){
    DrawHistograms(this_weight, muons, MET, "zpeak");
    if(vetoelectrons.size() == 0) DrawHistograms(this_weight, muons, MET, "zpeak_noelectrons");
  }

  return;
}// End of execute event loop
  
void ExampleAnalyzer::DrawHistograms(double this_weight, std::vector<snu::KMuon> muons, snu::KParticle MET, TString prefix){

  FillHist(prefix + "_mu1_pt", muons.at(0).Pt(), this_weight, 0., 200., 200);
  FillHist(prefix + "_mu1_eta", muons.at(0).Eta(), this_weight, -3., 3., 60);
  FillHist(prefix + "_mu1_dxy", muons.at(0).dXY(), this_weight, -1., 1., 200);
  FillHist(prefix + "_mu1_dz", muons.at(0).dZ(), this_weight, -1., 1., 200);
  FillHist(prefix + "_mu1_reliso", muons.at(0).RelIso04(), this_weight, 0., 1., 200);

  FillHist(prefix + "_mu2_pt", muons.at(1).Pt(), this_weight, 0., 200., 200);
  FillHist(prefix + "_mu2_eta", muons.at(1).Eta(), this_weight, -3., 3., 60);
  FillHist(prefix + "_mu2_dxy", muons.at(1).dXY(), this_weight, -1., 1., 200);
  FillHist(prefix + "_mu2_dz", muons.at(1).dZ(), this_weight, -1., 1., 200);
  FillHist(prefix + "_mu2_reliso", muons.at(1).RelIso04(), this_weight, 0., 1., 200);

  FillHist(prefix + "_mumu_mass", (muons.at(0)+muons.at(1)).M(), this_weight, 0., 500., 500);
  FillHist(prefix + "_mumu_pt", (muons.at(0)+muons.at(1)).Pt(), this_weight, 0., 200., 200);

  FillHist(prefix + "_met", MET.Pt(), this_weight, 0., 200., 200);

  return;
}

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



