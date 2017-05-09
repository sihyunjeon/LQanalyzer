// $Id: CFRateCalculator_FR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQCFRateCalculator_FR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "CFRateCalculator_FR.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (CFRateCalculator_FR);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
CFRateCalculator_FR::CFRateCalculator_FR() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("CFRateCalculator_FR");
  
  Message("In CFRateCalculator_FR constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void CFRateCalculator_FR::InitialiseAnalysis() throw( LQError ) {
  
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


void CFRateCalculator_FR::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  if(!PassMETFilter()) return;     /// Initial event cuts :
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex

  CFvalidation();

  return;
}// End of execute event loop
  


void CFRateCalculator_FR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void CFRateCalculator_FR::BeginCycle() throw( LQError ){
  
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

CFRateCalculator_FR::~CFRateCalculator_FR() {
  
  Message("In CFRateCalculator_FR Destructor" , INFO);
  
}


void CFRateCalculator_FR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void CFRateCalculator_FR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this CFRateCalculator_FRCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void CFRateCalculator_FR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


void CFRateCalculator_FR::GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index ){

  for( int i = it+1 ; i < truthColl.size(); i ++ ){
    if( truthColl.at(i).IndexMother() == it && truthColl.at(i).PdgId() == truthColl.at(it).PdgId() ){
      index.push_back(i);
    }
  }
  return;

}



void CFRateCalculator_FR::CFvalidation(void){

  bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if( !pass_trig ) return;

  double Z_mass = 91.1876;

  // define electron and muon Colls
  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false, false, "ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false, false, "ELECTRON_HN_TIGHT");

  if(electronLooseColl.size() != 2) return;
  if(isData){ if(electronTightColl.size() == 2) return;}
  else{ if(electronTightColl.size() != 2) return;}

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_VETO", false);
  if( muonLooseColl.size() != 0) return;

  FillHist("[CHECK]n_of_muons", muonLooseColl.size(), 1., 0., 5., 5);//check no muons
  FillHist("[CHECK]n_of_electrons", electronTightColl.size(), 1., 0., 5., 5);//check two electrons

  // define leptons and give Pt, MET cuts
  snu::KParticle lep[2];
  lep[0] = electronLooseColl.at(0);
  lep[1] = electronLooseColl.at(1);

  if( (lep[0].Charge() != lep[1].Charge()) ) return;

  if( lep[0].Pt() < 25 || lep[1].Pt() < 25 ) return;

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  if(METPt > 30) return;

  FillHist("[CHECK]MET", METPt, 1., 0., 100., 100);

  snu::KParticle Z_candidate_SHIFT0;//define Z candidate before energy shift
  Z_candidate_SHIFT0 = (lep[0] + lep[1]);
  bool Z_selection_SHIFT0 = (fabs(Z_candidate_SHIFT0.M() - Z_mass) < 10.);

  for(int i=0; i<2; i++){
    lep[i] = ShiftEnergy(lep[i], 1./0.981);
    FillHist("[CHECK]lepton_E_shift_up", lep[i].E(), 1., 0., 400., 400);
  }

  snu::KParticle Z_candidate;//define Z candidate after energy shift
  Z_candidate = (lep[0] + lep[1]);
  bool Z_selection = (fabs(Z_candidate.M() - Z_mass) < 10.);

  if( !Z_selection_SHIFT0 && !Z_selection ) return;

  double fake_weight = -999., fake_weight_err = -999.;

  if(isData){
    fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 0, electronLooseColl, "ELECTRON_HN_TIGHT", 2);
//    fake_weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 0, electronLooseColl, "ELECTRON_HN_TIGHT", 2);
  }
  else{
    fake_weight = weight*MCweight;
//    fake_weight_err = 0.;
  }

  bool is_region[2][2] = {{false,},};
  if( (fabs(lep[0].Eta()) < 1.4442) )                                   is_region[0][0] = true;
  else if( (fabs(lep[0].Eta()) > 1.556) && (fabs(lep[0].Eta()) < 2.5) ) is_region[1][0] = true;
  if( (fabs(lep[1].Eta()) < 1.4442) )                                   is_region[0][1] = true;
  else if( (fabs(lep[1].Eta()) > 1.556) && (fabs(lep[1].Eta()) < 2.5) ) is_region[1][1] = true;

  TString region = "none";
  if((is_region[0][0] && is_region[0][1]))      region = "BB";
  else if((is_region[1][0] && is_region[1][1])) region = "EE";
  else                                          region = "BE";

  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN", 30., 2.4);
  TString njets = "";
  if( jetTightColl.size() == 0 ) njets = "JETS0_";
  if( jetTightColl.size() > 1 ) njets = "JETS_";

  if( Z_selection_SHIFT0 ){//Z selection before shifting down energy
    FillHist("observed_Z_mass_global", Z_candidate_SHIFT0.M(), fake_weight, 70., 110., 40);
    FillHist("observed_n_events_global", 0., fake_weight, 0., 1., 1);
    FillHist("observed_Z_mass_"+region, Z_candidate_SHIFT0.M(), fake_weight, 70., 110., 40);
    FillHist("observed_n_events_"+region, 0., fake_weight, 0., 1., 1);
    FillHist(njets+"observed_Z_mass_global", Z_candidate_SHIFT0.M(), fake_weight, 70., 110., 40);
    FillHist(njets+"observed_n_events_global", 0., fake_weight, 0., 1., 1);
  }

  if( Z_selection ){//Z selection after shifting down energy
    FillHist("observed_Z_mass_SHIFT_global", Z_candidate.M(), fake_weight, 70., 110., 40);
    FillHist("observed_n_events_SHIFT_global", 0., fake_weight, 0., 1., 1);
    FillHist("observed_Z_mass_SHIFT_"+region, Z_candidate.M(), fake_weight, 70., 110., 40);
    FillHist("observed_n_events_SHIFT_"+region, 0., fake_weight, 0., 1., 1);
    FillHist(njets+"observed_Z_mass_SHIFT_global", Z_candidate.M(), fake_weight, 70., 110., 40);
    FillHist(njets+"observed_n_events_SHIFT_global", 0., fake_weight, 0., 1., 1);
  }

  return;
}

snu::KParticle CFRateCalculator_FR::ShiftEnergy( snu::KParticle old_lep, double shift_rate ){

   double new_E, new_px, new_py, new_pz;
   double mass;
   double new_psum, old_psum;
   new_E = old_lep.E() * shift_rate;
   mass = 0.511e-3;

   new_psum = sqrt(new_E*new_E - mass*mass);  
   old_psum = sqrt(old_lep.Pt()*old_lep.Pt() + old_lep.Pz()*old_lep.Pz());

   double ratio = -999.;
   ratio = new_psum/old_psum;

   new_px = old_lep.Px() * ratio;
   new_py = old_lep.Py() * ratio;
   new_pz = old_lep.Pz() * ratio;

   snu::KParticle new_lep;
   new_lep.SetPxPyPzE(new_px,new_py,new_pz,new_E);

   return new_lep;
}

