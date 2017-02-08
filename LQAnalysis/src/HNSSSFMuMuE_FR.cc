// $Id: HNSSSFMuMuE_FR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNSSSFMuMuE_FR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNSSSFMuMuE_FR.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNSSSFMuMuE_FR);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNSSSFMuMuE_FR::HNSSSFMuMuE_FR() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNSSSFMuMuE_FR");
  
  Message("In HNSSSFMuMuE_FR constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNSSSFMuMuE_FR::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
 
  string lqdir = getenv("LQANALYZER_DIR");

  TFile* file_single_highdxy = new TFile( (lqdir+"/data/Fake/80X/13TeV_muon_FR_SingleMuonTrigger_HighdXY.root").c_str() );
  hist_single_highdxy = (TH2F*)file_single_highdxy->Get("SingleMuonTrigger_HighdXY_NEvents_F")->Clone(); 

  TFile* file_single_dijet = new TFile( (lqdir+"/data/Fake/80X/13TeV_muon_FR_SingleMuonTrigger_Dijet.root").c_str() );
  hist_single_dijet = (TH2F*)file_single_dijet->Get("SingleMuonTrigger_Dijet_NEvents_F")->Clone();

  return;
}


void HNSSSFMuMuE_FR::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) return;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  
  FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", weight);

  /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton
  
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                                              
  float pileup_reweight=(1.0);
    
  TString mumue_trigger="HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v";
  vector<TString> trignames;
  trignames.push_back(mumue_trigger);
  float weight_trigger = WeightByTrigger("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v", TargetLumi);

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_NOCUT");
  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
  int nbjet = NBJet(GetJets("JET_HN"));

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE", false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TRI_TIGHT", false);

  std::vector<snu::KElectron> electronLooseColl = GetElectrons("ELECTRON_HN_FAKELOOSE", false);
  std::vector<snu::KElectron> electronTightColl = GetElectrons("ELECTRON_HN_TIGHT", false);

  bool trig_pass=PassTrigger("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v");
  if(!trig_pass) return;
  CorrectMuonMomentum(muonLooseColl);
   
  double ev_weight = weight;

  if( (muonLooseColl.size() != 2) || (electronLooseColl.size() != 1) ) return;
  if( (muonTightColl.size() == 2) && (electronTightColl.size() == 1) ) return;

  RAWmu[0] = muonLooseColl.at(0);
  RAWmu[1] = muonLooseColl.at(1);
  RAWel = electronLooseColl.at(0);

  if( RAWmu[0].Charge() != RAWmu[1].Charge() ) return;
  if( RAWmu[0].Charge() == RAWel.Charge() ) return;
  if( RAWmu[1].Charge() == RAWel.Charge() ) return;

  if( RAWmu[0].Pt() < 10 || RAWmu[1].Pt() < 10 || RAWel.Pt() < 10 ) return;

  cout << GetFakeRate("HighdXY", RAWmu[0]) << "  pt : " << RAWmu[0].Pt() << "  eta : " << RAWmu[0].Eta() << endl;

/*  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();

  snu::KParticle MET, RAWnu[2];
  MET.SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), 0., METPt);

  vector<double> fakerate;
  fakerate.clear();


  // HN mass divided in 4 classes
  // =============================================================================
  // =========================================
  // == CLASS 1
  // ========== 5 10 20 30 40 50
  //
  // =========================================
  // == CLASS 2
  // ========== 60 70
  //
  // =========================================
  // == CLASS 3
  // ========== 90 100 150 200
  //
  // =========================================
  // == CLASS4
  // ========== 300 400 500 700 1000
  // =============================================================================
 
  snu::KParticle W_lepton_lowmass, W_lepton_highmass;

  W_lepton_lowmass = (RAWmu[0] + RAWmu[1] + RAWel);
  W_lepton_highmass = RAWel;

  double nuPz = 999.;
  
  // =============================================================================
  // == LOW MASS REGION===========================================================
  // ====================================
  nuPz = 999.;
  nuPz = CalculateNuPz(W_lepton_lowmass, MET, 1);
  RAWnu[0].SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), nuPz, TMath::Sqrt( METPt*METPt + nuPz*nuPz ));
  nuPz = CalculateNuPz(W_lepton_lowmass, MET, -1);
  RAWnu[1].SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), nuPz, TMath::Sqrt( METPt*METPt + nuPz*nuPz ));
  if( fabs(RAWnu[0].Pz()) < fabs(RAWnu[1].Pz()) ){
    RECOnu_lowmass = RAWnu[0];
  }
  else RECOnu_lowmass = RAWnu[1];

  RECOW_pri_lowmass = RAWmu[0] + RAWmu[1] + RAWel + RECOnu_lowmass;
  RECOW_sec_lowmass = RAWel + RECOnu_lowmass;

  cout << RECOW_pri_lowmass.M() << endl;

  // ========== CLASS 1 =====================================
  EventSelectionStudy(RAWmu, RAWel, 1);// RECO particles output
  RECOHN[0] = RECOmu[1] + RECOel + RECOnu_lowmass;

  // ========== CLASS 2 =====================================
  EventSelectionStudy(RAWmu, RAWel, 2);
  RECOHN[1] = RECOmu[1] + RECOel + RECOnu_lowmass;


  // ==============================================================================
  // == HIGH MASS REGION===========================================================
  // ====================================
  nuPz = 999.;
  nuPz = CalculateNuPz(W_lepton_highmass, MET, 1);
  RAWnu[0].SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), nuPz, TMath::Sqrt( METPt*METPt + nuPz*nuPz ));
  nuPz = CalculateNuPz(W_lepton_highmass, MET, -1);
  RAWnu[1].SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), nuPz, TMath::Sqrt( METPt*METPt + nuPz*nuPz ));
  if( fabs(RAWnu[0].Pz()) < fabs(RAWnu[1].Pz()) ){
    RECOnu_highmass = RAWnu[0];
  }
  else RECOnu_highmass = RAWnu[1];

  RECOW_pri_highmass = RAWmu[0] + RAWmu[1] + RAWel + RECOnu_highmass;
  RECOW_sec_highmass = RAWel + RECOnu_highmass;

//cout << RECOW_sec_highmass.M() << endl;

  // ========== CLASS 3 =====================================
  EventSelectionStudy(RAWmu, RAWel, 3);// RECO particles output
  RECOHN[2] = RECOmu[1] + RECOel + RECOnu_highmass;

  // ========== CLASS 4 =====================================
  EventSelectionStudy(RAWmu, RAWel, 4);// RECO particles output
  RECOHN[3] = RECOmu[1] + RECOel + RECOnu_highmass;


  DrawHistograms("cut0", weight);
  FillCLHist(sssf_mumue, "cut0", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);

  ////DrawHistograms("cut", weight);
  ////FillCLHist(~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~);
  ////if ( some boolean cut )
  //////DrawHistograms("cut", weight);
*/

  return;

}// End of execute event loop
  


void HNSSSFMuMuE_FR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNSSSFMuMuE_FR::BeginCycle() throw( LQError ){
  
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

HNSSSFMuMuE_FR::~HNSSSFMuMuE_FR() {
  
  Message("In HNSSSFMuMuE_FR Destructor" , INFO);
  
}


void HNSSSFMuMuE_FR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNSSSFMuMuE_FR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNSSSFMuMuE_FRCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNSSSFMuMuE_FR::DrawHistograms(TString suffix, double weight){

  FillHist("number_of_events_"+suffix, 0., weight, 0., 1., 1);
  FillHist("W_primary_lowmass_"+suffix, RECOW_pri_lowmass.M(), weight, 0., 1000., 2000);
  FillHist("W_secondary_lowmass_"+suffix, RECOW_sec_lowmass.M(), weight, 0., 1000., 1000);
  FillHist("W_primary_highmass_"+suffix, RECOW_pri_highmass.M(), weight, 0., 1000., 2000);
  FillHist("W_secondary_highmass_"+suffix, RECOW_sec_highmass.M(), weight, 0., 1000., 1000);
  FillHist("HN_mass_class1_"+suffix, RECOHN[0].M(), weight, 0., 200., 200);
  FillHist("HN_mass_class2_"+suffix, RECOHN[1].M(), weight, 0., 200., 200);
  FillHist("HN_mass_class3_"+suffix, RECOHN[2].M(), weight, 0., 800., 800);
  FillHist("HN_mass_class4_"+suffix, RECOHN[3].M(), weight, 0., 1500., 1500);

  return;

}


void HNSSSFMuMuE_FR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
 
  RAWmu[0].SetPxPyPzE(0,0,0,0); RAWmu[1].SetPxPyPzE(0,0,0,0); RAWel.SetPxPyPzE(0,0,0,0); RAWnu[0].SetPxPyPzE(0,0,0,0); RAWnu[1].SetPxPyPzE(0,0,0,0);
  RECOmu[0].SetPxPyPzE(0,0,0,0); RECOmu[1].SetPxPyPzE(0,0,0,0); RECOel.SetPxPyPzE(0,0,0,0); RECOnu_lowmass.SetPxPyPzE(0,0,0,0); RECOnu_highmass.SetPxPyPzE(0,0,0,0); RECOW_pri_lowmass.SetPxPyPzE(0,0,0,0); RECOW_sec_lowmass.SetPxPyPzE(0,0,0,0); RECOW_pri_highmass.SetPxPyPzE(0,0,0,0); RECOW_sec_highmass.SetPxPyPzE(0,0,0,0); RECOHN[0].SetPxPyPzE(0,0,0,0); RECOHN[1].SetPxPyPzE(0,0,0,0); RECOHN[2].SetPxPyPzE(0,0,0,0); RECOHN[3].SetPxPyPzE(0,0,0,0);
 
  out_muons.clear();
  out_electrons.clear();
}


void HNSSSFMuMuE_FR::EventSelectionStudy( snu::KParticle RAWmu[], snu::KParticle RAWel, int signal_class ){

  if( signal_class == 1 || signal_class == 3 ){
    if(RAWmu[0].Pt() > RAWmu[1].Pt()){
      RECOmu[0] = RAWmu[0];
      RECOmu[1] = RAWmu[1];
    }
    else{
      RECOmu[0] = RAWmu[1];
      RECOmu[1] = RAWmu[0];
    }
  } 
  else if( signal_class == 2 || signal_class == 4 ){
    if(RAWmu[0].Pt() > RAWmu[1].Pt()){
      RECOmu[0] = RAWmu[1];
      RECOmu[1] = RAWmu[0];
    }
    else{
      RECOmu[0] = RAWmu[0];
      RECOmu[1] = RAWmu[1];
    }
  }

  RECOel = RAWel;

  return;

}


double HNSSSFMuMuE_FR::GetFakeRate( TString method, snu::KParticle ptl ){

  //                     0    1      2    3    4
  double etaarray [] = {0.0, 0.8, 1.479, 2.0, 2.5};
  //                    0   1   2   3   4   5   6   7   8    9
  //                      1   2   3   4   5   6   7   8   9
  double ptarray [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};

  double ptlPt = ptl.Pt();
  double ptlEta = fabs(ptl.Eta());

  int ptbin = -999; int etabin = -999;


  if( ptlPt > 100. ) ptbin = 9;
  else{
    for( int i = 0; i < 9 ; i++ ){
      if( (ptarray[i] < ptlPt) && (ptlPt < ptarray[i+1]) ){
        ptbin = i+1;
        break;
      }
    }
  }

  for( int i = 0; i < 4; i++ ){
    if( (etaarray[i] < ptlEta) && (ptlEta < etaarray[i+1]) ){
      etabin = i+1;
      break;
    }
  }

  double fakerate = -999.;
  if( method == "HighdXY" ){
    fakerate = hist_single_highdxy->GetBinContent(ptbin, etabin);
  }
  if( method == "Dijet" ){
    fakerate = hist_single_dijet->GetBinContent(ptbin, etabin);
  }       

  return fakerate;
}
