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
  
  Message("Making clever hists for Z ->ll test code", INFO);
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:


  return;
}


void HNSSSFMuMuE_FR::ExecuteEvents()throw( LQError ){

  // ========== Apply the gen weight ====================
  if(!isData) weight*=MCweight;
  // ================================================================================
  

  // ========== No cut ====================
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  FillCutFlow("NoCut", weight);
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  // ================================================================================


  // ========== MET filter cut ====================
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", weight);
  // ================================================================================


  // ========== Primary vertex cut ====================
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                                              
  // ================================================================================


  // ========== Trigger cut ====================
  TString mumu_trigger="HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v";
  vector<TString> trignames;
  trignames.push_back(mumu_trigger);

  bool trig_pass=PassTrigger("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v");
  if(!trig_pass) return;
  // ================================================================================


  // ========== Get Objects (muon, electron, jet) ====================
  bool isKeepFake = std::find(k_flags.begin(), k_flags.end(), "keepfake") != k_flags.end();

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TRI_TIGHT",false);

  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false,false,"ELECTRON16_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,false,"ELECTRON16_HN_TIGHT");

  if(isKeepFake){
    muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",true);
    muonTightColl = GetMuons("MUON_HN_TRI_TIGHT",true);

    electronLooseColl = GetElectrons(true,true,"ELECTRON16_HN_FAKELOOSE");
    electronTightColl = GetElectrons(true,true,"ELECTRON16_HN_TIGHT");
  }

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_NOCUT");
  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
  // ================================================================================


  // ========== Rochester Correction ====================
  CorrectMuonMomentum(muonLooseColl);

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();

  METPt = CorrectedMETRochester(muonLooseColl, METPt, METPhi, true);
  METPhi = CorrectedMETRochester(muonLooseColl, METPt, METPhi, false);

  MET.SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), 0., METPt);
  // ================================================================================


  // ========== Pileup reweight ====================
  float pileup_reweight=(1.0);
  if(!k_isdata){ pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);} 
  // ================================================================================


  // ========== Trigger reweight ====================
  float weight_trigger = WeightByTrigger("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v", TargetLumi);
  // ================================================================================


  // ========== Muon tracking efficiency ====================  
  float muon_trkeff = mcdata_correction->MuonTrackingEffScaleFactor(muonLooseColl);
  // ================================================================================


  // ========== Electron ID scalefactor ====================
  float electron_idsf = mcdata_correction->ElectronScaleFactor("ELECTRON16_HN_FAKELOOSE", electronLooseColl);
  // ================================================================================
  

  // ========== Electron RECO scalefactor ====================
  float electron_recosf = mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);
  // ================================================================================


  // ========== Reweight ====================
  if(!isData){
    weight *= weight_trigger;
    weight *= pileup_reweight;
    weight *= muon_trkeff;
    weight *= electron_idsf;
    weight *= electron_recosf;
  }
  // ================================================================================


  // ========== Others ====================
  int period_index = 0;
  period_index = GetPeriodIndex();
  nbjet = NBJet(jetLooseColl, snu::KJet::CSVv2, snu::KJet::Medium, period_index);
  // ================================================================================





  /*####################################################################################################
  ##		        Analysis Code 								      ##
  ##				For SameSign MuMuE Channel Analysis				      ##
  ####################################################################################################*/

  if( !((muonLooseColl.size() == 2) && (electronLooseColl.size() == 1)) ) return;
  if( ((muonTightColl.size() == 2) && (electronTightColl.size() == 1)) ) return;

  if(isData){
    weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON16_HN_TIGHT", 1);
    weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON16_HN_TIGHT", 1);
  }

  RAWmu[0] = muonLooseColl.at(0);
  RAWmu[1] = muonLooseColl.at(1);
  RAWel = electronLooseColl.at(0);

  if( RAWmu[0].Charge() != RAWmu[1].Charge() ) return;
  if( RAWmu[0].Charge() == RAWel.Charge() ) return;
  if( RAWmu[1].Charge() == RAWel.Charge() ) return;

  if( RAWmu[0].Pt() < 10 || RAWmu[1].Pt() < 10 || RAWel.Pt() < 10 ) return;

  if( ((RAWmu[0]+RAWmu[1]).M() < 4) || ((RAWmu[0]+RAWel).M() < 4) || ((RAWmu[1]+RAWel).M() < 4) ) return;


  // ================================================================================
  // HN mass divided in 4 classes
  // ================================================================================
  // ====== CLASS 1
  // ============== 5 10 20 30 40 50
  // ====== CLASS 2
  // ============== 60 70
  // ====== CLASS 3
  // ============== 90 100 150 200
  // ====== CLASS 4
  // ============== 300 400 500 700 1000
  // ================================================================================
 
  snu::KParticle W_lepton_lowmass, W_lepton_highmass;
  W_lepton_lowmass = (RAWmu[0] + RAWmu[1] + RAWel);
  W_lepton_highmass = RAWel;

  double nuPz;
 
  // ================================================================================
  // ====== LOW MASS REGION
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

  // ========== CLASS 1 =====================================
  EventSelectionStudy(RAWmu, RAWel, 1);// RECO particles output
  RECOHN[0] = RECOmu[1] + RECOel + RECOnu_lowmass;

  // ========== CLASS 2 =====================================
  EventSelectionStudy(RAWmu, RAWel, 2);// RECO particles output
  RECOHN[1] = RECOmu[1] + RECOel + RECOnu_lowmass;


  // ================================================================================
  // ====== HIGH MASS REGION
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

  // ========== CLASS 3 =====================================
  EventSelectionStudy(RAWmu, RAWel, 3);// RECO particles output
  RECOHN[2] = RECOmu[1] + RECOel + RECOnu_highmass;

  // ========== CLASS 4 =====================================
  EventSelectionStudy(RAWmu, RAWel, 4);// RECO particles output
  RECOHN[3] = RECOmu[1] + RECOel + RECOnu_highmass;





  DrawHistograms("cut0", weight);
  FillCLHist(sssf_mumue, "cut0", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);
  FillCLHist(sssf_mumue, "cut0_up", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight+weight_err);
  FillCLHist(sssf_mumue, "cut0_down", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight-weight_err);

  // Low mass region cuts
  if( RECOW_pri_lowmass.M() < 150. ){
    DrawHistograms("cutW150", weight);
    FillCLHist(sssf_mumue, "cutW150", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);
    FillCLHist(sssf_mumue, "cutW150_up", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight+weight_err);
    FillCLHist(sssf_mumue, "cutW150_down", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight-weight_err);
  }

  // High mass region cuts
  if( METPt > 20. ){
    DrawHistograms("cutMET20", weight);
    FillCLHist(sssf_mumue, "cutMET20", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);
    FillCLHist(sssf_mumue, "cutMET20_up", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight+weight_err);
    FillCLHist(sssf_mumue, "cutMET20_down", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight-weight_err);
  }
  //virtual W mass cut can also be used



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

  FillUpDownHist("number_of_events_"+suffix, 0., weight, weight_err, 0., 1., 1);
  FillUpDownHist("W_primary_lowmass_"+suffix, RECOW_pri_lowmass.M(), weight, weight_err, 0., 1000., 2000);
  FillUpDownHist("W_secondary_lowmass_"+suffix, RECOW_sec_lowmass.M(), weight, weight_err, 0., 1000., 1000);
  FillUpDownHist("W_primary_highmass_"+suffix, RECOW_pri_highmass.M(), weight, weight_err, 0., 1000., 2000);
  FillUpDownHist("W_secondary_highmass_"+suffix, RECOW_sec_highmass.M(), weight, weight_err, 0., 1000., 1000);
  FillUpDownHist("HN_mass_class1_"+suffix, RECOHN[0].M(), weight, weight_err, 0., 500., 500);
  FillUpDownHist("HN_mass_class2_"+suffix, RECOHN[1].M(), weight, weight_err, 0., 500., 500);
  FillUpDownHist("HN_mass_class3_"+suffix, RECOHN[2].M(), weight, weight_err, 0., 800., 800);
  FillUpDownHist("HN_mass_class4_"+suffix, RECOHN[3].M(), weight, weight_err, 0., 1500., 1500);
  FillUpDownHist("NBjets_"+suffix, nbjet, weight, weight_err, 0., 5., 5);
  FillUpDownHist("[SignalStudy]deltaR_elMET_"+suffix, RECOel.DeltaR(MET), weight, weight_err, 0., 5., 100);
  FillUpDownHist("[SignalStudy]transversemass_elMET_"+suffix, MT(MET,RECOel), weight, weight_err, 0., 1000., 2000);

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
  MET.SetPxPyPzE(0,0,0,0);
  nbjet = 0; 

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


int HNSSSFMuMuE_FR::DefineClass(){

  if( k_sample_name.Contains("1000")
    || k_sample_name.Contains("700") 
    || k_sample_name.Contains("500") 
    || k_sample_name.Contains("400")
    || k_sample_name.Contains("300") ) return 4;


  else if( k_sample_name.Contains("200")
    || k_sample_name.Contains("150")
    || k_sample_name.Contains("100")
    || k_sample_name.Contains("90") ) return 3;

  else if( k_sample_name.Contains("70")
    || k_sample_name.Contains("60") ) return 2;

  else if( k_sample_name.Contains("50")
    || k_sample_name.Contains("40")
    || k_sample_name.Contains("30")
    || k_sample_name.Contains("20")
    || k_sample_name.Contains("10")
    || k_sample_name.Contains("5") ) return 1;

  else return 0;

}


int HNSSSFMuMuE_FR::GetPeriodIndex(void){
  if( isData ){
    if( k_sample_name.Contains("B") ||  k_sample_name.Contains("C") || k_sample_name.Contains("D") || k_sample_name.Contains("E") || k_sample_name.Contains("F") ){
      return 1;
    }
    else return 7;
  }
  else return 0;
}

