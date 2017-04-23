// $Id: HNSSSFMuMuE_FR_ntuple.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNSSSFMuMuE_FR_ntuple Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNSSSFMuMuE_FR_ntuple.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNSSSFMuMuE_FR_ntuple);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNSSSFMuMuE_FR_ntuple::HNSSSFMuMuE_FR_ntuple() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNSSSFMuMuE_FR_ntuple");
  
  Message("In HNSSSFMuMuE_FR_ntuple constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNSSSFMuMuE_FR_ntuple::InitialiseAnalysis() throw( LQError ) {
  
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


void HNSSSFMuMuE_FR_ntuple::ExecuteEvents()throw( LQError ){

  // ========== Define RelIso ====================
  double this_reliso_mu = 0.4;
  double this_reliso_el = 0.5; 
  // ================================================================================


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
  std::vector<TString> triggerlist;
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  if(!PassTriggerOR(triggerlist)) return;
  // ================================================================================


  // ========== Get Objects (muon, electron, jet) ====================
  std::vector<snu::KMuon> muonVLooseColl = GetMuons("MUON_HN_TRI_LOOSE_lowestPtCut",false);

  std::vector<snu::KElectron> electronVLooseColl = GetElectrons(false,false,"ELECTRON16_POG_FAKELOOSE_CC_d0_lowestPtCut");

  std::vector<snu::KJet> jetVLooseColl = GetJets("JET_HN", 20., 2.4);
  // ================================================================================


  // ========== Rochester Correction ====================
//  CorrectMuonMomentum(muonLooseColl);

  snu::KEvent event = eventbase->GetEvent();

  double METPt = event.MET();
  double METPhi = event.METPhi();
  // ================================================================================


  // ========== Pileup reweight ====================
  float pileup_reweight=(1.0), pileup_reweight_up=(1.0), pileup_reweight_down=(1.0);
  if(!k_isdata){
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    pileup_reweight_up = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),-1);
    pileup_reweight_down = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),+1);
  } 
  // ================================================================================


  // ========== Trigger reweight ====================
  float weight_trigger = WeightByTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", TargetLumi);
  // ================================================================================


  // ========== Reweight ====================
  if(!isData){
    weight *= weight_trigger;
//    weight *= pileup_reweight;
//    weight *= muon_trkeff;
//    weight *= electron_idsf;
//    weight *= electron_reco;
  }
  // ================================================================================


  int N_sys = (2*7+1);
  for(int it_sys = 0; it_sys<N_sys; it_sys++){//N_sys; it_sys++){

    double this_weight = weight;
    TString this_syst;

    if(it_sys==0){
      this_syst = "MuonEn_up";////MET = event.PFMETShifted(snu::KEvent::MuonEn, snu::KEvent::up);
    }
    else if(it_sys==1){
      this_syst = "MuonEn_down";////MET = event.PFMETShifted(snu::KEvent::MuonEn, snu::KEvent::down);
    }
    else if(it_sys==2){
      this_syst = "ElectronEn_up";//
      METPt = event.PFMETShifted(snu::KEvent::ElectronEn, snu::KEvent::up);
    }
    else if(it_sys==3){
      this_syst = "ElectronEn_down";//
      METPt = event.PFMETShifted(snu::KEvent::ElectronEn, snu::KEvent::down);
    }
    else if(it_sys==4){
      this_syst = "JetEn_up";//
      METPt = event.PFMETShifted(snu::KEvent::JetEn, snu::KEvent::up);
    }
    else if(it_sys==5){
      this_syst = "JetEn_down";//
      METPt = event.PFMETShifted(snu::KEvent::JetEn, snu::KEvent::down);
    }
    else if(it_sys==6){
      this_syst = "JetRes_up";//
//      METPt = event.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::up);
    }
    else if(it_sys==7){
      this_syst = "JetRes_down";//
//      METPt = event.PFMETShifted(snu::KEvent::JetRes, snu::KEvent::down);
    }
    else if(it_sys==8){//
      this_syst = "Unclustered_up";
      METPt = event.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::up);
    }
    else if(it_sys==9){//
      this_syst = "Unclustered_down";
      METPt = event.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::down);
    }
    else if(it_sys==10){//
      this_syst = "Central";
    }
    else if(it_sys==11){//
      this_syst = "MuonIDSF_up";
    }
    else if(it_sys==12){//
      this_syst = "MuonIDSF_down";
    }
/*    else if(it_sys==13){
      this_syst = "MuonTrkEff_up";
    }
    else if(it_sys==14){
      this_syst = "MuonTrkEff_down";
    }
    else if(it_sys==15){//FIXME
      this_syst = "ElectronIDSF_up";
    }
    else if(it_sys==16){//FIXME
      this_syst = "ElectronIDSF_down";
    }
    else if(it_sys==17){//FIXME
      this_syst = "ElectronRecoSF_up";
    }
    else if(it_sys==18){//FIXME
      this_syst = "ElectronRecoSF_down";
    }*/
    else if(it_sys==13){
      this_syst = "PU_up";
    }
    else if(it_sys==14){
      this_syst = "PU_down";
    }
    else{
      Message("it_sys out of range!" , INFO);
      return;
    }

    // ========== Muon Energy systematics ====================
    std::vector<snu::KMuon> muonLooseColl, muonTightColl;
    if(this_syst == "MuonEn_up"){
      for(unsigned int i=0; i<muonVLooseColl.size(); i++){
        snu::KMuon this_muon = muonVLooseColl.at(i);
        this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedUp(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
        double new_reliso = this_muon.RelIso04()/this_muon.PtShiftedUp();
        if( this_muon.Pt() >= 10. && new_reliso < this_reliso_mu ) muonLooseColl.push_back( this_muon );
      }
    }
    else if(this_syst == "MuonEn_down"){
      for(unsigned int i=0; i<muonVLooseColl.size(); i++){
        snu::KMuon this_muon = muonVLooseColl.at(i);
        this_muon.SetPtEtaPhiM( this_muon.Pt()*this_muon.PtShiftedDown(), this_muon.Eta(), this_muon.Phi(), this_muon.M() );
        double new_reliso = this_muon.RelIso04()/this_muon.PtShiftedDown();
        if( this_muon.Pt() >= 10. && new_reliso < this_reliso_mu ) muonLooseColl.push_back( this_muon );
      }
    }
    else{
      for(unsigned int i=0; i<muonVLooseColl.size(); i++){
        snu::KMuon this_muon = muonVLooseColl.at(i);
        if( this_muon.Pt() >= 10. && this_muon.RelIso04() < this_reliso_mu ) muonLooseColl.push_back( this_muon );
      }
    }
    for(unsigned int i=0; i<muonLooseColl.size(); i++){
      if(eventbase->GetMuonSel()->MuonPass(muonLooseColl.at(i), "MUON_HN_TRI_TIGHT")) muonTightColl.push_back( muonLooseColl.at(i) );
    }

    // ========== Electron Energy systematics ====================
    std::vector<snu::KElectron> electronLooseColl, electronTightColl;
    if(this_syst == "ElectronEn_up"){
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedUp(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        double new_reliso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedUp();
        if( this_electron.Pt() >= 10. && new_reliso < this_reliso_el ) electronLooseColl.push_back( this_electron );
      }
    }
    else if(this_syst == "ElectronEn_down"){
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedDown(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        double new_reliso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedDown();
        if( this_electron.Pt() >= 10. && new_reliso < this_reliso_el ) electronLooseColl.push_back( this_electron );
      }
    }
    else{
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        if( this_electron.Pt() >= 10. && this_electron.PFRelIso(0.3) < this_reliso_el ) electronLooseColl.push_back( this_electron );
      }
    }
    for(unsigned int i=0; i<electronLooseColl.size(); i++){
      if(eventbase->GetElectronSel()->ElectronPass(electronLooseColl.at(i), "ELECTRON16_POG_FAKELOOSE_CC_d0")) electronTightColl.push_back( electronLooseColl.at(i) );
    }

    // ========== Jet Energy systematics ====================
    std::vector<snu::KJet> jetTightColl;
    if(this_syst == "JetEn_up"){
      for(unsigned int i=0; i<jetVLooseColl.size(); i++){
        snu::KJet this_jet = jetVLooseColl.at(i);
        double this_E = this_jet.E()*this_jet.ScaledUpEnergy();
        double this_3p = sqrt(this_E*this_E-this_jet.M()*this_jet.M());
        double this_3p_sf = this_3p/this_jet.P();
        this_jet.SetPxPyPzE( this_3p_sf*this_jet.Px(), this_3p_sf*this_jet.Py(), this_3p_sf*this_jet.Pz(), this_E);
        if(this_jet.Pt() >= 30.) jetTightColl.push_back(this_jet);
      }
    }
    else if(this_syst == "JetEn_down"){
      for(unsigned int i=0; i<jetVLooseColl.size(); i++){
        snu::KJet this_jet = jetVLooseColl.at(i);
        double this_E = this_jet.E()*this_jet.ScaledDownEnergy();
        double this_3p = sqrt(this_E*this_E-this_jet.M()*this_jet.M());
        double this_3p_sf = this_3p/this_jet.P();
        this_jet.SetPxPyPzE( this_3p_sf*this_jet.Px(), this_3p_sf*this_jet.Py(), this_3p_sf*this_jet.Pz(), this_E);
        if(this_jet.Pt() >= 30.) jetTightColl.push_back(this_jet);
      }
    }
    else{
      for(unsigned int i=0; i<jetVLooseColl.size(); i++){
        snu::KJet this_jet = jetVLooseColl.at(i);
        if(this_jet.Pt() >= 30.) jetTightColl.push_back(this_jet);
      }
    }
    n_bjets=0;
    for(int j=0; j<jetTightColl.size(); j++){
      if(jetTightColl.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) n_bjets++;
    }
    float btag_sf = 1.0;//BTagScaleFactor_1a_Weighted(jetTightColl, snu::KJet::CSVv2, snu::KJet::Medium);//FIXME

    // ========== Muon ID Scalefactor systematics ====================
    double muon_id_iso_sf = 1.0;
    if(this_syst=="MuonIDSF_up"){
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonLooseColl, 1.); 
    }
    else if(this_syst=="MuonIDSF_down"){
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonLooseColl, -1.);
    }
    else{
      muon_id_iso_sf = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonLooseColl, 0);
    }

    // ========== Muon Tracking Efficiency systematics ====================
    float muon_trk_eff;
    muon_trk_eff = mcdata_correction->MuonTrackingEffScaleFactor(muonLooseColl);

    // ========== Electron ID Scalefactor systematics ====================FIXME
    double electron_id_iso_sf;
    electron_id_iso_sf = mcdata_correction->ElectronScaleFactor("ELECTRON16_FR_POG_TIGHT_CC", electronLooseColl);

    // ========== Electron RECO Scalefactor systematics ====================FIXME
    double electron_reco_sf;
    electron_reco_sf  = mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);

    // FIXME Btagging scalefactor needed

    // ========== MET ====================
    METPt = CorrectedMETRochester(muonLooseColl, METPt, METPhi, true);
    METPhi = CorrectedMETRochester(muonLooseColl, METPt, METPhi, false);
    MET.SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), 0., METPt);
 
    // ========== Pileup systematics ====================
    if(this_syst == "PU_up"){
      this_weight *= pileup_reweight_up;
    }
    else if(this_syst == "PU_down"){
      this_weight *= pileup_reweight_down;
    }
    else{
      this_weight *= pileup_reweight;
    }

    // ========== weight ====================
    this_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON16_FR_POG_TIGHT_CC", 1);
    double weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON16_FR_POG_TIGHT_CC", 1);

    /*####################################################################################################
    ##		        Analysis Code 								        ##
    ##				For SameSign MuMuE Channel Analysis				        ##
    ####################################################################################################*/

    if( !((muonLooseColl.size() == 2) && (electronLooseColl.size() == 1)) ) return;

    RAWmu[0] = muonLooseColl.at(0);
    RAWmu[1] = muonLooseColl.at(1);
    RAWel = electronLooseColl.at(0);

    if( RAWmu[0].Charge() != RAWmu[1].Charge() ) return;
    if( RAWmu[0].Charge() == RAWel.Charge() ) return;
    if( RAWmu[1].Charge() == RAWel.Charge() ) return;

    if( RAWmu[0].Pt() < 20 || RAWmu[1].Pt() < 10 || RAWel.Pt() < 10 ) return;

    if(n_bjets> 0) return;


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

    double pt0(0.), pt1(0.), pt2(0.);

    if(electronLooseColl.at(0).Pt()>muonLooseColl.at(0).Pt()){
      pt0 = electronLooseColl.at(0).Pt();
      pt1 = muonLooseColl.at(0).Pt();
      pt2 = muonLooseColl.at(1).Pt();
    }
    else if(electronLooseColl.at(0).Pt()>muonLooseColl.at(1).Pt()){
      pt0 = muonLooseColl.at(0).Pt();
      pt1 = electronLooseColl.at(0).Pt();
      pt2 = muonLooseColl.at(1).Pt();
    }
    else{ 
      pt0 = muonLooseColl.at(0).Pt();
      pt1 = muonLooseColl.at(1).Pt();
      pt2 = electronLooseColl.at(0).Pt();
    }

    double cutop[100];
    cutop[0] = pt0;
    cutop[1] = pt1;
    cutop[2] = pt2;
    cutop[3] = RECOHN[0].M();
    cutop[4] = RECOHN[1].M();
    cutop[5] = RECOHN[2].M();
    cutop[6] = RECOHN[3].M();
    cutop[7] = RECOW_pri_lowmass.M();
    cutop[8] = RECOW_pri_highmass.M();
    cutop[9] = RECOW_sec_lowmass.M();
    cutop[10] = RECOW_sec_highmass.M();
    cutop[11] = this_weight;
    cutop[12] = 0.; //weight_err
    cutop[13] = METPt;
    cutop[14] = n_bjets;

    FillNtp("ntuple_"+this_syst,cutop);


  }
  return;

}// End of execute event loop
  


void HNSSSFMuMuE_FR_ntuple::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNSSSFMuMuE_FR_ntuple::BeginCycle() throw( LQError ){
  
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

HNSSSFMuMuE_FR_ntuple::~HNSSSFMuMuE_FR_ntuple() {
  
  Message("In HNSSSFMuMuE_FR_ntuple Destructor" , INFO);
  
}


void HNSSSFMuMuE_FR_ntuple::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNSSSFMuMuE_FR_ntuple::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNSSSFMuMuE_FR_ntupleCore::MakeHistograms() to make new hists for your analysis
   **/

  MakeNtp("ntuple_MuonEn_up",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_MuonEn_down",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_ElectronEn_up",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_ElectronEn_down",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");

  MakeNtp("ntuple_JetEn_up",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_JetEn_down",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_JetRes_up",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_JetRes_down",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_Unclustered_up",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_Unclustered_down",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_Central",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_MuonIDSF_up",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_MuonIDSF_down",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_PU_up",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
  MakeNtp("ntuple_PU_down",  "first_pt:second_pt:third_pt:HN_1_mass:HN_2_mass:HN_3_mass:HN_4_mass:W_pri_lowmass_mass:W_pri_highmass_mass:W_sec_lowmass_mass:W_sec_highmass_mass:weight:weight_err:PFMET:nbjets");
 
}


void HNSSSFMuMuE_FR_ntuple::DrawHistograms(TString suffix, double weight){

  FillHist("number_of_events_"+suffix, 0., weight, 0., 1., 1);
  FillHist("W_primary_lowmass_"+suffix, RECOW_pri_lowmass.M(), weight, 0., 1000., 2000);
  FillHist("W_secondary_lowmass_"+suffix, RECOW_sec_lowmass.M(), weight, 0., 1000., 1000);
  FillHist("W_primary_highmass_"+suffix, RECOW_pri_highmass.M(), weight, 0., 1000., 2000);
  FillHist("W_secondary_highmass_"+suffix, RECOW_sec_highmass.M(), weight, 0., 1000., 1000);
  FillHist("HN_mass_class1_"+suffix, RECOHN[0].M(), weight, 0., 500., 500);
  FillHist("HN_mass_class2_"+suffix, RECOHN[1].M(), weight, 0., 500., 500);
  FillHist("HN_mass_class3_"+suffix, RECOHN[2].M(), weight, 0., 800., 800);
  FillHist("HN_mass_class4_"+suffix, RECOHN[3].M(), weight, 0., 1500., 1500);
  FillHist("NBjets_"+suffix, n_bjets, weight, 0., 5., 5);
  FillHist("[SignalStudy]deltaR_elMET_"+suffix, RECOel.DeltaR(MET), weight, 0., 5., 100);
  FillHist("[SignalStudy]transversemass_elMET_"+suffix, MT(MET,RECOel), weight, 0., 1000., 2000);

  return;

}


void HNSSSFMuMuE_FR_ntuple::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
 
  RAWmu[0].SetPxPyPzE(0,0,0,0); RAWmu[1].SetPxPyPzE(0,0,0,0); RAWel.SetPxPyPzE(0,0,0,0); RAWnu[0].SetPxPyPzE(0,0,0,0); RAWnu[1].SetPxPyPzE(0,0,0,0);
  RECOmu[0].SetPxPyPzE(0,0,0,0); RECOmu[1].SetPxPyPzE(0,0,0,0); RECOel.SetPxPyPzE(0,0,0,0); RECOnu_lowmass.SetPxPyPzE(0,0,0,0); RECOnu_highmass.SetPxPyPzE(0,0,0,0); RECOW_pri_lowmass.SetPxPyPzE(0,0,0,0); RECOW_sec_lowmass.SetPxPyPzE(0,0,0,0); RECOW_pri_highmass.SetPxPyPzE(0,0,0,0); RECOW_sec_highmass.SetPxPyPzE(0,0,0,0); RECOHN[0].SetPxPyPzE(0,0,0,0); RECOHN[1].SetPxPyPzE(0,0,0,0); RECOHN[2].SetPxPyPzE(0,0,0,0); RECOHN[3].SetPxPyPzE(0,0,0,0);
  MET.SetPxPyPzE(0,0,0,0);
  n_bjets = 0; 

  out_muons.clear();
  out_electrons.clear();
}


void HNSSSFMuMuE_FR_ntuple::EventSelectionStudy( snu::KParticle RAWmu[], snu::KParticle RAWel, int signal_class ){

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


int HNSSSFMuMuE_FR_ntuple::DefineClass(){

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


int HNSSSFMuMuE_FR_ntuple::GetPeriodIndex(void){
  if( isData ){
    if( k_sample_name.Contains("B") ||  k_sample_name.Contains("C") || k_sample_name.Contains("D") || k_sample_name.Contains("E") || k_sample_name.Contains("F") ){
      return 1;
    }
    else return 7;
  }
  else return 0;
}

TString HNSSSFMuMuE_FR_ntuple::GetSystematicString( int it_sys ){

  TString this_syst;

  if(it_sys==0){
    this_syst = "MuonEn_up";////MET = event.PFMETShifted(snu::KEvent::MuonEn, snu::KEvent::up);
  }
  else if(it_sys==1){
    this_syst = "MuonEn_down";////MET = event.PFMETShifted(snu::KEvent::MuonEn, snu::KEvent::down);
  }
  else if(it_sys==2){
    this_syst = "ElectronEn_up";//
  }
  else if(it_sys==3){
    this_syst = "ElectronEn_down";//
  }
  else if(it_sys==4){
    this_syst = "JetEn_up";//
  }
  else if(it_sys==5){
    this_syst = "JetEn_down";//
  }
  else if(it_sys==6){
    this_syst = "JetRes_up";//
  }
  else if(it_sys==7){
    this_syst = "JetRes_down";//
  }
  else if(it_sys==8){//
    this_syst = "Unclustered_up";
  }
  else if(it_sys==9){//
    this_syst = "Unclustered_down";
  }
  else if(it_sys==10){//
    this_syst = "Central";
  }
  else if(it_sys==11){//
    this_syst = "MuonIDSF_up";
  }
  else if(it_sys==12){//
    this_syst = "MuonIDSF_down";
  }
  else if(it_sys==13){
    this_syst = "PU_up";
  }
  else if(it_sys==14){
    this_syst = "PU_down";
  }

  return this_syst;

}
