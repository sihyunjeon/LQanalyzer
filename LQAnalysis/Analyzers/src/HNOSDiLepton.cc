// $Id: HNOSDiLepton.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNOSDiLepton Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNOSDiLepton.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNOSDiLepton);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNOSDiLepton::HNOSDiLepton() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNOSDiLepton");
  
  Message("In HNOSDiLepton constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNOSDiLepton::InitialiseAnalysis() throw( LQError ) {
  
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


void HNOSDiLepton::ExecuteEvents()throw( LQError ){

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
  std::vector<TString> triggerlist_mm, triggerlist_me, triggerlist_em1, triggerlist_em2, triggerlist_ee;
  triggerlist_mm.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist_mm.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  triggerlist_em1.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_em2.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
//  triggerlist_em2.push_back("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");

//  triggerlist_me1.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");
//  triggerlist_me2.push_back("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  triggerlist_ee.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  if(!(PassTriggerOR(triggerlist_mm)||PassTriggerOR(triggerlist_em1)||PassTriggerOR(triggerlist_em2)||PassTriggerOR(triggerlist_ee))) return;
  bool Pass_Trigger_mm = PassTriggerOR(triggerlist_mm);
  bool Pass_Trigger_em = PassTriggerOR(triggerlist_em1) && PassTriggerOR(triggerlist_em2);
  bool Pass_Trigger_ee = PassTriggerOR(triggerlist_ee);
  // ================================================================================


  std::vector<bool> isT;
  // ========== Get Objects (muon, electron, jet) ====================
  std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_HN_VETO", false);
  std::vector<snu::KMuon> muons;
  std::vector<snu::KMuon> muonTightColl;
  muonTightColl.clear(); muons.clear();
  int muonVetoN = muonVetoColl.size(), muonTightN = 0, muonsN = 0;
  for(int i=0;i<muonVetoN;i++){
    if(PassID(muonVetoColl.at(i), "MUON_HN_LOOSE")){
      muons.push_back(muonVetoColl.at(i));
      muonsN ++;
      if(PassID(muonVetoColl.at(i), "MUON_HN_TIGHT")){
        muonTightColl.push_back(muonVetoColl.at(i));
        muonTightN++;
        isT.push_back(true);
      }
      else isT.push_back(false);
    }
  }
  if(!k_running_nonprompt) if(muonsN != muonTightN || muonsN != muonVetoN) return;
  if(k_running_nonprompt) if(muonsN == muonTightN || muonsN != muonVetoN) return;

  std::vector<snu::KElectron> electronVetoColl = GetElectrons(false,false,"ELECTRON_HN_VETO");
  std::vector<snu::KElectron> electrons;
  std::vector<snu::KElectron> electronTightColl;
  electronTightColl.clear(); electrons.clear();
  int electronVetoN = electronVetoColl.size(), electronTightN = 0, electronsN = 0;
  for(int i=0;i<electronVetoN;i++){
    if(PassID(electronVetoColl.at(i), "ELECTRON_HN_FAKELOOSE")){
      electrons.push_back(electronVetoColl.at(i));
      electronsN ++;
      if(PassID(electronVetoColl.at(i), "ELECTRON_HN_TIGHTv4")){
        electronTightColl.push_back(electronVetoColl.at(i));
        electronTightN++;
        isT.push_back(true);
      }
      else isT.push_back(false);
    }
  }
  if(!k_running_nonprompt) if(electronsN != electronTightN || electronsN != electronVetoN) return;
  if(k_running_nonprompt) if(electronsN == electronTightN || electronsN != electronVetoN) return;

  std::vector<snu::KJet> jets = GetJets("JET_HN");

  double H_T = 0.;
  for(int i=0; i<jets.size(); i++){
    H_T += jets.at(i).Pt();
  }

  int bjetsN=0;
  for(int j=0; j<jets.size(); j++){
    if(jets.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) bjetsN++;
  }
  // ================================================================================


  // ========== Rochester Correction ====================
  CorrectMuonMomentum(muons);

  double MET = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();

  MET = CorrectedMETRochester(muons, true);
  // ================================================================================

  // ========== Pileup reweight ====================
  float pileup_reweight=(1.0);
  if(!k_isdata){ pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);}
  // ================================================================================


  // ========== Trigger reweight ====================
  float weight_trigger = WeightByTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", TargetLumi);
  // ================================================================================

  // ========== Muon tracking efficiency ====================  
  float muon_trkeff = mcdata_correction->MuonTrackingEffScaleFactor(muons);
  // ================================================================================


  // ========== Electron ID scalefactor ====================
  float electron_idsf = mcdata_correction->ElectronScaleFactor("ELECTRON_HN_TIGHTv4", electrons);
  // ================================================================================


  // ========== Electron RECO scalefactor ====================
  float electron_recosf = mcdata_correction->ElectronRecoScaleFactor(electrons);
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

  double fake_weight_mm, fake_weight_em, fake_weight_ee;

  TString LeptonConfig_mm = "NULL", ChargeConfig_mm = "NULL";
  bool Pass_Pt_mm = false;
  if(Pass_Trigger_mm){
    if(muonsN == 2 && electronsN == 0) LeptonConfig_mm = "2mu0el";
    if(muonsN == 3 && electronsN == 0) LeptonConfig_mm = "3mu0el";
    if(muonsN == 2 && electronsN == 1) LeptonConfig_mm = "2mu1el";

    if(muonsN > 1) Pass_Pt_mm = (muons.at(0).Pt() > 20 && muons.at(1).Pt() >10);

    if(muonsN == 2){
      if(muons.at(0).Charge() == muons.at(1).Charge()) ChargeConfig_mm = "SSSF";
      else ChargeConfig_mm = "OSSF";
    }
    else ChargeConfig_mm = "OSSF";
  }

  TString LeptonConfig_ee = "NULL", ChargeConfig_ee = "NULL";
  bool Pass_Pt_ee = false;
  if(Pass_Trigger_ee){
    if(muonsN == 0 && electronsN == 2) LeptonConfig_ee = "0mu2el";
    if(muonsN == 0 && electronsN == 3) LeptonConfig_ee = "0mu3el";
    if(muonsN == 1 && electronsN == 2) LeptonConfig_ee = "1mu2el";

    if(electronsN > 1) Pass_Pt_ee = (electrons.at(0).Pt() > 25 && electrons.at(1).Pt() >15);

    if(electronsN == 2){
      if(electrons.at(0).Charge() == electrons.at(1).Charge()) ChargeConfig_ee = "SSSF";
      else ChargeConfig_ee = "OSSF";
    }
    else ChargeConfig_ee = "OSSF";
  }

  TString LeptonConfig_em = "NULL", ChargeConfig_em = "NULL";
  bool Pass_Pt_em;
  if(Pass_Trigger_em){
    if(muonsN == 1 && electronsN == 1) LeptonConfig_em = "1mu1el";

    if(muonsN > 0 && electronsN > 0) Pass_Pt_em = (muons.at(0).Pt() > 10 && electrons.at(0).Pt() > 25);

    if(muonsN == 1 && electronsN == 1){
      if(muons.at(0).Charge() != electrons.at(0).Charge()) ChargeConfig_em = "OSOF";
      else ChargeConfig_em = "SSOF";
    }
  }
  if(LeptonConfig_mm == "NULL" && LeptonConfig_ee == "NULL" && LeptonConfig_em == "NULL") return;

  std::vector<KLepton> lep;
  for(int i=0; i<muonsN; i++){
    lep.push_back(muons.at(i));
  }
  for(int i=0; i<electronsN; i++){
    lep.push_back(electrons.at(i));
  }

  TString Define_Region = "NULL";
  if(Pass_Trigger_mm && LeptonConfig_mm == "2mu0el" && Pass_Pt_mm && ChargeConfig_mm == "SSSF") Define_Region = "DiMu_SS";
  if(Pass_Trigger_mm && LeptonConfig_mm == "2mu0el" && Pass_Pt_mm && ChargeConfig_mm == "OSSF") Define_Region = "DiMu_OS";
  if(Pass_Trigger_ee && LeptonConfig_ee == "0mu2el" && Pass_Pt_ee && ChargeConfig_ee == "SSSF") Define_Region = "DiEl_SS";
  if(Pass_Trigger_ee && LeptonConfig_ee == "0mu2el" && Pass_Pt_ee && ChargeConfig_ee == "OSSF") Define_Region = "DiEl_OS";
  int const DiMu_OS_cutN = 2;

  if(Define_Region == "DiMu_SS"){
    ;
  }
  if(Define_Region == "DiMu_OS"){
    for(int DiMu_OS_cut=0; DiMu_OS_cut<DiMu_OS_cutN; DiMu_OS_cut++){
      if(GetCuts(Define_Region, GetCuts_suffix(Define_Region, DiMu_OS_cut), lep, jets, bjetsN, MET)){ 
        FillHist(Define_Region+"_LeadingLepton_Pt_"+GetCuts_suffix(Define_Region, DiMu_OS_cut), lep.at(0).Pt(), weight, 0., 400., 400);
        FillHist(Define_Region+"_LeadingLepton_Eta_"+GetCuts_suffix(Define_Region, DiMu_OS_cut), lep.at(0).Eta(), weight, -3., 3., 60);
        FillHist(Define_Region+"_LeadingLepton_dXY_"+GetCuts_suffix(Define_Region, DiMu_OS_cut), lep.at(0).dXY(), weight, 0., 0.01, 100);
        FillHist(Define_Region+"_SubleadingLepton_Pt_"+GetCuts_suffix(Define_Region, DiMu_OS_cut), lep.at(1).Pt(), weight, 0., 400., 400);
        FillHist(Define_Region+"_SubleadingLepton_Eta_"+GetCuts_suffix(Define_Region, DiMu_OS_cut), lep.at(1).Eta(), weight, -3., 3., 60);
        FillHist(Define_Region+"_SubleadingLepton_dXY_"+GetCuts_suffix(Define_Region, DiMu_OS_cut), lep.at(1).dXY(), weight, 0., 0.01, 100);
        FillHist(Define_Region+"_DiLepton_mass_"+GetCuts_suffix(Define_Region, DiMu_OS_cut), (lep.at(0)+lep.at(1)).M(), weight, 0., 400., 400);
        FillHist(Define_Region+"_N_of_events_"+GetCuts_suffix(Define_Region, DiMu_OS_cut), 0., weight, 0., 1., 1);
      }
    }
  }


  double fake_weight;
  if(k_running_nonprompt){
//    fake_weight = get_eventweight(muons, electrons, isT);
  }
   
   return;
}// End of execute event loop
  
bool HNOSDiLepton::GetCuts(TString region, TString cut, std::vector<KLepton> lep, std::vector<snu::KJet> jets, int bjetsN, double MET){

  if(region == "DiMu_OS"){
    if((lep.at(0).Charge() == lep.at(1).Charge())) return false;

    if(cut == "SR_Preselection"){
      if((lep.at(0)+lep.at(1)).M() < 15) return false;
      if(jets.size() < 2) return false;
      if(jets.at(0).Pt() < 30) return false;

      return true;
    }

    if(cut == "SR_HighMass"){
      if((lep.at(0)+lep.at(1)).M() < 15) return false;
      if(((lep.at(0)+lep.at(1)).M() > 70) && ((lep.at(0)+lep.at(1)).M() < 110)) return false;
      if(MET < 35) return false;
      if(jets.size() < 2) return false;
      if(bjetsN != 0) return false;
      if(jets.at(0).Pt() < 30) return false;
      snu::KParticle W_candidate, W_selection;
      W_selection.SetPxPyPzE(0,0,0,0);
      for(int i=0; i<jets.size(); i++){
        for(int j=i+1; j<jets.size(); j++){
          W_candidate = jets.at(i)+jets.at(j);
          if(fabs(W_candidate.M() - 80.4) < fabs(W_selection.M() - 80.4)) W_selection = W_candidate;
        }
      }
      if((W_selection.M() < 50) || (W_selection.M() > 110)) return false;

      return true;
    }
  }  
}

TString HNOSDiLepton::GetCuts_suffix(TString region, int cut){

  if(region == "DiMu_OS") if(cut == 0) return "SR_HighMass";
  if(region == "DiMu_OS") if(cut == 1) return "SR_Preselection";

}

void HNOSDiLepton::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNOSDiLepton::BeginCycle() throw( LQError ){
  
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

HNOSDiLepton::~HNOSDiLepton() {
  
  Message("In HNOSDiLepton Destructor" , INFO);
  
}


void HNOSDiLepton::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNOSDiLepton::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNOSDiLeptonCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNOSDiLepton::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



