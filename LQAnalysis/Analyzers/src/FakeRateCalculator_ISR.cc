// $Id: FakeRateCalculator_ISR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQFakeRateCalculator_ISR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "FakeRateCalculator_ISR.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (FakeRateCalculator_ISR);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
FakeRateCalculator_ISR::FakeRateCalculator_ISR() :  AnalyzerCore(), out_muons(0)  {
  
  // To have the correct name in the log:       
  SetLogName("FakeRateCalculator_ISR");
  
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void FakeRateCalculator_ISR::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  

  return;

}


void FakeRateCalculator_ISR::ExecuteEvents()throw( LQError ){

  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex       

  std::vector<TString> trigger_list_mu; trigger_list_mu.clear();
  std::vector<TString> trigger_list_el; trigger_list_el.clear();
  std::vector<TString> trigger_list_mumu; trigger_list_mumu.clear();
  std::vector<TString> trigger_list_elel; trigger_list_elel.clear();

  trigger_list_mu.push_back("HLT_Mu8_TrkIsoVVL_v");
  trigger_list_mu.push_back("HLT_Mu17_TrkIsoVVL_v");
  trigger_list_mumu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
  trigger_list_mumu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
  trigger_list_mumu.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  trigger_list_mumu.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  trigger_list_el.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
  trigger_list_el.push_back("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
  trigger_list_el.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
  trigger_list_elel.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");

  if(!PassTriggerOR(trigger_list_mu) && !PassTriggerOR(trigger_list_el) && !PassTriggerOR(trigger_list_mumu) && !PassTriggerOR(trigger_list_elel)) return;
  FillHist("#CutFlowCheck", 0., 1., 0., 5., 5);

  std::vector<snu::KMuon> LooseMu; LooseMu.clear();
  std::vector<snu::KElectron> LooseEl; LooseEl.clear();
  std::vector<snu::KMuon> LooseMu_BeforeMatching; LooseMu_BeforeMatching.clear();
  std::vector<snu::KElectron> LooseEl_BeforeMatching; LooseEl_BeforeMatching.clear();

  LooseMu_BeforeMatching = GetMuons("MUON_POG_FAKETIGHT", true, 10., 2.4);
  LooseEl_BeforeMatching = GetElectrons(true, true, "ELECTRON_POG_MEDIUM_FAKELOOSE", 15., 2.5);

  bool Is_QCDMC_or_Data = ((k_sample_name.Contains("QCD") || k_sample_name.Contains("qcd")) || isData);
  bool Is_Wjets_or_QCD = (k_sample_name.Contains("QCD") || k_sample_name.Contains("qcd") || k_sample_name.Contains("WJets"));

  std::vector<bool> is_Prompt; is_Prompt.clear();

  std::vector<snu::KTruth> truthColl; truthColl.clear();
  if(!Is_QCDMC_or_Data) truthColl=eventbase->GetTruth();

  for(unsigned int i=0; i<LooseMu_BeforeMatching.size(); i++){
    snu::KMuon this_muon = LooseMu_BeforeMatching.at(i);
    int lepton_type_mu = 9999;
    if(!Is_QCDMC_or_Data) lepton_type_mu = GetLeptonType(this_muon, truthColl);
    FillHist("#LeptonType_MUON_BeforeMatching", lepton_type_mu, 1., -6.5, 5.5, 12);

    if(lepton_type_mu > 0){
      LooseMu.push_back(this_muon);
      FillHist("#LeptonType_MUON_AfterMatching", lepton_type_mu, 1., -6.5, 5.5, 12);
      is_Prompt.push_back(true);
    }
    else is_Prompt.push_back(false);
  }
  for(unsigned int i=0; i<LooseEl_BeforeMatching.size(); i++){
    snu::KElectron this_electron = LooseEl_BeforeMatching.at(i);
    int lepton_type_el = 9999;
    if(!Is_QCDMC_or_Data) lepton_type_el = GetLeptonType(this_electron, truthColl);
    FillHist("#LeptonType_ELECTRON_BeforeMatching", lepton_type_el, 1., -6.5, 5.5, 12);

    if(lepton_type_el > 0){
      LooseEl.push_back(this_electron);
      FillHist("#LeptonType_ELECTRON_AfterMatching", lepton_type_el, 1., -6.5, 5.5, 12);
      is_Prompt.push_back(true);
    }
    else is_Prompt.push_back(false);

  }

  double this_weight = 1.;
  if(!isData){
    this_weight = weight;
    this_weight *= MCweight;
//    this_weight *= mcdata_correction->MuonTrackingEffScaleFactor(LooseMu);
    this_weight *= mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    this_weight *= mcdata_correction->ElectronRecoScaleFactor(LooseEl);
    this_weight *= GetKFactor();
  }

  CorrectMuonMomentum(LooseMu);
  CorrectedMETRochester(LooseMu);

  if(LooseMu.size() > 0 && LooseEl.size() == 0){
    FillHist("#Number_of_Muons", LooseMu.size(), 1., 0., 5., 5);
    GetFakeRatesAndCorrectionFactors(LooseMu, this_weight);
  }
  else if(LooseMu.size() == 0 && LooseEl.size() > 0){
    GetFakeRatesAndCorrectionFactors(LooseEl, this_weight);
  }
  if(PassTriggerOR(trigger_list_mumu) || PassTriggerOR(trigger_list_elel)){
    if(PassTriggerOR(trigger_list_mumu)){
      if(LooseMu.size() == 2 && LooseEl.size() == 0){
        if(Is_Wjets_or_QCD) DoMCClosureTests(LooseMu, is_Prompt, this_weight);
      }
    }
    if(PassTriggerOR(trigger_list_elel)){
      if(LooseMu.size() == 0 && LooseEl.size() == 2){
        if(Is_Wjets_or_QCD) DoMCClosureTests(LooseEl, is_Prompt, this_weight);
      }
    }
  }


  return;


}// End of execute event loop

void FakeRateCalculator_ISR::GetFakeRatesAndCorrectionFactors(std::vector<snu::KMuon> muons, double this_weight){

  std::vector<snu::KJet> jets_noveto; jets_noveto.clear();
  std::vector<snu::KJet> jets; jets.clear();

  jets_noveto = GetJets("JET_NOLEPTONVETO", 30, 2.4);

  for(unsigned int i=0; i<jets_noveto.size(); i++){
    for(unsigned int j=0; j<muons.size(); j++){
      if(muons.at(j).DeltaR(jets_noveto.at(i)) > 0.4){
        jets.push_back(jets_noveto.at(i));
      }
    }
  }

  if(muons.size() == 1){

    if(jets.size() == 0) return;

    KLepton Probe = muons.at(0);
    bool Probe_passing_Tight = PassID(muons.at(0), "MUON_POG_TIGHT");

    GetFakeRates("muon", Probe, jets, Probe_passing_Tight, this_weight);

  }

  if(muons.size() == 2){
    std::vector<KLepton> leptons;

    bool both_pass_Tight = true;
    for(unsigned int i=0; i<muons.size(); i++){
      if(!PassID(muons.at(i), "MUON_POG_TIGHT")) both_pass_Tight = false;
      leptons.push_back(muons.at(i));
    }

    GetCorrectionFactors("muon", leptons, jets, this_weight, both_pass_Tight);
    if(both_pass_Tight) GetCorrectionFactors("muon", leptons, jets, this_weight, both_pass_Tight);

  }

  return;
}

void FakeRateCalculator_ISR::GetFakeRatesAndCorrectionFactors(std::vector<snu::KElectron> electrons, double this_weight){

  std::vector<snu::KJet> jets_noveto; jets_noveto.clear();
  std::vector<snu::KJet> jets; jets.clear();

  jets_noveto = GetJets("JET_NOLEPTONVETO", 30, 2.5);

  for(unsigned int i=0; i<jets_noveto.size(); i++){
    for(unsigned int j=0; j<electrons.size(); j++){
      if(electrons.at(j).DeltaR(jets_noveto.at(i)) > 0.4){
        jets.push_back(jets_noveto.at(i));
      }
    }
  }

  if(jets.size() == 0) return;
  if(jets.at(0).Pt() < 50.) return;


  if(electrons.size() == 1){

    KLepton Probe = electrons.at(0);
    bool Probe_passing_Tight = PassID(electrons.at(0), "ELECTRON_POG_MEDIUM");

    GetFakeRates("electron", Probe, jets, Probe_passing_Tight, this_weight);

  }

  if(electrons.size() == 2){

    std::vector<KLepton> leptons;

    bool both_pass_Tight = true;
    for(unsigned int i=0; i<electrons.size(); i++){
      if(!PassID(electrons.at(i), "ELECTRON_POG_MEDIUM")) both_pass_Tight = false;
      leptons.push_back(electrons.at(i));
    }

    GetCorrectionFactors("electron", leptons, jets, this_weight, both_pass_Tight);
    if(both_pass_Tight) GetCorrectionFactors("electron", leptons, jets, this_weight, both_pass_Tight);

  }

  return;
}

void FakeRateCalculator_ISR::DoMCClosureTests(std::vector<snu::KMuon> muons, std::vector<bool> is_Prompt, double this_weight){

  if(muons.size() != 2) return;

  std::vector<KLepton> leptons; leptons.clear();
  std::vector<bool> is_Tight; is_Tight.clear();

  for(unsigned int i=0; i<muons.size(); i++){
    if(PassID(muons.at(i), "MUON_POG_TIGHT")) is_Tight.push_back(true);
    else is_Tight.push_back(false);
    leptons.push_back(muons.at(i));
  }

  DoMCClosureTests("muon", leptons, is_Tight, is_Prompt, this_weight);

  return;
}

void FakeRateCalculator_ISR::DoMCClosureTests(std::vector<snu::KElectron> electrons, std::vector<bool> is_Prompt, double this_weight){

  if(electrons.size() != 2) return;

  std::vector<KLepton> leptons; leptons.clear();
  std::vector<bool> is_Tight; is_Tight.clear();

  for(unsigned int i=0; i<electrons.size(); i++){
    if(PassID(electrons.at(i), "ELECTRON_POG_MEDIUM")) is_Tight.push_back(true);
    else is_Tight.push_back(false);
    leptons.push_back(electrons.at(i));
  }

  DoMCClosureTests("electron", leptons, is_Tight, is_Prompt, this_weight);

  return;
}


void FakeRateCalculator_ISR::GetFakeRates(TString flavour, KLepton Probe, std::vector<snu::KJet> jets, bool Probe_passing_Tight, double this_weight){

  double trig_psweight = -999.;

  if(flavour == "muon"){

    if(Probe.Pt() > 20.){
      if(PassTrigger("HLT_Mu17_TrkIsoVVL_v")) trig_psweight = 217.553;
      else return;
    }
    else if(Probe.Pt() > 10.){
      if(PassTrigger("HLT_Mu8_TrkIsoVVL_v")) trig_psweight = 7.832*1.33;
      else return;
    }
    else return;

  }
  else if(flavour == "electron"){

    if(Probe.Pt() > 25.){
      if(PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) trig_psweight = 63.046;
      else return;
    }
    else if(Probe.Pt() > 20.){
      if(PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) trig_psweight = 58.896;
      else return;
    }
    else if(Probe.Pt() > 15.){
      if(PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) trig_psweight = 14.888;
      else return;
    }
    else return;

  }
  else return;

  if(!isData) this_weight *= trig_psweight;

  snu::KJet Tag;
  bool tagprobe_found = false;

  double LeptonJet_DeltaPhi = 2.5;

  for(unsigned int i=0; i<jets.size(); i++){
    if(jets.at(i).DeltaPhi(Probe) > LeptonJet_DeltaPhi){
      Tag = jets.at(i);
      tagprobe_found = true;
      break;
    }
  }
  if(!tagprobe_found) return;

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  snu::KParticle MET;
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);
  double MT = GetTransverseMass(Probe, MET);

  std::map<TString, bool> map_cut;
  map_cut["nocut"] = true;
  map_cut["ptbal0p7met80cut"] = ((Tag.Pt()/Probe.Pt()) > 0.7) && (MET.Pt() < 80.);
  map_cut["ptbal0p7mt80cut"] = ((Tag.Pt()/Probe.Pt()) > 0.7) && (MT < 80.);
  map_cut["met80mt80cut"] = (MET.Pt() < 80.) && (MT < 80.);
  map_cut["ptbal0p7met80mt80cut"] = ((Tag.Pt()/Probe.Pt()) > 0.7) && (MET.Pt() < 80.) && (MT < 80.);

  for(std::map<TString,bool>::iterator it_map_cut = map_cut.begin(); it_map_cut != map_cut.end(); it_map_cut++){
    TString this_cut_string = it_map_cut->first;
    bool this_cut_bool = it_map_cut->second;

    if(this_cut_bool){
      GetFakeRates_FillHistograms(Probe, Tag, MET, this_cut_string, "loose", flavour, this_weight);
      if(Probe_passing_Tight) GetFakeRates_FillHistograms(Probe, Tag, MET, this_cut_string, "tight", flavour, this_weight);
    }
  }


  return;

}

  
void FakeRateCalculator_ISR::GetFakeRates_FillHistograms(KLepton Probe, snu::KJet Tag, snu::KParticle MET, TString cut, TString workingpoint, TString flavour, double this_weight){

  int nvertex = eventbase->GetEvent().nVertices();
  double MT = GetTransverseMass(Probe, MET);

  Float_t ptbins[9] = { -999, 20., 25., 30., 35., 40., 50., 60., 100.};
  Float_t etabins[4] = { 0., 0.8, 1.479, 999};
  if(flavour == "muon"){ptbins[0] = 10.; etabins[3] = 2.4;}
  if(flavour == "electron"){ptbins[0] = 15.; etabins[3] = 2.5;}

  //if(k_sample_name.Contains("QCD") || k_sample_name.Contains("qcd")) flavour = "QCD_" + flavour;
  TString hist_prefix = flavour + "_" + cut + "_" + workingpoint;
  TString hist_etabin = "NULL";
  for(unsigned int i = 0; i < 4; i ++){
    if(i==3){
      hist_etabin = "all";
      FillHist("FAKERATE_" + hist_prefix + "_probe_pt_eta_" + hist_etabin, Probe.Pt(), fabs(Probe.Eta()), this_weight, ptbins, 8, etabins, 3);
    }
    else{
      if(etabins[i] <= fabs(Probe.Eta()) && fabs(Probe.Eta()) < etabins[i+1]) hist_etabin = "etabin" + TString::Itoa(i, 10);
      else continue;
    }

    FillHist("FAKERATE_" + hist_prefix + "_nevents_" + hist_etabin, 0., this_weight, 0., 1., 1);
    FillHist("FAKERATE_" + hist_prefix + "_probe_pt_" + hist_etabin, Probe.Pt(), this_weight, 0., 200., 200);
    FillHist("FAKERATE_" + hist_prefix + "_probe_eta_" + hist_etabin, Probe.Eta(), this_weight, -3., 3., 60);
    FillHist("FAKERATE_" + hist_prefix + "_probe_phi_" + hist_etabin, Probe.Phi(), this_weight, -3.5, 3.5, 70);
    FillHist("FAKERATE_" + hist_prefix + "_probe_dxy_" + hist_etabin, Probe.dXY(), this_weight, -1., 1., 200);
    FillHist("FAKERATE_" + hist_prefix + "_probe_dz_" + hist_etabin, Probe.dZ(), this_weight, -1., 1., 200);
    FillHist("FAKERATE_" + hist_prefix + "_probe_reliso_" + hist_etabin, Probe.RelIso(), this_weight, 0., 0.6, 60);
    FillHist("FAKERATE_" + hist_prefix + "_tag_pt_" + hist_etabin, Tag.Pt(), this_weight, 0., 200., 200);
    FillHist("FAKERATE_" + hist_prefix + "_tag_eta_" + hist_etabin, Tag.Eta(), this_weight, -3., 3., 60);
    FillHist("FAKERATE_" + hist_prefix + "_tag_phi_" + hist_etabin, Tag.Phi(), this_weight, -3.5, 3.5, 70);
    FillHist("FAKERATE_" + hist_prefix + "_met_" + hist_etabin, MET.Pt(), this_weight, 0., 200., 200);
    FillHist("FAKERATE_" + hist_prefix + "_tagprobe_mass_" + hist_etabin, (Tag+Probe).M(), this_weight, 0., 200., 200);
    FillHist("FAKERATE_" + hist_prefix + "_tagprobe_pt_" + hist_etabin, (Tag+Probe).Pt(), this_weight, 0., 200., 200);
    FillHist("FAKERATE_" + hist_prefix + "_tagprobe_ptbalance_" + hist_etabin, (Tag.Pt()/Probe.Pt()), this_weight, 0., 5., 50);
    FillHist("FAKERATE_" + hist_prefix + "_mt_" + hist_etabin, MT, this_weight, 0., 250., 250);  
    FillHist("FAKERATE_" + hist_prefix + "_nvertex_" + hist_etabin, nvertex, this_weight, 0., 100., 100);
  }

  return;
}

void FakeRateCalculator_ISR::GetCorrectionFactors(TString flavour, std::vector<KLepton> leptons, std::vector<snu::KJet> jets, double this_weight, bool both_pass_Tight){

  if(leptons.at(0).Charge() == leptons.at(1).Charge()) return;
  if((leptons.at(0)+leptons.at(1)).M() < 15.) return;

  std::vector<TString> trigger_list_mu; trigger_list_mu.clear();
  std::vector<TString> trigger_list_el; trigger_list_el.clear();
  std::vector<TString> trigger_list_mumu; trigger_list_mumu.clear();
  std::vector<TString> trigger_list_elel; trigger_list_elel.clear();

  std::vector<TString> trigger_list; trigger_list.clear();

  trigger_list_mu.push_back("HLT_Mu8_TrkIsoVVL_v");
  trigger_list_mu.push_back("HLT_Mu17_TrkIsoVVL_v");

  trigger_list_el.push_back("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
  trigger_list_el.push_back("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v");
  trigger_list_el.push_back("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v");

  std::vector<double> trigger_prescale_mu; trigger_prescale_mu.clear();
  std::vector<double> trigger_prescale_el; trigger_prescale_el.clear();

  trigger_prescale_mu.push_back(7.832*1.33);
  trigger_prescale_mu.push_back(217.553);

  trigger_prescale_el.push_back(14.888);
  trigger_prescale_el.push_back(58.896);
  trigger_prescale_el.push_back(63.046);

  std::vector<double> trigger_prescale; trigger_prescale.clear();

  if(flavour == "muon"){
    if(!((leptons.at(0).Pt() > 20.) && (leptons.at(1).Pt() > 10.))) return;
    for(unsigned int i=0; i<trigger_list_mu.size(); i++){
      if(PassTrigger(trigger_list_mu.at(i))){
        trigger_list.push_back(trigger_list_mu.at(i));
	if(!isData) trigger_prescale.push_back(trigger_prescale_mu.at(i));
        else trigger_prescale.push_back(1.);
      }
    }
  }
  if(flavour == "electron"){
    if(!((leptons.at(0).Pt() > 25.) && (leptons.at(1).Pt() > 15.))) return;
    for(unsigned int i=0; i<trigger_list_el.size(); i++){
      if(PassTrigger(trigger_list_el.at(i))){
        trigger_list.push_back(trigger_list_el.at(i));
        if(!isData)trigger_prescale.push_back(trigger_prescale_el.at(i));
        else trigger_prescale.push_back(1.);
      }
    }
  }

  for(unsigned int i=0; i<trigger_list.size(); i++){

    if(SelectZpeak(leptons, 91.1876, 15)){
      GetCorrectionFactors_FillHistograms(flavour, leptons, jets, "loose", trigger_list.at(i), (this_weight*trigger_prescale.at(i)));
      if(both_pass_Tight) GetCorrectionFactors_FillHistograms(flavour, leptons, jets, "tight", trigger_list.at(i), (this_weight*trigger_prescale.at(i)));  
    }

  }

  return;
}

void FakeRateCalculator_ISR::GetCorrectionFactors_FillHistograms(TString flavour, std::vector<KLepton> leptons, std::vector<snu::KJet> jets, TString workingpoint, TString trigger, double this_weight){

  int nvertex = eventbase->GetEvent().nVertices();
  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  snu::KParticle MET;
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);

  TString hist_prefix = flavour + "_" + trigger + "_" + workingpoint;

  FillHist("CORRFACTORS_" + hist_prefix + "_ll_mass", (leptons.at(0)+leptons.at(1)).M(), this_weight, 0., 200., 200);
  FillHist("CORRFACTORS_" + hist_prefix + "_ll_mass_zoomedin", (leptons.at(0)+leptons.at(1)).M(), this_weight, 91.1876-20., 91.1876+20, 40);

  FillHist("CORRFACTORS_" + hist_prefix + "_ll_pt", (leptons.at(0)+leptons.at(1)).Pt(), this_weight, 0., 100., 100);
  FillHist("CORRFACTORS_" + hist_prefix + "_leadinglepton_pt", leptons.at(0).Pt(), this_weight, 0., 200., 200);
  FillHist("CORRFACTORS_" + hist_prefix + "_leadinglepton_eta", leptons.at(0).Eta(), this_weight, -3., 3., 60);
  FillHist("CORRFACTORS_" + hist_prefix + "_subleadinglepton_pt", leptons.at(1).Pt(), this_weight, 0., 200., 200);
  FillHist("CORRFACTORS_" + hist_prefix + "_subleadinglepton_eta", leptons.at(1).Eta(), this_weight, -3., 3., 60);

  FillHist("CORRFACTORS_" + hist_prefix + "_njets", jets.size(), this_weight, 0., 10., 10);
  if(jets.size()!=0) FillHist("CORRFACTORS_" + hist_prefix + "_leadingjet_pt", jets.at(0).Pt(), this_weight, 0., 200., 200);
  FillHist("CORRFACTORS_" + hist_prefix + "_met", METPt, this_weight, 0., 200., 200);
  FillHist("CORRFACTORS_" + hist_prefix + "_nvertex", nvertex, this_weight, 0., 100., 100);
  FillHist("CORRFACTORS_" + hist_prefix + "_nevents", 0., this_weight, 0., 1., 1);

  return;
}

void FakeRateCalculator_ISR::DoMCClosureTests(TString flavour, std::vector<KLepton> leptons, std::vector<bool> is_Tight, std::vector<bool> is_Prompt, double this_weight){

  this_weight = 1.;

  if(leptons.size() != 2) return;
  if(flavour == "muon"){
    if(leptons.at(0).Pt() < 20. || leptons.at(1).Pt() < 10.) return;
  }
  if(flavour == "electron"){
    if(leptons.at(0).Pt() < 25. || leptons.at(1).Pt() < 15.) return;
  }

  bool is_OS = false;
  if(leptons.at(0).Charge() != leptons.at(1).Charge()) is_OS = true;

  if(is_Tight.size() != is_Prompt.size() || is_Tight.size() != leptons.size()) FillHist("#DoMCClosureTests_ERROR", 0., 1., 0., 1., 1);  
/*
  if(!(is_Tight.at(0) && is_Tight.at(1))){
    TString hist_prefix = flavour + "_expected";
    double temp_weight = this_weight * GetFakeRateWeightings(flavour, leptons, is_Tight);

    FillHist("MCCLOSURE_" + hist_prefix + "_leadinglepton_pt", leptons.at(0).Pt(), temp_weight, 0., 1., 1);
    FillHist("MCCLOSURE_" + hist_prefix + "_leadinglepton_eta", leptons.at(0).Eta(), temp_weight, -3., 3., 60);
    FillHist("MCCLOSURE_" + hist_prefix + "_subleadinglepton_pt", leptons.at(1).Pt(), temp_weight, 0., 1., 1);
    FillHist("MCCLOSURE_" + hist_prefix + "_subleadinglepton_eta", leptons.at(1).Eta(), temp_weight, -3., 3., 60);
    FillHist("MCCLOSURE_" + hist_prefix + "_nevents", 0., temp_weight, 0., 1., 1);

  }
  if(!(is_Prompt.at(0) && is_Prompt.at(1))){
    TString hist_prefix = flavour + "_observed";

    FillHist("MCCLOSURE_" + hist_prefix + "_leadinglepton_pt", leptons.at(0).Pt(), this_weight, 0., 1., 1);
    FillHist("MCCLOSURE_" + hist_prefix + "_leadinglepton_eta", leptons.at(0).Eta(), this_weight, -3., 3., 60);
    FillHist("MCCLOSURE_" + hist_prefix + "_subleadinglepton_pt", leptons.at(1).Pt(), this_weight, 0., 1., 1);
    FillHist("MCCLOSURE_" + hist_prefix + "_subleadinglepton_eta", leptons.at(1).Eta(), this_weight, -3., 3., 60);
    FillHist("MCCLOSURE_" + hist_prefix + "_nevents", 0., this_weight, 0., 1., 1);

  }
*/
  return;
}

bool FakeRateCalculator_ISR::SelectZpeak(std::vector<KLepton> leptons, double zmass, double window){

  if(leptons.size() != 2) return false;
  if(fabs((leptons.at(0) + leptons.at(1)).M() - zmass) > window) return false;

  return true;
}

double FakeRateCalculator_ISR::GetTransverseMass(KLepton Lepton, snu::KParticle MET){

  double dphi = Lepton.DeltaPhi(MET);
  double MT = TMath::Sqrt( 2. * Lepton.Pt() * MET.Pt() * ( 1. - TMath::Cos(dphi) ) );
  return MT;

}

double FakeRateCalculator_ISR::GetWeightFromFakeRate(TString flavour, std::vector<KLepton> leptons, std::vector<bool> is_Tight){

  if(leptons.size() != 2.) return 0.;
  if(is_Tight.at(0) && is_Tight.at(1)) return 0.;

  double fake_weight = -1.;

  for(unsigned int i_lep = 0; i_lep<2; i_lep++){
    if(!is_Tight.at(i_lep)){
      double fake_rate = 0.;
      fake_rate = GetFakeRateForWeight(flavour, leptons.at(i_lep));
      fake_weight *= -(fake_rate/(1.-fake_rate));
    }
  }

  return fake_weight;
}

double FakeRateCalculator_ISR::GetFakeRateForWeight(TString flavour, KLepton Lepton){

  double lep_pt = Lepton.Pt(), lep_eta = fabs(Lepton.Eta());
  if(Lepton.Pt() > 100.) lep_pt= 99.;

  if(10. < lep_pt && lep_pt < 20.){
	if(0. < lep_eta && lep_eta < 0.8)	return 0.175053;
	else if(0.8 < lep_eta && lep_eta < 1.479)	return 0.212175;
	else if(1.479 < lep_eta && lep_eta < 2.4)	return 0.27744;
  }
  else if(20. < lep_pt && lep_pt < 25.){
	if(0. < lep_eta && lep_eta < 0.8)	return 0.135246;
	else if(0.8 < lep_eta && lep_eta < 1.479)	return 0.174551;
	else if(1.479 < lep_eta && lep_eta < 2.4)	return 0.249973;
  }
  else if(25. < lep_pt && lep_pt < 30.){
	if(0. < lep_eta && lep_eta < 0.8)	return 0.127072;
	else if(0.8 < lep_eta && lep_eta < 1.479)	return 0.16721;
	else if(1.479 < lep_eta && lep_eta < 2.4)	return 0.244609;
  }
  else if(30. < lep_pt && lep_pt < 35.){
	if(0. < lep_eta && lep_eta < 0.8)	return 0.120696;
	else if(0.8 < lep_eta && lep_eta < 1.479)	return 0.157759;
	else if(1.479 < lep_eta && lep_eta < 2.4)	return 0.244579;
  }
  else if(35. < lep_pt && lep_pt < 40.){
	if(0. < lep_eta && lep_eta < 0.8)	return 0.119853;
	else if(0.8 < lep_eta && lep_eta < 1.479)	return 0.156102;
	else if(1.479 < lep_eta && lep_eta < 2.4)	return 0.23302;
  }
  else if(40. < lep_pt && lep_pt < 50.){
	if(0. < lep_eta && lep_eta < 0.8)	return 0.124877;
	else if(0.8 < lep_eta && lep_eta < 1.479)	return 0.15659;
	else if(1.479 < lep_eta && lep_eta < 2.4)	return 0.242115;
  }
  else if(50. < lep_pt && lep_pt < 60.){
	if(0. < lep_eta && lep_eta < 0.8)	return 0.101897;
	else if(0.8 < lep_eta && lep_eta < 1.479)	return 0.130986;
	else if(1.479 < lep_eta && lep_eta < 2.4)	return 0.251526;
  }
  else if(60. < lep_pt && lep_pt < 100.){
	if(0. < lep_eta && lep_eta < 0.8)	return 0.117231;
	else if(0.8 < lep_eta && lep_eta < 1.479)	return 0.103267;
	else if(1.479 < lep_eta && lep_eta < 2.4)	return 0.241473;
  }
  else return -99999999.;
}

void FakeRateCalculator_ISR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void FakeRateCalculator_ISR::BeginCycle() throw( LQError ){
  
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

FakeRateCalculator_ISR::~FakeRateCalculator_ISR() {
  
  Message("In FakeRateCalculator_ISR Destructor" , INFO);
  
}


void FakeRateCalculator_ISR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void FakeRateCalculator_ISR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this FakeRateCalculator_ISRCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void FakeRateCalculator_ISR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



