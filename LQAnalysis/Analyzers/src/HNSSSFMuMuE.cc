// $Id: HNSSSFMuMuE.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNSSSFMuMuE Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNSSSFMuMuE.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNSSSFMuMuE);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNSSSFMuMuE::HNSSSFMuMuE() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNSSSFMuMuE");
  
  Message("In HNSSSFMuMuE constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNSSSFMuMuE::InitialiseAnalysis() throw( LQError ) {
  
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


void HNSSSFMuMuE::ExecuteEvents()throw( LQError ){
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
  if(!k_isdata){ pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);} 
  TString mumu_trigger="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";
  vector<TString> trignames;
  trignames.push_back(mumu_trigger);

  float weight_trigger = WeightByTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", TargetLumi);

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_NOCUT");
  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
  int nbjet = NBJet(GetJets("JET_HN"));

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TRI_TIGHT",false);

  std::vector<snu::KElectron> electronLooseColl = GetElectrons(false,false,"ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronTightColl = GetElectrons(false,false,"ELECTRON_HN_TIGHT");

  bool trig_pass=PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
  if(!trig_pass) return;

  CorrectMuonMomentum(muonLooseColl);
  float electron_idsf = mcdata_correction->ElectronScaleFactor("ELECTRON_HN_FAKELOOSE", electronLooseColl);
  float electron_reco = mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);

  if(!isData){
    weight *= weight_trigger;
    weight *= pileup_reweight;
    weight *= electron_idsf;
    weight *= electron_reco;
  }


  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();

  MET.SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), 0., METPt);

  if( !((electronLooseColl.size() == 1) && (electronTightColl.size() == 1)) ) return;
  if( !((muonLooseColl.size() == 2) && (muonTightColl.size() == 2)) ) return;

  RAWmu[0] = muonLooseColl.at(0);
  RAWmu[1] = muonLooseColl.at(1);
  RAWel = electronLooseColl.at(0);

  if( RAWmu[0].Charge() != RAWmu[1].Charge() ) return;
  if( RAWmu[0].Charge() == RAWel.Charge() ) return;
  if( RAWmu[1].Charge() == RAWel.Charge() ) return;

  if( RAWmu[0].Pt() < 20 || RAWmu[1].Pt() < 10 || RAWel.Pt() < 10 ) return;

  if( ((RAWmu[0]+RAWmu[1]).M() < 4) || ((RAWmu[0]+RAWel).M() < 4) || ((RAWmu[1]+RAWel).M() < 4) ) return;

  if( k_sample_name.Contains( "HN_SSSF_" ) ){

    GENSignalStudy(false);

    bool electron_matched = DoMatchingBydR( GENel, RAWel );
    int muon_matched = DoMatchingBydR( GENmu, RAWmu );

    if( !(electron_matched) ) FillHist("ElectronMatching", 0., 1., 0., 2., 2);

    else{

      FillHist("ElectronMatching", 1., 1., 0., 2., 2);

      if( muon_matched == 1 ){
        FillHist("MuonMatching", 1., 1., 0., 2., 2);
      }
      else if( muon_matched == -1 ){
        snu::KParticle TEMPmu;
        TEMPmu = RAWmu[0];
        RAWmu[0] = RAWmu[1];
        RAWmu[1] = TEMPmu;

        FillHist("MuonMatching", 1., 1., 0., 2., 2);
      }
      else if( muon_matched == 0 ){
        FillHist("MuonMatching", 0., 1., 0., 2., 2);
        return;
      }

    }

    GENSignalStudy(true);

  }


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

  // ========== CLASS 3 =====================================
  EventSelectionStudy(RAWmu, RAWel, 3);// RECO particles output
  RECOHN[2] = RECOmu[1] + RECOel + RECOnu_highmass;

  // ========== CLASS 4 =====================================
  EventSelectionStudy(RAWmu, RAWel, 4);// RECO particles output
  RECOHN[3] = RECOmu[1] + RECOel + RECOnu_highmass;


  DrawHistograms("cut0", weight);
  FillCLHist(sssf_mumue, "cut0", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);

  // Low mass region cuts
  if( RECOW_pri_lowmass.M() < 150. ){
    DrawHistograms("cutW150", weight);
    FillCLHist(sssf_mumue, "cutW150", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);
  }




  // High mass region cuts
  if( METPt > 20. ){
    DrawHistograms("cutMET20", weight);
    FillCLHist(sssf_mumue, "cutMET20", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);
  }
  //virtual W mass cut can also be used



  return;

}// End of execute event loop
  


void HNSSSFMuMuE::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNSSSFMuMuE::BeginCycle() throw( LQError ){
  
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

HNSSSFMuMuE::~HNSSSFMuMuE() {
  
  Message("In HNSSSFMuMuE Destructor" , INFO);
  
}


void HNSSSFMuMuE::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNSSSFMuMuE::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNSSSFMuMuECore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNSSSFMuMuE::DrawHistograms(TString suffix, double weight){

  FillHist("number_of_events_"+suffix, 0., weight, 0., 1., 1);
  FillHist("transversemass_elMET_"+suffix, MT(MET,RECOel), weight, 0., 1000., 2000);
  FillHist("W_primary_lowmass_"+suffix, RECOW_pri_lowmass.M(), weight, 0., 1000., 2000);
  FillHist("W_secondary_lowmass_"+suffix, RECOW_sec_lowmass.M(), weight, 0., 1000., 1000);
  FillHist("W_primary_highmass_"+suffix, RECOW_pri_highmass.M(), weight, 0., 1000., 2000);
  FillHist("W_secondary_highmass_"+suffix, RECOW_sec_highmass.M(), weight, 0., 1000., 1000);
  FillHist("HN_mass_class1_"+suffix, RECOHN[0].M(), weight, 0., 500., 500);
  FillHist("HN_mass_class2_"+suffix, RECOHN[1].M(), weight, 0., 500., 500);
  FillHist("HN_mass_class3_"+suffix, RECOHN[2].M(), weight, 0., 800., 800);
  FillHist("HN_mass_class4_"+suffix, RECOHN[3].M(), weight, 0., 1500., 1500);

  return;

}


void HNSSSFMuMuE::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
 
  GENmu[0].SetPxPyPzE(0,0,0,0); GENmu[1].SetPxPyPzE(0,0,0,0); GENel.SetPxPyPzE(0,0,0,0); GENnu.SetPxPyPzE(0,0,0,0); GENHN.SetPxPyPzE(0,0,0,0);
  RAWmu[0].SetPxPyPzE(0,0,0,0); RAWmu[1].SetPxPyPzE(0,0,0,0); RAWel.SetPxPyPzE(0,0,0,0); RAWnu[0].SetPxPyPzE(0,0,0,0); RAWnu[1].SetPxPyPzE(0,0,0,0);
  RECOmu[0].SetPxPyPzE(0,0,0,0); RECOmu[1].SetPxPyPzE(0,0,0,0); RECOel.SetPxPyPzE(0,0,0,0); RECOnu_lowmass.SetPxPyPzE(0,0,0,0); RECOnu_highmass.SetPxPyPzE(0,0,0,0); RECOW_pri_lowmass.SetPxPyPzE(0,0,0,0); RECOW_sec_lowmass.SetPxPyPzE(0,0,0,0); RECOW_pri_highmass.SetPxPyPzE(0,0,0,0); RECOW_sec_highmass.SetPxPyPzE(0,0,0,0); RECOHN[0].SetPxPyPzE(0,0,0,0); RECOHN[1].SetPxPyPzE(0,0,0,0); RECOHN[2].SetPxPyPzE(0,0,0,0); RECOHN[3].SetPxPyPzE(0,0,0,0);
  MET.SetPxPyPzE(0,0,0,0);
 
  out_muons.clear();
  out_electrons.clear();
}


void HNSSSFMuMuE::GENSignalStudy( bool doGENEventSelection ){

  if( doGENEventSelection ){
    GENEventSelectionStudy(GENmu, GENel, GENnu, GENHN);
    return;
  }

  bool is_MuMuE = false, is_MuEMu = false;
  if( k_sample_name.Contains( "HN_SSSF_MuMuE" ) ) is_MuMuE = true;
  if( k_sample_name.Contains( "HN_SSSF_MuEMu" ) ) is_MuEMu = true;


  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

  int electron_ID = 11, muon_ID = 13, neutrino_el_ID = 12, neutrino_mu_ID = 14, W_ID = 24, HN_ID = 9900012;

/*  cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
  cout << "size : " << truthColl.size() << endl;
  cout << "index\tID\tm_index\tm_ID" << endl;
  for(int i=2; i<truthColl.size(); i++){
    if(truthColl.at(i).PdgId() != 21 && truthColl.at(i).PdgId() != 22)
    cout << i << '\t' << truthColl.at(i).PdgId() << '\t' << truthColl.at(i).IndexMother() << '\t' << truthColl.at(truthColl.at(i).IndexMother()).PdgId() << endl;
  }*/


  int max = truthColl.size();
  vector<int> HN_index, onshell_W_index, mu1_index, mu2_index, el_index, nu_index;
  HN_index.clear(); onshell_W_index.clear(); mu1_index.clear(); mu2_index.clear(); el_index.clear(); nu_index.clear();


  // Look for Heavy Neutrino using PdgId
  for( int i = 2 ; i < max ; i++ ){
    if( abs(truthColl.at(i).PdgId()) == HN_ID ){
      HN_index.push_back(i);
      GENFindDecayIndex( truthColl, i, HN_index );
      break;
    }
  }
  if( HN_index.size() == 0 ){
    FillHist("[GEN]HN_not_found", 0., 1., 0., 1., 1);
    return;
  }
  int HN_mother_index = truthColl.at( HN_index.at(0) ).IndexMother(); // Save HN mother index for muon1

  // Look for muon1 using PdgId and mother index (sister : HN)
  for( int i = 2 ; i < max ; i++ ){
    if( abs(truthColl.at(i).PdgId()) == muon_ID && truthColl.at(i).IndexMother() == HN_mother_index ){
      mu1_index.push_back(i);
      GENFindDecayIndex( truthColl, i, mu1_index);
      break;
    }
  }
  if( mu1_index.size () == 0 ){
    FillHist("[GEN]mu1_not_found", 0., 1., 0., 1., 1);
    return;
  }


  if(is_MuMuE){

    // Look for muon2 using PdgId and mother index (mother : HN)
    for( int index = 0 ; index < HN_index.size() ; index++ ){ // One of the HN index decays into muon
      for( int i = 2 ; i < max ; i++ ){
        if( abs(truthColl.at(i).PdgId()) == muon_ID && truthColl.at(i).IndexMother() == HN_index.at(index) ){
          mu2_index.push_back(i);
          GENFindDecayIndex( truthColl, i, mu2_index);
          break;
        }
      }
    }
    if( mu2_index.size () == 0 ){
      FillHist("[GEN]mu2_not_found", 0., 1., 0., 1., 1);
      return;
    }

    // Look for neutrino using PdgId
    for( int i = 2 ; i < max ; i++ ){
      if( abs(truthColl.at(i).PdgId()) == neutrino_el_ID ){
        nu_index.push_back(i);
        GENFindDecayIndex( truthColl, i, nu_index);
        break;
      }
    }
    int nu_mother_index = truthColl.at( nu_index.at(0) ).IndexMother();

    // Look for electron using PdgId and mother index (sister : neutrino)
    for( int i = 2 ; i < max ; i++ ){
      if( abs(truthColl.at(i).PdgId()) == electron_ID && truthColl.at(i).IndexMother() == nu_mother_index ){
        el_index.push_back(i);
        GENFindDecayIndex( truthColl, i, el_index);
        break;
      }
    }
  }


  if(is_MuEMu){
    // Look for electron using PdgId and mother index (mother : HN)
    for( int index = 0 ; index < HN_index.size() ; index++ ){ // One of the HN index decays into muon
      for( int i = 2 ; i < max ; i++ ){
        if( abs(truthColl.at(i).PdgId()) == electron_ID && truthColl.at(i).IndexMother() == HN_index.at(index) ){
          el_index.push_back(i);
          GENFindDecayIndex( truthColl, i, el_index);
          break;
        }
      }
    }
    if( el_index.size () == 0 ){
      FillHist("[GEN]el_not_found", 0., 1., 0., 1., 1);
      return;
    }

    // Look for neutrino using PdgId
    for( int i = 2 ; i < max ; i++ ){
      if( abs(truthColl.at(i).PdgId()) == neutrino_mu_ID ){
        nu_index.push_back(i);
        GENFindDecayIndex( truthColl, i, nu_index);
        break;
      }
    }
    int nu_mother_index = truthColl.at( nu_index.at(0) ).IndexMother();

    // Look for muon2 using PdgId and mother index (sister : neutrino)
    for( int i = 2 ; i < max ; i++ ){
      if( abs(truthColl.at(i).PdgId()) == muon_ID && truthColl.at(i).IndexMother() == nu_mother_index ){
        mu2_index.push_back(i);
        GENFindDecayIndex( truthColl, i, mu2_index);
        break;
      }
    }
  }



  GENmu[0] = truthColl.at( mu1_index.back() );
  GENmu[1] = truthColl.at( mu2_index.back() );
  GENel = truthColl.at( el_index.back() );
  GENnu = truthColl.at( nu_index.back() ); 
  GENHN = truthColl.at( HN_index.back() );

  return;

}


void HNSSSFMuMuE::GENEventSelectionStudy( snu::KParticle GENmu[], snu::KParticle GENel , snu::KParticle GENnu, snu::KParticle GENHN ){

  FillHist("[GEN]mu1_Pt", GENmu[0].Pt(), 1., 0., 500., 500);
  FillHist("[GEN]mu2_Pt", GENmu[1].Pt(), 1., 0., 500., 500);
  FillHist("[GEN]el_Pt", GENel.Pt(), 1., 0., 500., 500);
  FillHist("[GEN]HN_mass", GENHN.M(), 1., 0., 1500., 1500);
  FillHist("[GEN]nu_Pt", GENnu.Pt(), 1., 0., 500., 500);
  FillHist("[GEN]nu_Pz", GENnu.Pz(), 1., -500., 500., 1000);

  if( GENmu[0].Pt() > GENmu[1].Pt() ) FillHist("[GEN]Pt_mu1_bigger_than_mu2", 1, 1., 0., 2., 2);
  else  FillHist("[GEN]Pt_mu1_bigger_than_mu2", 0., 1., 0., 2., 2);

  if( GENmu[0].DeltaR(GENel) > GENmu[1].DeltaR(GENel) ) FillHist("[GEN]dR_el_mu1_bigger_than_mu2", 1, 1., 0., 2., 2);
  else FillHist("[GEN]dR_el_mu1_bigger_than_mu2", 0, 1., 0., 2., 2);

  FillHist("[GEN]HN_mass_after_decay", (GENmu[1]+GENel+GENnu).M(), 1., 0., 1500., 1500);
  
  return;

}


void HNSSSFMuMuE::EventSelectionStudy( snu::KParticle RAWmu[], snu::KParticle RAWel, int signal_class ){

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


void HNSSSFMuMuE::GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index ){

  for( int i = it+1 ; i < truthColl.size(); i ++ ){
    if( truthColl.at(i).IndexMother() == it && truthColl.at(i).PdgId() == truthColl.at(it).PdgId() ){
      index.push_back(i);
    }
  }
  return;

}


int HNSSSFMuMuE::DefineClass(){

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
