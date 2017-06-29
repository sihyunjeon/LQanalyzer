// $Id: CFRateCalculator_syst.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQCFRateCalculator_syst Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "CFRateCalculator_syst.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (CFRateCalculator_syst);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
CFRateCalculator_syst::CFRateCalculator_syst() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("CFRateCalculator_syst");
  
  Message("In CFRateCalculator_syst constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  TFile* file_madgraph = new TFile("/home/shjeon/CATanalyzer_v807/data/Fake/80X/ChargeFlip_madgraph_v807.root");
  TFile* file_powheg = new TFile("/home/shjeon/CATanalyzer_v807/data/Fake/80X/ChargeFlip_powheg_v807.root");

  HNTIGHT_CF_hist_madgraph  = (TH2F*)file_madgraph ->Get("Pt_eta_global_CF_HNTIGHT_PU")->Clone();
  HNTIGHT_CF_sampleA_hist_madgraph  = (TH2F*)file_madgraph ->Get("Pt_eta_global_CF_HNTIGHT_PU_sampleA")->Clone();
  HNTIGHT_CF_sampleB_hist_madgraph  = (TH2F*)file_madgraph ->Get("Pt_eta_global_CF_HNTIGHT_PU_sampleB")->Clone();
  MVATIGHT_CF_hist_madgraph  = (TH2F*)file_madgraph ->Get("Pt_eta_global_CF_MVATIGHT_PU")->Clone();
  MVATIGHT_CF_sampleA_hist_madgraph  = (TH2F*)file_madgraph ->Get("Pt_eta_global_CF_MVATIGHT_PU_sampleA")->Clone();
  MVATIGHT_CF_sampleB_hist_madgraph  = (TH2F*)file_madgraph ->Get("Pt_eta_global_CF_MVATIGHT_PU_sampleB")->Clone();
  
  HNTIGHT_CF_hist_powheg  = (TH2F*)file_powheg ->Get("Pt_eta_global_CF_HNTIGHT_PU")->Clone();
  HNTIGHT_CF_sampleA_hist_powheg  = (TH2F*)file_powheg ->Get("Pt_eta_global_CF_HNTIGHT_PU_sampleA")->Clone();
  HNTIGHT_CF_sampleB_hist_powheg  = (TH2F*)file_powheg ->Get("Pt_eta_global_CF_HNTIGHT_PU_sampleB")->Clone();
  MVATIGHT_CF_hist_powheg  = (TH2F*)file_powheg ->Get("Pt_eta_global_CF_MVATIGHT_PU")->Clone();
  MVATIGHT_CF_sampleA_hist_powheg  = (TH2F*)file_powheg ->Get("Pt_eta_global_CF_MVATIGHT_PU_sampleA")->Clone();
  MVATIGHT_CF_sampleB_hist_powheg  = (TH2F*)file_powheg ->Get("Pt_eta_global_CF_MVATIGHT_PU_sampleB")->Clone();

}


void CFRateCalculator_syst::InitialiseAnalysis() throw( LQError ) {
  
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


void CFRateCalculator_syst::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  if(!PassMETFilter()) return;     /// Initial event cuts :
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex

  // ========== Pileup reweight ====================
  float pileup_reweight=(1.0);
  if(!k_isdata){
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    weight *= pileup_reweight;
  }
  // ================================================================================



  CFvalidation();
  if(isData)  return;
 
  bool trig_pass = (PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v") || PassTrigger("HLT_Ele17_CaloIdL_GsfTrkIdVL_v"));
  TString s_trigger = "single";
  if(!trig_pass) return;

  TString sample_suffix = "";
  if(k_sample_name.Contains("DYtoEE")) sample_suffix = "_powheg";
  else if(k_sample_name.Contains("DYJets_MG")) sample_suffix = "_madgraph";
  else return;

  TString s_sample = "none", s_h_sample = "", s_hh_sample = "";
  int event_number = eventbase->GetEvent().EventNumber();
  int ran_num = event_number%2;
  if(ran_num == 0){
    s_sample = "_sampleA";
    s_h_sample = "sampleB";
    s_hh_sample = "sampleA";
  }
  else if(ran_num==1){
    s_sample = "_sampleB";
    s_h_sample = "sampleA";
    s_hh_sample = "sampleB";
  }

  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN", 30., 2.4);
  int Njets = jetTightColl.size();
  double HT = 0.;
  for(int i=0; i<Njets; i++){
    HT += jetTightColl.at(i).Pt();
  }

  double MET = eventbase->GetEvent().MET();

  for(int aa=0; aa<1; aa++){
    if(aa==0) weight = pileup_reweight;
    if(aa==1) weight = 1.;

    for(int aaa=0; aaa<1; aaa++){

      std::vector<snu::KElectron> electronTightColl, electronPromptColl;
      electronPromptColl.clear();

      TString IDsuffix = "";
      TString el_ID = "";

      if(aaa == 0){
        el_ID = "ELECTRON_HN_TIGHT";
        IDsuffix = "_HNTIGHT";
      }
      if(aaa == 1){
        el_ID = "ELECTRON_MVA_TIGHT";
        IDsuffix = "_MVATIGHT";
      }

      electronTightColl = GetElectrons(true, false, el_ID);

      if(aa == 0){IDsuffix += "_PU";}
      if(aa == 1){IDsuffix += "_PU0";}

      IDsuffix += sample_suffix;

      double LT = 0.;
      for(int i=0; i<electronTightColl.size(); i++){
        snu::KElectron this_lep;
        this_lep = electronTightColl.at(i);
        if((this_lep.Pt() < 25)) continue;

        FillHist("[CHECK]lepton_eta", this_lep.Eta(), weight, -5., 5., 1000);
	FillHist("[CHECK]lepton_SCeta", this_lep.SCEta(), weight, -5., 5., 1000);

        if((fabs(this_lep.Eta()) < 1.556) && (fabs(this_lep.Eta()) > 1.4442)) continue;
        if((this_lep.MCIsPrompt())){
          electronPromptColl.push_back(this_lep);
          LT += this_lep.Pt();
        }
      }

      if(electronPromptColl.size() == 0) continue;

      int is_region = 0;
      bool is_CF = false;
      bool is_CONV0 = false;

      float etaarray [] = {0.0, 0.9, 1.4442, 1.556, 2.5};
      float ptarray [] = {20., 40., 60., 80., 100., 200., 500.};

      for(int i=0; i<electronPromptColl.size(); i++){
        is_region = 0;
        is_CF = false;
        is_CONV0 = false;
        snu::KElectron this_lep;
        this_lep = electronPromptColl.at(i);

        //return objects : eta, pt, nonprompt
        if( (fabs(this_lep.Eta()) < 0.9) )                                        is_region = 1;
        else if( (fabs(this_lep.Eta()) < 1.4442) )                                is_region = 2;
        else if( (fabs(this_lep.Eta()) > 1.556) && (fabs(this_lep.Eta()) < 2.5) ) is_region = 3;
	else continue;

        TString s_region = "";
        if( is_region == 1 ) s_region = "_region1";
        else if( is_region == 2 ) s_region = "_region2";
        else if( is_region == 3 ) s_region = "_region3";
        else FillHist("[WARNING]suffix_not_defined", 0., 1., 0., 1., 1);

        if( (MCIsCF(this_lep)) ) is_CF = true;
        if( !(this_lep.MCIsFromConversion()) ) is_CONV0 = true;

        FillHist("Pt_eta_global"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
        if(Njets == 0) FillHist("Pt_eta_global_JETS0"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
        if(Njets != 0) FillHist("Pt_eta_global_JETS"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
        FillHist("n_events_global"+IDsuffix, 0., weight, 0., 1., 1);
        FillHist("Pt_global"+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
        FillHist("invPt_global"+IDsuffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist("eta_global"+IDsuffix, this_lep.Eta(), weight, -3., 3., 120);
        FillHist("dXY_global"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
        FillHist("HT_global"+IDsuffix, HT, weight, 0., 1000., 1000);
        FillHist("MET_global"+IDsuffix, MET, weight, 0., 1000., 1000);
        FillHist("LT_global"+IDsuffix, LT, weight, 0., 1000., 1000);
        FillHist("energy_global"+IDsuffix, this_lep.E(), weight, 0., 500., 500);

        FillHist("n_events"+s_region+IDsuffix, 0., weight, 0., 1., 1);
        FillHist("Pt"+s_region+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
        FillHist("invPt"+s_region+IDsuffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist("eta"+s_region+IDsuffix, this_lep.Eta(), weight, -3., 3., 120);
        FillHist("dXY"+s_region+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
        FillHist("HT"+s_region+IDsuffix, HT, weight, 0., 1000., 1000);
        FillHist("MET"+s_region+IDsuffix, MET, weight, 0., 1000., 1000);
        FillHist("LT"+s_region+IDsuffix, LT, weight, 0., 1000., 1000);
        FillHist("energy"+s_region+IDsuffix, this_lep.E(), weight, 0., 500., 500);

        FillHist("HALFTEST_n_events_global"+IDsuffix+s_sample, 0., weight, 0., 1., 1);
	FillHist("HALFTEST_MET_global"+IDsuffix+s_sample, MET, weight, 0., 100., 10);
        FillHist("HALFTEST_n_jets_global"+IDsuffix+s_sample, Njets, weight, 0., 5., 5);
	FillHist("HALFTEST_HT_global"+IDsuffix+s_sample, HT, weight, 0., 150., 10);
	FillHist("Pt_eta_global"+IDsuffix+s_sample, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
	if(Njets == 0) FillHist("Pt_eta_global_JETS0"+IDsuffix+s_sample, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
        if(Njets != 0) FillHist("Pt_eta_global_JETS"+IDsuffix+s_sample, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);

        if( is_CF ){
          FillHist("Pt_eta_global_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
          if(Njets == 0) FillHist("Pt_eta_global_JETS0_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
          if(Njets != 0) FillHist("Pt_eta_global_JETS_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
          FillHist("n_events_global_CF"+IDsuffix, 0., weight, 0., 1., 1);
          FillHist("Pt_global_CF"+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
          FillHist("invPt_global_CF"+IDsuffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist("eta_global_CF"+IDsuffix, this_lep.Eta(), weight, -3., 3., 120);
          FillHist("dXY_global_CF"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
          FillHist("HT_global_CF"+IDsuffix, HT, weight, 0., 1000., 1000);
          FillHist("MET_global_CF"+IDsuffix, MET, weight, 0., 1000., 1000);
          FillHist("LT_global_CF"+IDsuffix, LT, weight, 0., 1000., 1000);
          FillHist("energy_global_CF"+IDsuffix, this_lep.E(), weight, 0., 500., 500);

          FillHist("n_events"+s_region+"_CF"+IDsuffix, 0., weight, 0., 1., 1);
          FillHist("Pt"+s_region+"_CF"+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
          FillHist("invPt"+s_region+"_CF"+IDsuffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist("eta"+s_region+"_CF"+IDsuffix, this_lep.Eta(), weight, -3., 3., 120);
          FillHist("dXY"+s_region+"_CF"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
          FillHist("HT"+s_region+"_CF"+IDsuffix, HT, weight, 0., 1000., 1000);
          FillHist("MET"+s_region+"_CF"+IDsuffix, MET, weight, 0., 1000., 1000);
          FillHist("LT"+s_region+"_CF"+IDsuffix, LT, weight, 0., 1000., 1000);
          FillHist("energy"+s_region+"_CF"+IDsuffix, this_lep.E(), weight, 0., 500., 500);

          FillHist("HALFTEST_n_events_global_CF"+IDsuffix+s_sample, 0., weight, 0., 1., 1);
          FillHist("HALFTEST_MET_global_CF"+IDsuffix+s_sample, MET, weight, 0., 100., 10);
          FillHist("HALFTEST_n_jets_global_CF"+IDsuffix+s_sample, Njets, weight, 0., 5., 5);
          FillHist("HALFTEST_HT_global_CF"+IDsuffix+s_sample, HT, weight, 0., 150., 10);
          FillHist("Pt_eta_global_CF"+IDsuffix+s_sample, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
          if(Njets == 0) FillHist("Pt_eta_global_JETS0_CF"+IDsuffix+s_sample, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
          if(Njets != 0) FillHist("Pt_eta_global_JETS_CF"+IDsuffix+s_sample, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);


          if( is_CONV0 ){
            FillHist("Pt_eta_global_CONV0_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
            if(Njets == 0) FillHist("Pt_eta_global_JETS0_CONV0_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
            if(Njets != 0) FillHist("Pt_eta_global_JETS_CONV0_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
  	    FillHist("n_events_global_CONV0_CF"+IDsuffix, 0., weight, 0., 1., 1);
            FillHist("Pt_global_CONV0_CF"+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
            FillHist("invPt_global_CONV0_CF"+IDsuffix, this_lep.Pt(), weight, 0., 0.04, 40);
            FillHist("eta_global_CONV0_CF"+IDsuffix, this_lep.Eta(), weight, -3., 3., 120);
            FillHist("dXY_global_CONV0_CF"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
            FillHist("HT_global_CONV0_CF"+IDsuffix, HT, weight, 0., 1000., 1000);
            FillHist("MET_global_CONV0_CF"+IDsuffix, MET, weight, 0., 1000., 1000);
            FillHist("LT_global_CONV0_CF"+IDsuffix, LT, weight, 0., 1000., 1000);
            FillHist("energy_global_CONV0_CF"+IDsuffix, this_lep.E(), weight, 0., 500., 500);

            FillHist("n_events"+s_region+"_CONV0_CF"+IDsuffix, 0., weight, 0., 1., 1);
            FillHist("Pt"+s_region+"_CONV0_CF"+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
	    FillHist("invPt"+s_region+"_CONV0_CF"+IDsuffix, this_lep.Pt(), weight, 0., 0.04, 40);
            FillHist("eta"+s_region+"_CONV0_CF"+IDsuffix, this_lep.Eta(), weight, -3., 3., 120);
            FillHist("dXY"+s_region+"_CONV0_CF"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
            FillHist("HT"+s_region+"_CONV0_CF"+IDsuffix, HT, weight, 0., 1000., 1000);
            FillHist("MET"+s_region+"_CONV0_CF"+IDsuffix, MET, weight, 0., 1000., 1000);
            FillHist("LT"+s_region+"_CONV0_CF"+IDsuffix, LT, weight, 0., 1000., 1000);
            FillHist("energy"+s_region+"_CONV0_CF"+IDsuffix, this_lep.E(), weight, 0., 500., 500);

            FillHist("HALFTEST_n_events_global_CONV0_CF"+IDsuffix+s_sample, 0., weight, 0., 1., 1);
            FillHist("HALFTEST_MET_global_CONV0_CF"+IDsuffix+s_sample, MET, weight, 0., 100., 10);
            FillHist("HALFTEST_n_jets_global_CONV0_CF"+IDsuffix+s_sample, Njets, weight, 0., 5., 5);
            FillHist("HALFTEST_HT_global_CONV0_CF"+IDsuffix+s_sample, HT, weight, 0., 150., 10);
            FillHist("Pt_eta_global_CONV0_CF"+IDsuffix+s_sample, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
            if(Njets == 0) FillHist("Pt_eta_global_JETS0_CONV0_CF"+IDsuffix+s_sample, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
            if(Njets != 0) FillHist("Pt_eta_global_JETS_CONV0_CF"+IDsuffix+s_sample, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);

          }
        }
      }

/*      if((electronPromptColl.size() != 0)){
        int Nel = electronPromptColl.size();
	double h_cf_rate = -999., h_cf_weight = 0.;
	double hh_cf_rate = -999., hh_cf_weight = 0.;
        for(int hhh=0; hhh<Nel; hhh++){
	  snu::KParticle h_el = electronPromptColl.at(hhh);

          h_cf_rate = Get2DCFRates(false, h_el.Pt(), fabs(h_el.Eta()), el_ID, s_sample, "");
          h_cf_weight += (h_cf_rate / (1 - h_cf_rate));

	  hh_cf_rate = Get2DCFRates(false, h_el.Pt(), fabs(h_el.Eta()), el_ID, s_sample, "");
          hh_cf_weight += (hh_cf_rate / (1 - hh_cf_rate));

	}
        FillHist("HALFTEST_n_events_global_CF"+IDsuffix+s_sample+"TO"+s_h_sample, 0., h_cf_weight, 0., 1., 1);
        FillHist("HALFTEST_HT_global_CF"+IDsuffix+s_sample+"TO"+s_h_sample, HT, h_cf_weight, 0., 150., 10);
        FillHist("HALFTEST_MET_global_CF"+IDsuffix+s_sample+"TO"+s_h_sample, MET, h_cf_weight, 0., 100., 10);
        FillHist("HALFTEST_n_jets_global_CF"+IDsuffix+s_sample+"TO"+s_h_sample, Njets, h_cf_weight, 0., 5., 5);

        FillHist("HALFTEST_n_events_global_CF"+IDsuffix+s_sample+"TO"+s_hh_sample, 0., hh_cf_weight, 0., 1., 1);
        FillHist("HALFTEST_HT_global_CF"+IDsuffix+s_sample+"TO"+s_hh_sample, HT, hh_cf_weight, 0., 150., 10);
        FillHist("HALFTEST_MET_global_CF"+IDsuffix+s_sample+"TO"+s_hh_sample, MET, hh_cf_weight, 0., 100., 10);
        FillHist("HALFTEST_n_jets_global_CF"+IDsuffix+s_sample+"TO"+s_hh_sample, Njets, hh_cf_weight, 0., 5., 5);
      }*/

      if(electronPromptColl.size() == 2){

        std::vector<snu::KTruth> truthColl;
        eventbase->GetTruthSel()->Selection(truthColl);      
        vector<int> el_index_1, el_index_2;
        el_index_1.clear(); el_index_2.clear();

        for( int i = 2 ; i < truthColl.size() ; i++ ){
          if( (truthColl.at(i).PdgId() == 11) && (truthColl.at((truthColl.at(i).IndexMother())).PdgId() == 23) ){
            el_index_1.push_back(i);
            GENFindDecayIndex( truthColl, i, el_index_1);
            break;
          }
        }
        for( int i = 2 ; i < truthColl.size() ; i++ ){
          if( (truthColl.at(i).PdgId() == -11) && (truthColl.at((truthColl.at(i).IndexMother())).PdgId() == 23) ){
            el_index_2.push_back(i);
            GENFindDecayIndex( truthColl, i, el_index_2);
            break;
          }
        }

        if( (el_index_1.size() == 0) || (el_index_2.size() == 0) ) continue;

        snu::KTruth TRUTHel[2];
        TRUTHel[0] = truthColl.at(el_index_1.back());
        TRUTHel[1] = truthColl.at(el_index_2.back());
        snu::KElectron RECOel[2];
        RECOel[0] = electronPromptColl.at(0);
        RECOel[1] = electronPromptColl.at(1);

        if(is_CF){
          if(RECOel[0].Charge() == RECOel[1].Charge()){

	    int CF_el_index = -999;
	    int mcmatch_it = 0;
	    bool CF_el_index_found = false;
	    for(int ii=0; ii<2; ii++){
	      for(int jj=0; jj<2; jj++){
	        if((TRUTHel[ii].DeltaR(RECOel[jj]) < 0.1)){
		  mcmatch_it++;
		  if(((TRUTHel[ii].PdgId() * RECOel[jj].Charge()) > 0)){
		    CF_el_index = jj;
		    if(CF_el_index_found) FillHist("[WARNING]chargeflip_doubly_found", 0., 1., 0., 1., 1);
		    CF_el_index_found = true;
		  }
	        }
	      }
	    }
	    if(mcmatch_it != 2){ FillHist("[WARNING]more_than_two_pairs_mcmatched", mcmatch_it , 1., 0., 5., 5); continue; }
	    if(!CF_el_index_found){ FillHist("[WARNING]no_chargeflip_electron_found", 0., 1., 0., 1., 1); continue; }

//            FillHist("RECO_El_energy_SS"+IDsuffix, RECOel[CF_el_index].E(), weight, 0., 500., 500);
//            FillHist("RECO_El_Pt_SS"+IDsuffix, RECOel[CF_el_index].Pt(), weight, 0., 500., 500);
	    FillHist("RECO_Z_mass_SS"+IDsuffix, (RECOel[0]+RECOel[1]).M(), weight, 70., 110., 40);
//            FillHist("TRUTH_El_energy_SS"+IDsuffix, TRUTHel[CF_el_index].E(), weight, 0., 500., 500);
//            FillHist("TRUTH_El_Pt_SS"+IDsuffix, TRUTHel[CF_el_index].Pt(), weight, 0., 500., 500);
//            FillHist("TRUTH_Z_mass_SS"+IDsuffix, (TRUTHel[0]+TRUTHel[1]).M(), weight, 70., 110., 40);
//            FillHist("RECO_div_TRUTH_El_energy_SS"+IDsuffix, (RECOel[CF_el_index].E()/TRUTHel[CF_el_index].E()), weight, 0., 2., 2000);
//            FillHist("RECO_div_TRUTH_El_Pt_SS"+IDsuffix, (RECOel[CF_el_index].Pt()/TRUTHel[CF_el_index].Pt()), weight, 0., 2., 2000);
//            FillHist("RECO_div_TRUTH_Z_mass_SS"+IDsuffix, ((RECOel[0]+RECOel[1]).M()/(TRUTHel[0]+TRUTHel[1]).M()), weight, 0., 2., 2000);

	    double CFrate_onCFel = Get2DCFRates(false, RECOel[CF_el_index].Pt(), fabs(RECOel[CF_el_index].Eta()), el_ID, "", "");
	    double cf_onCFel_weight = (CFrate_onCFel/(1-CFrate_onCFel));
	    FillHist("RECO_Z_mass_SS_CFrateToCFEl"+IDsuffix, ((RECOel[0]+RECOel[1]).M()), cf_onCFel_weight, 70., 110., 40);
	  }
        }
        else{
          if(RECOel[0].Charge() != RECOel[1].Charge()){
          
            int Q_el_index = 0;
            double Q_check = -1.;
            if(((TRUTHel[0].PdgId() * RECOel[0].Charge()) < 0) && ((TRUTHel[1].PdgId() * RECOel[1].Charge()) < 0)){
              Q_el_index = 1;
              Q_check *= -1.;
            }
            if(((TRUTHel[0].PdgId() * RECOel[1].Charge()) < 0) && ((TRUTHel[1].PdgId() * RECOel[0].Charge()) < 0)){
              Q_el_index = -1;
              Q_check *= -1.;
            }
            if(Q_check < 0){ FillHist("[WARNING]truthmatching_method_needs_check", 0., 1., 0., 1., 1); continue;}

//            FillHist("RECO_El_energy_OS"+IDsuffix, RECOel[0].E(), weight, 0., 500., 500);
//            FillHist("RECO_El_energy_OS"+IDsuffix, RECOel[1].E(), weight, 0., 500., 500);
//            FillHist("RECO_El_Pt_OS"+IDsuffix, RECOel[0].Pt(), weight, 0., 500., 500);
//            FillHist("RECO_El_Pt_OS"+IDsuffix, RECOel[1].Pt(), weight, 0., 500., 500);
            FillHist("RECO_Z_mass_OS"+IDsuffix, (RECOel[0]+RECOel[1]).M(), weight, 70., 110., 40);
//            FillHist("TRUTH_El_energy_OS"+IDsuffix, TRUTHel[0].E(), weight, 0., 500., 500);
//            FillHist("TRUTH_El_energy_OS"+IDsuffix, TRUTHel[1].E(), weight, 0., 500., 500);
//            FillHist("TRUTH_El_Pt_OS"+IDsuffix, TRUTHel[0].Pt(), weight, 0., 500., 500);
//            FillHist("TRUTH_El_Pt_OS"+IDsuffix, TRUTHel[1].Pt(), weight, 0., 500., 500);
//            FillHist("TRUTH_Z_mass_OS"+IDsuffix, (TRUTHel[0]+TRUTHel[1]).M(), weight, 70., 110., 40);
/*            if(Q_el_index>0){
	      FillHist("RECO_div_TRUTH_El_energy_OS"+IDsuffix, (RECOel[0].E()/TRUTHel[0].E()), weight, 0., 2., 2000);
              FillHist("RECO_div_TRUTH_El_energy_OS"+IDsuffix, (RECOel[1].E()/TRUTHel[1].E()), weight, 0., 2., 2000);
              FillHist("RECO_div_TRUTH_El_Pt_OS"+IDsuffix, (RECOel[0].Pt()/TRUTHel[0].Pt()), weight, 0., 2., 2000);
              FillHist("RECO_div_TRUTH_El_Pt_OS"+IDsuffix, (RECOel[1].Pt()/TRUTHel[1].Pt()), weight, 0., 2., 2000);
	    }
	    if(Q_el_index<0){
              FillHist("RECO_div_TRUTH_El_energy_OS"+IDsuffix, (RECOel[0].E()/TRUTHel[1].E()), weight, 0., 2., 2000);
              FillHist("RECO_div_TRUTH_El_energy_OS"+IDsuffix, (RECOel[1].E()/TRUTHel[0].E()), weight, 0., 2., 2000);
              FillHist("RECO_div_TRUTH_El_Pt_OS"+IDsuffix, (RECOel[0].Pt()/TRUTHel[1].Pt()), weight, 0., 2., 2000);
              FillHist("RECO_div_TRUTH_El_Pt_OS"+IDsuffix, (RECOel[1].Pt()/TRUTHel[0].Pt()), weight, 0., 2., 2000);
	    }*/
	    if(Q_el_index==0){ FillHist("[WARNING]truthmatching_method_needs_check", 0., 1., 0., 1., 1); continue;}
            FillHist("RECO_div_TRUTH_Z_mass_OS"+IDsuffix, ((RECOel[0]+RECOel[1]).M()/(TRUTHel[0]+TRUTHel[1]).M()), weight, 0., 2., 2000);

            double CFrate[2] = {0.,};
//            CFrate[0] = Get2DCFRates(false, RECOel[0].Pt(), fabs(RECOel[0].Eta()), el_ID, "", "");
//            CFrate[1] = Get2DCFRates(false, RECOel[1].Pt(), fabs(RECOel[1].Eta()), el_ID, "", "");
            CFrate[0] = GetCFRates(0, RECOel[0].Pt(), RECOel[0].Eta(), el_ID);
            CFrate[1] = GetCFRates(0, RECOel[1].Pt(), RECOel[1].Eta(), el_ID);


	    double cf_weight = (CFrate[0] / (1 - CFrate[0])) + (CFrate[1] / (1 - CFrate[1]));

            for(int E_it = 0; E_it < 51; E_it++){
              snu::KParticle SHIFT_RECOel[2];
              TString s_shift = "_SHIFT_"+TString::Itoa(E_it, 10)+"div1000";

              double shift_rate = 0.;
              shift_rate = 1.-0.1*E_it/100.;

              SHIFT_RECOel[0] = ShiftEnergy( RECOel[0], shift_rate );
              SHIFT_RECOel[1] = ShiftEnergy( RECOel[1], shift_rate );

              FillHist("RECO_Z_mass_OS"+IDsuffix+s_shift, (SHIFT_RECOel[0]+SHIFT_RECOel[1]).M(), cf_weight, 70., 110., 40);
//              FillHist("RECO_El_energy_OS"+IDsuffix+s_shift, SHIFT_RECOel[0].E(), cf_weight, 0., 500., 500);
//              FillHist("RECO_El_energy_OS"+IDsuffix+s_shift, SHIFT_RECOel[1].E(), cf_weight, 0., 500., 500);
//              FillHist("RECO_El_Pt_OS"+IDsuffix+s_shift, SHIFT_RECOel[0].Pt(), cf_weight, 0., 500., 500);
//              FillHist("RECO_El_Pt_OS"+IDsuffix+s_shift, SHIFT_RECOel[1].Pt(), cf_weight, 0., 500., 500);

            }//E shift
          }//OS charge
        }//!is_CF
      }//prompt size == 2
    }//for different tight id iteration
  }

  return;
}// End of execute event loop
  


void CFRateCalculator_syst::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void CFRateCalculator_syst::BeginCycle() throw( LQError ){
  
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

CFRateCalculator_syst::~CFRateCalculator_syst() {
  
  Message("In CFRateCalculator_syst Destructor" , INFO);
  
}


void CFRateCalculator_syst::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void CFRateCalculator_syst::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this CFRateCalculator_systCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void CFRateCalculator_syst::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


void CFRateCalculator_syst::GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index ){

  for( int i = it+1 ; i < truthColl.size(); i ++ ){
    if( truthColl.at(i).IndexMother() == it && truthColl.at(i).PdgId() == truthColl.at(it).PdgId() ){
      index.push_back(i);
    }
  }
  return;

}



void CFRateCalculator_syst::CFvalidation(void){

  bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if( !pass_trig ) return;

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_VETO", false);
  if( muonLooseColl.size() != 0) return;

  double Z_mass = 91.1876;
  double this_reliso_el = 0.6;

  TString el_ID = "ELECTRON_HN_TIGHT";
  TString IDsuffix = "_HNTIGHT";
  std::vector<snu::KElectron> electronVLooseColl = GetElectrons(false, false, "ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KJet> jetVLooseColl = GetJets("JET_HN");

  snu::KEvent event = eventbase->GetEvent();

  int N_sys = (2*2+1);
  for(int it_sys = 0; it_sys<N_sys; it_sys++){

    double METPt = event.MET();
    double METPhi = event.METPhi();

    double this_weight = weight;
    TString this_syst;

    if(it_sys==0){
      this_syst = "ElectronE_Up";//
      METPt = event.PFMETShifted(snu::KEvent::ElectronEn, snu::KEvent::up);
    }
    else if(it_sys==1){
      this_syst = "ElectronE_Down";//
      METPt = event.PFMETShifted(snu::KEvent::ElectronEn, snu::KEvent::down);
    }
    else if(it_sys==2){//
      this_syst = "Unclustered_Up";
      METPt = event.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::up);
    }
    else if(it_sys==3){//
      this_syst = "Unclustered_Down";
      METPt = event.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::down);
    }
    else if(it_sys==4){//
      this_syst = "Central";
    }
    else{
      Message("it_sys out of range!" , INFO);
      return;
    }

    std::vector<snu::KElectron> electronLooseColl, electronTightColl;
    if(this_syst == "ElectronE_Up"){
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedUp(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        double new_reliso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedUp();
        if( this_electron.Pt() >= 25. && new_reliso < this_reliso_el ) electronLooseColl.push_back( this_electron );
      }
    }
    else if(this_syst == "ElectronE_Down"){
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedDown(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        double new_reliso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedDown();
        if( this_electron.Pt() >= 25. && new_reliso < this_reliso_el ) electronLooseColl.push_back( this_electron );
      }
    }
    else{
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        if( this_electron.Pt() >= 25. && this_electron.PFRelIso(0.3) < this_reliso_el ) electronLooseColl.push_back( this_electron );
      }
    }
    for(unsigned int i=0; i<electronLooseColl.size(); i++){
      if(PassID(electronLooseColl.at(i), "ELECTRON_HN_TIGHT")) electronTightColl.push_back( electronLooseColl.at(i) );
    }

    if(METPt > 30) continue;

    TString CFsample="_powheg";

    if(electronLooseColl.size() != 2) continue;
    if(electronTightColl.size() != 2) continue;

    FillHist("[CHECK]n_of_muons"+IDsuffix+CFsample, muonLooseColl.size(), 1., 0., 5., 5);//check no muons
    FillHist("[CHECK]n_of_electrons"+IDsuffix+CFsample, electronTightColl.size(), 1., 0., 5., 5);//check two electrons

    // define leptons and give Pt, MET cuts
    snu::KParticle lep[2];
    lep[0] = electronTightColl.at(0);
    lep[1] = electronTightColl.at(1);

    bool is_SS = false;
    if( (lep[0].Charge() == lep[1].Charge()) ) is_SS = true;

    FillHist("[CHECK]MET"+IDsuffix+CFsample+"_"+this_syst, METPt, 1., 0., 100., 100);

    double cf_weight[3] = {GetCFweight(0, electronTightColl, false),GetCFweight(1, electronTightColl, false),GetCFweight(-1, electronTightColl, false)};//0 : central, 1 : up, 2 : down
    double sf_cf_weight[3] = {GetCFweight(0, electronTightColl, true),GetCFweight(1, electronTightColl, true),GetCFweight(-1, electronTightColl, true)};

    TString region = "";
    if( (fabs(lep[0].Eta()) < 0.9) )                                      region += "iB";
    else if( (fabs(lep[0].Eta()) < 1.4442) )				    region += "oB";
    else if( (fabs(lep[0].Eta()) > 1.556) && (fabs(lep[0].Eta()) < 2.5) ) region += "E";
    else continue;

    if( (fabs(lep[1].Eta()) < 0.9) )                                      region += "iB";
    else if( (fabs(lep[1].Eta()) < 1.4442) )                              region += "oB";
    else if( (fabs(lep[1].Eta()) > 1.556) && (fabs(lep[1].Eta()) < 2.5) ) region += "E";
    else continue;

    if((region == "iBE") || (region == "EiB") || (region == "oBE") || (region == "EoB"))  region = "BE";

      // reduce energy of leptons if OS because of photon radiation
      // increase energy of leptons if SS for better fitting using Gaussian
    double shiftrate=(1-0.009);
    for(int i=0; i<2; i++){
      if( !is_SS ){
        lep[i] = ShiftEnergy(lep[i], shiftrate);
        FillHist("[CHECK]lepton_E_shift_down"+IDsuffix+CFsample, lep[i].E(), 1., 0., 400., 400);
      }
      else{
	FillHist("[CHECK]lepton_E_shift_0"+IDsuffix+CFsample, lep[i].E(), 1., 0., 400., 400);
      }
    }
    snu::KParticle Z_candidate;//define Z candidate after shifting energy (SS : better fitting, energy scale up // OS : photon radiation E loss)
    Z_candidate = (lep[0] + lep[1]);
    bool Z_selection = (((Z_candidate.M() - Z_mass) < 10.) && ((Z_mass - Z_candidate.M()) < 10.));

    int Njets = jetVLooseColl.size();
    TString s_njets = "";
    if( Njets == 0 ) s_njets = "JETS0";
    if( Njets != 0 ) s_njets = "JETS";

    if( is_SS ){
      if( Z_selection ){
        FillHist("observed_Z_mass_global"+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), weight, 70., 110., 40);
        FillHist("observed_n_events_global"+IDsuffix+CFsample+"_"+this_syst, 0., weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+region+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), weight, 70., 110., 40);
        FillHist("observed_n_events_"+region+IDsuffix+CFsample+"_"+this_syst, 0., weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+s_njets+"_global"+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), weight, 70., 110., 40);
        FillHist("observed_n_events_"+s_njets+"_global"+IDsuffix+CFsample+"_"+this_syst, 0., weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+s_njets+"_"+region+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), weight, 70., 110., 40);
        FillHist("observed_n_events_"+s_njets+"_"+region+IDsuffix+CFsample+"_"+this_syst, 0., weight, 0., 1., 1);
        if(Njets == 0) FillHist("observed_n_jets_global"+IDsuffix+CFsample+"_"+this_syst, 0., weight, 0., 2., 2);
        if(Njets != 0) FillHist("observed_n_jets_global"+IDsuffix+CFsample+"_"+this_syst, 1., weight, 0., 2., 2);
        if(Njets == 1){
          FillHist("observed_Z_mass_JETS1_global"+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), weight, 70., 110., 40);
          FillHist("observed_n_events_JETS1_global"+IDsuffix+CFsample+"_"+this_syst, 0., weight, 0., 1., 1);	 
          FillHist("observed_Z_mass_JETS1_"+region+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), weight, 70., 110., 40);
          FillHist("observed_n_events_JETS1_"+region+IDsuffix+CFsample+"_"+this_syst, 0., weight, 0., 1., 1);
        }
      }// Z selection
    }// is_SS
    if( !is_SS ){
      if( Z_selection ){//Z selection after shifting down energy
	TString CF_syst ="";
        for(int iii=0; iii<3; iii++){
          if(iii==0) CF_syst = "";
	  if(iii==1) CF_syst = "_CFup";
          if(iii==2) CF_syst = "_CFdown";
          FillHist("predicted_Z_mass_global"+IDsuffix+CFsample+"_"+this_syst+CF_syst, Z_candidate.M(), cf_weight[iii], 70., 110., 40);
          FillHist("predicted_n_events_global"+IDsuffix+CFsample+"_"+this_syst+CF_syst, 0., cf_weight[iii], 0., 1., 1);
          FillHist("predicted_Z_mass_"+region+IDsuffix+CFsample+"_"+this_syst+CF_syst, Z_candidate.M(), cf_weight[iii], 70., 110., 40);
          FillHist("predicted_n_events_"+region+IDsuffix+CFsample+"_"+this_syst+CF_syst, 0., cf_weight[iii], 0., 1., 1);
          FillHist("predicted_Z_mass_"+s_njets+"_global"+IDsuffix+CFsample+"_"+this_syst+CF_syst, Z_candidate.M(), cf_weight[iii], 70., 110., 40);
          FillHist("predicted_n_events_"+s_njets+"_global"+IDsuffix+CFsample+"_"+this_syst+CF_syst, 0., cf_weight[iii], 0., 1., 1);
          FillHist("predicted_Z_mass_"+s_njets+"_"+region+IDsuffix+CFsample+"_"+this_syst+CF_syst, Z_candidate.M(), cf_weight[iii], 70., 110., 40);
          FillHist("predicted_n_events_"+s_njets+"_"+region+IDsuffix+CFsample+"_"+this_syst+CF_syst, 0., cf_weight[iii], 0., 1., 1);
          if(Njets == 0) FillHist("predicted_n_jets_global"+IDsuffix+CFsample+"_"+this_syst+CF_syst, 0., cf_weight[iii], 0., 2., 2);
          if(Njets != 0) FillHist("predicted_n_jets_global"+IDsuffix+CFsample+"_"+this_syst+CF_syst, 1., cf_weight[iii], 0., 2., 2);
          if(Njets == 1){
            FillHist("predicted_Z_mass_JETS1_global"+IDsuffix+CFsample+"_"+this_syst+CF_syst, Z_candidate.M(), cf_weight[iii], 70., 110., 40);
            FillHist("predicted_n_events_JETS1_global"+IDsuffix+CFsample+"_"+this_syst+CF_syst, 0., cf_weight[iii], 0., 1., 1);
            FillHist("predicted_Z_mass_JETS1_"+region+IDsuffix+CFsample+"_"+this_syst+CF_syst, Z_candidate.M(), cf_weight[iii], 70., 110., 40);
            FillHist("predicted_n_events_JETS1_"+region+IDsuffix+CFsample+"_"+this_syst+CF_syst, 0., cf_weight[iii], 0., 1., 1);
          }

          FillHist("predicted_Z_mass_global"+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, Z_candidate.M(), sf_cf_weight[iii], 70., 110., 40);
          FillHist("predicted_n_events_global"+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, 0., sf_cf_weight[iii], 0., 1., 1);
          FillHist("predicted_Z_mass_"+region+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, Z_candidate.M(), sf_cf_weight[iii], 70., 110., 40);
          FillHist("predicted_n_events_"+region+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, 0., sf_cf_weight[iii], 0., 1., 1);
          FillHist("predicted_Z_mass_"+s_njets+"_global"+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, Z_candidate.M(), sf_cf_weight[iii], 70., 110., 40);
          FillHist("predicted_n_events_"+s_njets+"_global"+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, 0., sf_cf_weight[iii], 0., 1., 1);
          FillHist("predicted_Z_mass_"+s_njets+"_"+region+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, Z_candidate.M(), sf_cf_weight[iii], 70., 110., 40);
          FillHist("predicted_n_events_"+s_njets+"_"+region+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, 0., sf_cf_weight[iii], 0., 1., 1);
          if(Njets == 0) FillHist("predicted_n_jets_global"+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, 0., sf_cf_weight[iii], 0., 2., 2);
          if(Njets != 0) FillHist("predicted_n_jets_global"+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, 1., sf_cf_weight[iii], 0., 2., 2);
          if(Njets == 1){
            FillHist("predicted_Z_mass_JETS1_global"+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, Z_candidate.M(), sf_cf_weight[iii], 70., 110., 40);
            FillHist("predicted_n_events_JETS1_global"+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, 0., sf_cf_weight[iii], 0., 1., 1);
            FillHist("predicted_Z_mass_JETS1_"+region+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, Z_candidate.M(), sf_cf_weight[iii], 70., 110., 40);
            FillHist("predicted_n_events_JETS1_"+region+IDsuffix+CFsample+"_SF"+"_"+this_syst+CF_syst, 0., sf_cf_weight[iii], 0., 1., 1);
          }
        }     
      } // Z selection
    } // !is_SS  
  }
  
  return;
}

snu::KParticle CFRateCalculator_syst::ShiftEnergy( snu::KParticle old_lep, double shift_rate ){

  double mass = 0.511e-3;
  snu::KParticle new_lep;
  new_lep.SetPtEtaPhiM((shift_rate*old_lep.Pt()), old_lep.Eta(), old_lep.Phi(), mass) ;
  return new_lep;


/*
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

  return new_lep;*/
}

double CFRateCalculator_syst::GetCFweight(int sys, std::vector<snu::KElectron> electrons, bool apply_sf){

  if(electrons.size() != 2) return 0.;

  snu::KParticle lep[2];
  lep[0] = electrons.at(0);
  lep[1] = electrons.at(1);

  TString el_ID = "ELECTRON_HN_TIGHT";
  double CFrate[2] = {0.,}, CFweight[2] = {0.,};
  CFrate[0] = GetCFRates(sys, lep[0].Pt(), lep[0].Eta(), el_ID);
  CFrate[1] = GetCFRates(sys, lep[1].Pt(), lep[1].Eta(), el_ID);

  CFweight[0] = CFrate[0] / (1-CFrate[0]);
  CFweight[1] = CFrate[1] / (1-CFrate[1]);

  double sf[2] = {1., 1.};
  if(apply_sf){
    for(int i=0; i<2; i++){
      if(fabs(lep[i].Eta()) < 0.9) sf[i] = 0.611579;
      else if (fabs(lep[i].Eta()) < 1.4442) sf[i] = 0.806748;
      else sf[i] = 0.803323;
    }
  }

  return (CFweight[0]*sf[0] + CFweight[1]*sf[1]);
}

double CFRateCalculator_syst::GetCFRates(int sys, double el_pt, double el_eta, TString el_ID){

  el_eta = fabs(el_eta);
  if(el_eta > 1.4442 && el_eta < 1.556) return 0.;

  double invPt = 1./el_pt;
  double a = 999., b= 999.;
  double da = 999., db = 999.;
  if(el_eta < 0.9){
    if(invPt< 0.022){
      a=(-0.00218306); da=(0.000680599);
      b=(5.79586e-05); db=(1.3127e-05);
    }
    else{
      a=(-4.69233e-05); da=(0.000235286);
      b=(1.05643e-05); db=(6.76982e-06);
    }
  }
  else if(el_eta < 1.4442){
    if(invPt < 0.006){
      a=(-0.641053); da=(0.430823);
      b=(0.00437416); db=(0.00219064);
    }
    else if(invPt< 0.021){
      a=(-0.0356571); da=(0.00374569);
      b=(0.00080504); db=(6.99218e-05);
    }
    else{
      a=(-0.00172847); da=(0.000675869);
      b=(9.17233e-05); db=(1.96677e-05);
    }
  }
  else{
    if(invPt< 0.011){
      a=(-0.398458); da=(0.0922853);
      b=(0.00623644); db=(0.000852791);
    }
    else if(invPt< 0.021){
      a=(-0.138186); da=(0.00989859);
      b=(0.00328402); db=(0.000183971);
    }
    else{
      a=(-0.0150409); da=(0.0015816);
      b=(0.000696152); db=(4.70148e-05);
    }
  }

  double rate = (a+sys*da)*invPt + (b+sys*db);
  if(rate < 0) rate = 0.;
  return rate;

}

double CFRateCalculator_syst::Get2DCFRates(bool apply_sf, double el_pt, double el_eta, TString el_ID, TString halfsample, TString CFsample){

  int N_pt = 7, N_eta = 5;

  double ptarray[7] = {20., 40., 60., 80., 100., 200., 500.};
  double etaarray[5] = {0.0, 0.9, 1.4442, 1.556, 2.5};

  if(el_pt < 20) el_pt = 21.;
  if(el_pt > 200) el_pt = 201.;

  if(el_eta > 1.4442 && el_eta < 1.556) return 0.;

  int this_pt_bin, this_eta_bin;
  for(int i=0; i<N_pt; i++){
    if( ptarray[i] <= el_pt && el_pt < ptarray[i+1] ){
      this_pt_bin = i+1;
      break;
    }
  }
  for(int i=0; i<N_eta; i++){
    if( etaarray[i] <= el_eta && el_eta < etaarray[i+1] ){
      this_eta_bin = i+1;
      break;
    }
  }

  double CFrate = -999.;
  
  if(k_sample_name.Contains("DYtoEE") || CFsample == "_powheg"){
    if(halfsample == ""){
      if(el_ID == "ELECTRON_HN_TIGHT") CFrate =  HNTIGHT_CF_hist_powheg->GetBinContent(this_eta_bin, this_pt_bin);
      if(el_ID == "ELECTRON_MVA_TIGHT") CFrate =  MVATIGHT_CF_hist_powheg->GetBinContent(this_eta_bin, this_pt_bin);
    }
    if(halfsample == "_sampleA"){
      if(el_ID == "ELECTRON_HN_TIGHT") CFrate =  HNTIGHT_CF_sampleA_hist_powheg->GetBinContent(this_eta_bin, this_pt_bin);
      if(el_ID == "ELECTRON_MVA_TIGHT") CFrate =  MVATIGHT_CF_sampleA_hist_powheg->GetBinContent(this_eta_bin, this_pt_bin);
    }
    if(halfsample == "_sampleB"){
      if(el_ID == "ELECTRON_HN_TIGHT") CFrate =  HNTIGHT_CF_sampleB_hist_powheg->GetBinContent(this_eta_bin, this_pt_bin);
      if(el_ID == "ELECTRON_MVA_TIGHT") CFrate =  MVATIGHT_CF_sampleB_hist_powheg->GetBinContent(this_eta_bin, this_pt_bin);
    }

    if(apply_sf){
      if(el_ID == "ELECTRON_HN_TIGHT"){
        if(el_eta < 1.4442) CFrate *= 0.609966;
        else CFrate *= 0.818954;
      }
      if(el_ID == "ELECTRON_MVA_TIGHT"){
        if(el_eta < 1.4442) CFrate *= 0.81493;
        else CFrate *= 1.06593;
      }
    }
  }
  if(k_sample_name.Contains("DYJets_MG") || CFsample == "_madgraph"){
    if(halfsample == ""){
      if(el_ID == "ELECTRON_HN_TIGHT") CFrate =  HNTIGHT_CF_hist_madgraph->GetBinContent(this_eta_bin, this_pt_bin);
      if(el_ID == "ELECTRON_MVA_TIGHT") CFrate =  MVATIGHT_CF_hist_madgraph->GetBinContent(this_eta_bin, this_pt_bin);
    }
    if(halfsample == "_sampleA"){
      if(el_ID == "ELECTRON_HN_TIGHT") CFrate =  HNTIGHT_CF_sampleA_hist_madgraph->GetBinContent(this_eta_bin, this_pt_bin);
      if(el_ID == "ELECTRON_MVA_TIGHT") CFrate =  MVATIGHT_CF_sampleA_hist_madgraph->GetBinContent(this_eta_bin, this_pt_bin);
    }
    if(halfsample == "_sampleB"){
      if(el_ID == "ELECTRON_HN_TIGHT") CFrate =  HNTIGHT_CF_sampleB_hist_madgraph->GetBinContent(this_eta_bin, this_pt_bin);
      if(el_ID == "ELECTRON_MVA_TIGHT") CFrate =  MVATIGHT_CF_sampleB_hist_madgraph->GetBinContent(this_eta_bin, this_pt_bin);
    }

    if(apply_sf){
      if(el_ID == "ELECTRON_HN_TIGHT"){
        if(el_eta < 1.4442) CFrate *= 0.611501;
        else CFrate *= 0.821459;
      }
      if(el_ID == "ELECTRON_MVA_TIGHT"){
        if(el_eta < 1.4442) CFrate *= 0.815791;
        else CFrate *= 1.06738;
      }
    }
  }


  return CFrate;
}

