// $Id: CFRateCalculator.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQCFRateCalculator Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "CFRateCalculator.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (CFRateCalculator);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
CFRateCalculator::CFRateCalculator() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("CFRateCalculator");
  
  Message("In CFRateCalculator constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

  TFile* file_1 = new TFile("/home/shjeon/CATanalyzer_v806/data/Fake/80X/CFRateHist_TIGHT1.root");
  TFile* file_2 = new TFile("/home/shjeon/CATanalyzer_v806/data/Fake/80X/CFRateHist_TIGHT2.root");
  TIGHT1_CF_hist = (TH2F*)file_1->Get("Pt_eta_global_CF_TIGHT1")->Clone();
  TIGHT2_CF_hist = (TH2F*)file_2->Get("Pt_eta_global_CF_TIGHT2")->Clone();

}


void CFRateCalculator::InitialiseAnalysis() throw( LQError ) {
  
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


void CFRateCalculator::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);

  if( isData ){//|| k_sample_name.Contains("DY") ){
    CFvalidation();
    if( isData ) return;
  }

  // ========== Pileup reweight ====================
  float pileup_reweight=(1.0);
  if(!k_isdata){ pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);}
  // ================================================================================

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", weight);

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                

  bool trig_pass = false;
  TString s_trigger = "";
  if(std::find(k_flags.begin(), k_flags.end(), "single") !=k_flags.end()){
    trig_pass = (PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v") || PassTrigger("HLT_Ele17_CaloIdL_GsfTrkIdVL_v"));
    s_trigger = "single";
    if(!trig_pass) return;
  }
  else if(std::find(k_flags.begin(), k_flags.end(), "double") !=k_flags.end()){
    trig_pass = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    s_trigger = "double";
    if(!trig_pass) return;
  }
  else{
    FillHist("[ERROR]trigger_not_defined", 0., 1., 0., 1., 1);
    trig_pass = (PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v") || PassTrigger("HLT_Ele17_CaloIdL_GsfTrkIdVL_v"));
    s_trigger = "single";
    if(!trig_pass) return;
  }

  weight = pileup_reweight;

  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN", 30., 2.4);
  int Njets = jetTightColl.size();
  double HT = 0.;
  for(int i=0; i<Njets; i++){
    HT += jetTightColl.at(i).Pt();
  }

  double MET = eventbase->GetEvent().MET();

//  if(MET > 30) return;

  for(int aaa=0; aaa<2; aaa++){

    std::vector<snu::KElectron> electronTightColl, electronPromptColl;
    electronPromptColl.clear();

    TString IDsuffix = "";

    if(aaa == 0){
      electronTightColl = GetElectrons(true,false,"ELECTRON_HN_TIGHT");
      IDsuffix = "_TIGHT1";
    }
    if(aaa == 1){
      electronTightColl = GetElectrons(true,false,"ELECTRON_HN_TIGHT2");
      IDsuffix = "_TIGHT2";
    }

    double LT = 0.;
    for(int i=0; i<electronTightColl.size(); i++){
      snu::KElectron this_lep;
      this_lep = electronTightColl.at(i);
      if((this_lep.Pt() < 25)) continue;
      if((fabs(this_lep.Eta()) < 1.556) && (fabs(this_lep.Eta()) > 1.4442)) continue;
      if((this_lep.MCIsPrompt())){
        electronPromptColl.push_back(this_lep);
        LT += this_lep.Pt();
      }
    }

    FillHist("SIZE_ELECTRON_PROMPTCOLL"+IDsuffix, electronPromptColl.size(), weight, 0., 5., 5);

    if(electronPromptColl.size() == 0) return;

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
      else{
        if(s_trigger == "single") continue;
      }

      if( (MCIsCF(this_lep)) ) is_CF = true;
      if( !(this_lep.MCIsFromConversion()) ) is_CONV0 = true;

      FillHist("Pt_eta_global"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
      if(Njets == 0) FillHist("Pt_eta_global_JETS0"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
      if(Njets > 1) FillHist("Pt_eta_global_JETS"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
      FillHist("n_events_global"+IDsuffix, 0., weight, 0., 1., 1);
      FillHist("dXY_global"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
      FillHist("eta_global"+IDsuffix, this_lep.Eta(), weight, -3., 3., 120);
      FillHist("invPt_global"+IDsuffix, (1./this_lep.Pt()), weight, 0., 0.05, 25);
      FillHist("HT_global"+IDsuffix, HT, weight, 0., 1000., 1000);
      FillHist("MET_global"+IDsuffix, MET, weight, 0., 1000., 1000);
      FillHist("LT_global"+IDsuffix, LT, weight, 0., 1000., 1000);

      if( is_CF ){
        FillHist("Pt_eta_global_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
        if(Njets == 0) FillHist("Pt_eta_global_JETS0_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
        if(Njets > 1) FillHist("Pt_eta_global_JETS_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
        FillHist("n_events_global_CF"+IDsuffix, 0., weight, 0., 1., 1);
        FillHist("dXY_global_CF"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
        FillHist("eta_global_CF"+IDsuffix, this_lep.Eta(), weight, -3., 3., 120);
        FillHist("invPt_global_CF"+IDsuffix, (1./this_lep.Pt()), weight, 0., 0.05, 25);
        FillHist("HT_global_CF"+IDsuffix, HT, weight, 0., 1000., 1000);
        FillHist("MET_global_CF"+IDsuffix, MET, weight, 0., 1000., 1000);
        FillHist("LT_global_CF"+IDsuffix, LT, weight, 0., 1000., 1000);

        if( is_CONV0 ){
          FillHist("Pt_eta_global_CONV0_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
          if(Njets == 0) FillHist("Pt_eta_global_JETS0_CONV0_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
          if(Njets > 1) FillHist("Pt_eta_global_JETS_CONV0_CF"+IDsuffix, fabs(this_lep.Eta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
  	  FillHist("n_events_global_CONV0_CF"+IDsuffix, 0., weight, 0., 1., 1);
          FillHist("dXY_global_CONV0_CF"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
          FillHist("eta_global_CONV0_CF"+IDsuffix, this_lep.Eta(), weight, -3., 3., 120);
	  FillHist("invPt_global_CONV0_CF"+IDsuffix, (1./this_lep.Pt()), weight, 0., 0.05, 25);
          FillHist("HT_global_CONV0_CF"+IDsuffix, HT, weight, 0., 1000., 1000);
          FillHist("MET_global_CONV0_CF"+IDsuffix, MET, weight, 0., 1000., 1000);
          FillHist("LT_global_CONV0_CF"+IDsuffix, LT, weight, 0., 1000., 1000);
        }
      }

      TString suffix = "";
      if( is_region == 1 ) suffix = "region1";
      else if( is_region == 2 ) suffix = "region2";
      else if( is_region == 3 ) suffix = "region3";
      else FillHist("[WARNING]suffix_not_defined", 0., 1., 0., 1., 1);

      for( int j=1; j<4; j++ ){

        if( is_region == j ){
          FillHist("n_events_"+suffix+IDsuffix, 0., weight, 0., 1., 1);
          FillHist("eta_"+suffix+IDsuffix, (this_lep.Eta()), weight, -3., 3., 30);
          FillHist("invPt_"+suffix+IDsuffix, (1./this_lep.Pt()), weight, 0., 0.05, 25);

          if( is_CF ){
            FillHist("n_events_"+suffix+"_CF"+IDsuffix, 0., weight, 0., 1., 1);
            FillHist("eta_"+suffix+"_CF"+IDsuffix, (this_lep.Eta()), weight, -3., 3., 30);
            FillHist("invPt_"+suffix+"_CF"+IDsuffix, (1./this_lep.Pt()), weight, 0., 0.05, 25);

	    if( is_CONV0 ){
              FillHist("n_events_"+suffix+"_CONV0_CF"+IDsuffix, 0., weight, 0., 1., 1);
              FillHist("eta_"+suffix+"_CONV0_CF"+IDsuffix, (this_lep.Eta()), weight, -3., 3., 30);
              FillHist("invPt_"+suffix+"_CONV0_CF"+IDsuffix, (1./this_lep.Pt()), weight, 0., 0.05, 25);

	    }
  	  }
        }
      }//Fill hists for different regions
    }
    if( is_CF ){
      if(electronPromptColl.size() == 2){
        if(electronPromptColl.at(0).Charge() == electronPromptColl.at(1).Charge()){

          std::vector<snu::KTruth> truthColl;
          eventbase->GetTruthSel()->Selection(truthColl);

          vector<int> el_index_1, el_index_2;
          el_index_1.clear(); el_index_2.clear();

          for( int i = 2 ; i < truthColl.size() ; i++ ){
            if( (truthColl.at(i).PdgId()) == 11 ){
              el_index_1.push_back(i);
              GENFindDecayIndex( truthColl, i, el_index_1);
              break;
            }
	  }
          for( int i = 2 ; i < truthColl.size() ; i++ ){
	    if( (truthColl.at(i).PdgId()) == -11 ){
	      el_index_2.push_back(i);
	      GENFindDecayIndex( truthColl, i, el_index_2);
  	      break;
  	    }
	  }

          snu::KTruth TRUTHel[2];
          TRUTHel[0] = truthColl.at(el_index_1.back());
          TRUTHel[1] = truthColl.at(el_index_2.back());

          if(!((truthColl.at(TRUTHel[0].IndexMother()).PdgId() == 23) && (truthColl.at(TRUTHel[1].IndexMother()).PdgId() == 23))) continue;

 	  snu::KParticle RECOel[2];
	  RECOel[0] = electronPromptColl.at(0);
	  RECOel[1] = electronPromptColl.at(1);

	  //PdgId = 11 & charge -1 -> not CF
	  //PdgId = 11 & charge 1 -> CF
	  int CF_el_index = -999;
	  double check = -1.;
	  snu::KParticle CFTRUTHel, CFRECOel;
	  if(((TRUTHel[0].PdgId() * RECOel[0].Charge()) < 0) && ((TRUTHel[1].PdgId() * RECOel[1].Charge()) > 0)){
	    CF_el_index = 1;
	    check *= -1.;
  	  }
	  if(((TRUTHel[0].PdgId() * RECOel[0].Charge()) > 0) && ((TRUTHel[1].PdgId() * RECOel[1].Charge()) < 0)){
            CF_el_index = 0;
	    check *= -1.;
          }
	  if(check < 0){ FillHist("[Warning]charge_flip_method_needs_check", 0., 1., 0., 1., 1); return;}

	  FillHist("RECO_energy_CF"+IDsuffix, TRUTHel[CF_el_index].E(), weight, 0., 500., 500);
	  FillHist("TRUTH_energy_CF"+IDsuffix, RECOel[CF_el_index].E(), weight, 0., 500., 500);
	  FillHist("RECO_div_TRUTH_energy_CF"+IDsuffix, (RECOel[CF_el_index].E()/TRUTHel[CF_el_index].E()), weight, 0., 2., 2000);
        }//Samesign electrons
      }//2 prompt electrons
    }// cf is true

  }//for different tight id iteration


  return;
}// End of execute event loop
  


void CFRateCalculator::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void CFRateCalculator::BeginCycle() throw( LQError ){
  
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

CFRateCalculator::~CFRateCalculator() {
  
  Message("In CFRateCalculator Destructor" , INFO);
  
}


void CFRateCalculator::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void CFRateCalculator::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this CFRateCalculatorCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void CFRateCalculator::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


void CFRateCalculator::GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index ){

  for( int i = it+1 ; i < truthColl.size(); i ++ ){
    if( truthColl.at(i).IndexMother() == it && truthColl.at(i).PdgId() == truthColl.at(it).PdgId() ){
      index.push_back(i);
    }
  }
  return;

}



void CFRateCalculator::CFvalidation(void){

  bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if( !pass_trig ) return;

  // define electron and muon Colls
  std::vector<snu::KElectron> electronTightColl;

  electronTightColl =  GetElectrons(true,false,"ELECTRON_HN_TIGHT");

  if(electronTightColl.size() != 2) return;

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE", false);
  if( muonLooseColl.size() != 0) return;

  // define leptons and give Pt, MET cuts
  snu::KParticle lep[2];
  lep[0] = electronTightColl.at(0);
  lep[1] = electronTightColl.at(1);

  bool is_SS = false;
  if( (lep[0].Charge() == lep[1].Charge()) ) is_SS = true;

  if( lep[0].Pt() < 25 || lep[1].Pt() < 15 ) return;

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();

  if(METPt > 30) return;

  FillHist("n_muons", muonLooseColl.size(), 1., 0., 5., 5);
  FillHist("n_electrons", electronTightColl.size(), 1., 0., 5., 5);

  // reduce energy of leptons if OS because of photon radiation
  bool energy_reduced = false;
  if( !is_SS ){
    for(int i=0; i<2; i++){
      FillHist("lepton_originalE", lep[i].E(), 1., 0., 400., 400);
      
      lep[i] = ReduceEnergy(lep[i]);
      energy_reduced = true;
    }
  }

  if( energy_reduced ){
    FillHist("lepton_reducedE", lep[0].E(), 1., 0., 400., 400);
    FillHist("lepton_reducedE", lep[1].E(), 1., 0., 400., 400);
  }

  // =============================== Get global Scale Factor from observed/predicted
  // define dilepton mass
  snu::KParticle Z_candidate;
  Z_candidate = (lep[0] + lep[1]);
  bool Z_selection = (fabs(Z_candidate.M() - 90) < 15);

  // observed same sign events
  if( is_SS ){
    if( Z_selection ){
      FillHist("observed_Z_mass_global", Z_candidate.M(), 1., 70., 110., 40);
      FillHist("observed_n_events_global", 0., 1., 0., 1., 1);
    }
  }

  double CFrate[2] = {-999.,};
  CFrate[0] = Get2DCFRates(false, lep[0].Pt(), fabs(lep[0].Eta())); // do not apply sf since this codes are for getting the sf
  CFrate[1] = Get2DCFRates(false, lep[1].Pt(), fabs(lep[1].Eta()));

  double cf_weight = (CFrate[0] / (1 - CFrate[0])) + (CFrate[1] / (1 - CFrate[1]));

  // predicteded same sign events
  if( !is_SS ){
    if(  Z_selection ){
      FillHist("predicted_Z_mass_global", Z_candidate.M(), cf_weight, 70., 110., 40);
      FillHist("predicted_n_events_global", 0., cf_weight, 0., 1., 1);
    }
  }
  // ================================================ end of global SF obtaining

  if(k_sample_name.Contains("DY")){

    if( (fabs(Z_candidate.M() - 90) < 10) ){
      if( is_SS ){
        FillHist("DY_observed_Z_mass_global", Z_candidate.M(), 1., 70., 110., 40);
        FillHist("DY_observed_n_events_global", 0., 1., 0., 1., 1);
      }
      if( !is_SS ){
        FillHist("DY_predicted_Z_mass_global", Z_candidate.M(), cf_weight, 70., 110., 40);
        FillHist("DY_predicted_n_events_global", 0., cf_weight, 0., 1., 1);
      }
    }

    return;
  }


  bool is_region[2][2] = {{false,},};
  if( (fabs(lep[0].Eta()) < 1.4442) )                              	is_region[0][0] = true;
  else if( (fabs(lep[0].Eta()) > 1.556) && (fabs(lep[0].Eta()) < 2.5) ) is_region[1][0] = true;

  if( (fabs(lep[1].Eta()) < 1.4442) )                            	is_region[0][1] = true;
  else if( (fabs(lep[1].Eta()) > 1.556) && (fabs(lep[1].Eta()) < 2.5) ) is_region[1][1] = true;

 
  TString region = "none";
  if((is_region[0][0] && is_region[0][1])) 	region = "BB";
  else if((is_region[1][0] && is_region[1][1])) region = "EE";
  else						region = "BE";
  
//  cout << lep[0].Charge() << " " << lep[1].Charge() << endl;

  if( region == "BB" ){
    if( Z_selection ){
      if( is_SS ){
	FillHist("observed_n_BB_EE_BE", 0., 1., 0., 3., 3);
        FillHist("observed_Z_mass_BB", Z_candidate.M(), 1., 70., 110., 40);
      }
      if( !is_SS ){
        FillHist("predicted_n_BB_EE_BE", 0., cf_weight, 0., 3., 3);
        FillHist("predicted_Z_mass_BB", Z_candidate.M(), cf_weight, 70., 110., 40);
      }
    }
  }
  else if( region == "EE" ){
    if( Z_selection ){
      if( is_SS ){
	FillHist("observed_n_BB_EE_BE", 1., 1., 0., 3., 3);
        FillHist("observed_Z_mass_EE", Z_candidate.M(), 1., 70., 110., 40);
      }
      if( !is_SS ){
        FillHist("predicted_n_BB_EE_BE", 1., cf_weight, 0., 3., 3);
        FillHist("predicted_Z_mass_EE", Z_candidate.M(), cf_weight, 70., 110., 40);
      }
    }
  }
  else if( region == "BE" ){
    if( Z_selection ){
      if( is_SS ){
        FillHist("observed_n_BB_EE_BE", 2., 1., 0., 3., 3);
        FillHist("observed_Z_mass_BE", Z_candidate.M(), 1., 70., 110., 40);
      }
      if( !is_SS ){
        FillHist("predicted_n_BB_EE_BE", 2., cf_weight, 0., 3., 3);
        FillHist("predicted_Z_mass_BE", Z_candidate.M(), cf_weight, 70., 110., 40);
      }
    }
  }
  else{
    FillHist("[WARNING]region_not_defined", 0., 1., 0., 1., 1);
    return;
  }


  CFrate[0] = Get2DCFRates(true, lep[0].Pt(), fabs(lep[0].Eta())); // do not apply sf since this codes are for getting the sf
  CFrate[1] = Get2DCFRates(true, lep[1].Pt(), fabs(lep[1].Eta()));

  cf_weight = (CFrate[0] / (1 - CFrate[0])) + (CFrate[1] / (1 - CFrate[1]));

  if( region == "BE" ){
    if( Z_selection ){
      if( is_SS ){
        FillHist("SF_observed_Z_mass_BE", Z_candidate.M(), 1., 70., 110., 40);
        FillHist("SF_observed_n_BE", 0., 1., 0., 1., 1);
      }
      if( !is_SS ){
	FillHist("SF_predicted_Z_mass_BE", Z_candidate.M(), cf_weight, 70., 110., 40);
	FillHist("SF_predicted_n_BE", 0., cf_weight, 0., 1., 1);
      }
    }
  }

  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN", 30., 2.4);
  if( Z_selection ){
    if( is_SS ){
      FillHist("SF_observed_Z_mass_global", Z_candidate.M(), 1., 70., 110., 40);
      FillHist("SF_observed_n_global", 0., 1., 0., 1., 1);
    }
    if( !is_SS ){
      FillHist("SF_predicted_Z_mass_global", Z_candidate.M(), cf_weight, 70., 110., 40);
      FillHist("SF_predicted_n_global", 0., cf_weight, 0., 1., 1);
    }
  }
  if( jetTightColl.size() == 0 ){
    if( Z_selection ){
      if( is_SS ){
        FillHist("SF_0jet_observed_Z_mass_global", Z_candidate.M(), 1., 70., 110., 40);
        FillHist("SF_0jet_observed_n_global", 0., 1., 0., 1., 1);
      }
      if( !is_SS ){
        FillHist("SF_0jet_predicted_Z_mass_global", Z_candidate.M(), cf_weight, 70., 110., 40);
        FillHist("SF_0jet_predicted_n_global", 0., cf_weight, 0., 1., 1);
      }
    }
  }
  if( jetTightColl.size() > 0 ){
    if( Z_selection ){
      if( is_SS ){
        FillHist("SF_1jet_observed_Z_mass_global", Z_candidate.M(), 1., 70., 110., 40);
        FillHist("SF_1jet_observed_n_global", 0., 1., 0., 1., 1);
      }
      if( !is_SS ){
        FillHist("SF_1jet_predicted_Z_mass_global", Z_candidate.M(), cf_weight, 70., 110., 40);
        FillHist("SF_1jet_predicted_n_global", 0., cf_weight, 0., 1., 1);
      }
    }
  }
  if( jetTightColl.size() > 1 ){
    if( Z_selection ){
      if( is_SS ){
        FillHist("SF_Njet_observed_Z_mass_global", Z_candidate.M(), 1., 70., 110., 40);
        FillHist("SF_Njet_observed_n_global", 0., 1., 0., 1., 1);
      }
      if( !is_SS ){
        FillHist("SF_Njet_predicted_Z_mass_global", Z_candidate.M(), cf_weight, 70., 110., 40);
        FillHist("SF_Njet_predicted_n_global", 0., cf_weight, 0., 1., 1);
      }
    }
  }

  return;

}

snu::KParticle CFRateCalculator::ReduceEnergy( snu::KParticle old_lep ){

   double new_E, new_px, new_py, new_pz;
   double mass;
   double new_psum, old_psum;
   new_E = old_lep.E() * 0.96;
   mass = old_lep.M();

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


double CFRateCalculator::Get2DCFRates(bool apply_sf, double el_pt, double el_eta){

  int N_pt = 7, N_eta = 5;

  double ptarray[7] = {20., 40., 60., 80., 100., 200., 500.};
  double etaarray[5] = {0.0, 0.9, 1.4442, 1.556, 2.5};

  if(el_pt < 20) el_pt = 21.;
  if(el_pt > 200) el_pt = 200.;

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
  CFrate =  TIGHT1_CF_hist->GetBinContent(this_eta_bin, this_pt_bin);

  if(apply_sf){
    if(el_eta < 1.4442) CFrate *= 0.6180;
    else CFrate *= 0.8518;
  }
  return CFrate;
}
