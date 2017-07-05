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

/*  TFile* file_madgraph = new TFile("/home/shjeon/CATanalyzer_v807/data/Fake/80X/ChargeFlip_madgraph_v807.root");
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
  MVATIGHT_CF_sampleB_hist_powheg  = (TH2F*)file_powheg ->Get("Pt_eta_global_CF_MVATIGHT_PU_sampleB")->Clone();*/

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
   
  if(!PassMETFilter()) return;     /// Initial event cuts :
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex

  // ========== Pileup reweight ====================
  float pileup_reweight=(1.0);
  if(!k_isdata){
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    weight *= pileup_reweight;
  }
  // ================================================================================



  if(isData){CFvalidation(); return;}
 
  bool trig_pass = (PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v") || PassTrigger("HLT_Ele17_CaloIdL_GsfTrkIdVL_v"));
  TString s_trigger = "single";
  if(!trig_pass) return;

  TString sample_suffix = "";
  if(k_sample_name.Contains("DYtoEE")) sample_suffix = "_powheg";
  else if(k_sample_name.Contains("DYJets_MG")) sample_suffix = "_madgraph";
  else return;

  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
  int Njets = jetTightColl.size();
  double HT = 0.;
  for(int i=0; i<Njets; i++){
    HT += jetTightColl.at(i).Pt();
  }

  double MET = eventbase->GetEvent().MET();

  for(int aa=0; aa<1; aa++){
    if(aa==0) weight = pileup_reweight;
    else return;

    for(int aaa=0; aaa<1; aaa++){

      std::vector<snu::KElectron> electronTightColl, electronPromptColl;
      electronPromptColl.clear();

      TString IDsuffix = "";
      TString el_ID = "";

      if(aaa == 0){
        el_ID = "ELECTRON_HN_TIGHTv4";
        IDsuffix = "_HNTIGHT";
      }
      else return;

      electronTightColl = GetElectrons(true, false, el_ID);

      if(aa == 0){IDsuffix += "_PU";}

      IDsuffix += sample_suffix;

      double LT = 0.;
      for(int i=0; i<electronTightColl.size(); i++){
        snu::KElectron this_lep;
        this_lep = electronTightColl.at(i);
        if((this_lep.Pt() < 25)) continue;

        FillHist("[CHECK]lepton_eta", this_lep.Eta(), weight, -5., 5., 1000);
	FillHist("[CHECK]lepton_SCeta", this_lep.SCEta(), weight, -5., 5., 1000);

        if((fabs(this_lep.SCEta()) < 1.556) && (fabs(this_lep.SCEta()) > 1.4442)) continue;
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
        if( (fabs(this_lep.SCEta()) < 0.9) )                                          is_region = 1;
        else if( (fabs(this_lep.SCEta()) < 1.4442) )                                  is_region = 2;
        else if( (fabs(this_lep.SCEta()) > 1.556) && (fabs(this_lep.SCEta()) < 2.5) ) is_region = 3;
	else continue;

        TString s_region = "";
        if( is_region == 1 ) s_region = "_region1";
        else if( is_region == 2 ) s_region = "_region2";
        else if( is_region == 3 ) s_region = "_region3";
        else FillHist("[WARNING]suffix_not_defined", 0., 1., 0., 1., 1);

        if( (MCIsCF(this_lep)) ) is_CF = true;
        if( !(this_lep.MCIsFromConversion()) ) is_CONV0 = true;

        FillHist("Pt_eta_global"+IDsuffix, fabs(this_lep.SCEta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
        if(Njets == 0) FillHist("Pt_eta_global_JETS0"+IDsuffix, fabs(this_lep.SCEta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
        if(Njets != 0) FillHist("Pt_eta_global_JETS"+IDsuffix, fabs(this_lep.SCEta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
        FillHist("n_events_global"+IDsuffix, 0., weight, 0., 1., 1);
        FillHist("Pt_global"+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
        FillHist("invPt_global"+IDsuffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist("eta_global"+IDsuffix, this_lep.SCEta(), weight, -3., 3., 120);
        FillHist("dXY_global"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
        FillHist("HT_global"+IDsuffix, HT, weight, 0., 1000., 1000);
        FillHist("MET_global"+IDsuffix, MET, weight, 0., 1000., 1000);
        FillHist("LT_global"+IDsuffix, LT, weight, 0., 1000., 1000);
        FillHist("energy_global"+IDsuffix, this_lep.E(), weight, 0., 500., 500);

        FillHist("n_events"+s_region+IDsuffix, 0., weight, 0., 1., 1);
        FillHist("Pt"+s_region+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
        FillHist("invPt"+s_region+IDsuffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist("eta"+s_region+IDsuffix, this_lep.SCEta(), weight, -3., 3., 120);
        FillHist("dXY"+s_region+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
        FillHist("HT"+s_region+IDsuffix, HT, weight, 0., 1000., 1000);
        FillHist("MET"+s_region+IDsuffix, MET, weight, 0., 1000., 1000);
        FillHist("LT"+s_region+IDsuffix, LT, weight, 0., 1000., 1000);
        FillHist("energy"+s_region+IDsuffix, this_lep.E(), weight, 0., 500., 500);

        if( is_CF ){
          FillHist("Pt_eta_global_CF"+IDsuffix, fabs(this_lep.SCEta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
          if(Njets == 0) FillHist("Pt_eta_global_JETS0_CF"+IDsuffix, fabs(this_lep.SCEta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
          if(Njets != 0) FillHist("Pt_eta_global_JETS_CF"+IDsuffix, fabs(this_lep.SCEta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
          FillHist("n_events_global_CF"+IDsuffix, 0., weight, 0., 1., 1);
          FillHist("Pt_global_CF"+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
          FillHist("invPt_global_CF"+IDsuffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist("eta_global_CF"+IDsuffix, this_lep.SCEta(), weight, -3., 3., 120);
          FillHist("dXY_global_CF"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
          FillHist("HT_global_CF"+IDsuffix, HT, weight, 0., 1000., 1000);
          FillHist("MET_global_CF"+IDsuffix, MET, weight, 0., 1000., 1000);
          FillHist("LT_global_CF"+IDsuffix, LT, weight, 0., 1000., 1000);
          FillHist("energy_global_CF"+IDsuffix, this_lep.E(), weight, 0., 500., 500);

          FillHist("n_events"+s_region+"_CF"+IDsuffix, 0., weight, 0., 1., 1);
          FillHist("Pt"+s_region+"_CF"+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
          FillHist("invPt"+s_region+"_CF"+IDsuffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist("eta"+s_region+"_CF"+IDsuffix, this_lep.SCEta(), weight, -3., 3., 120);
          FillHist("dXY"+s_region+"_CF"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
          FillHist("HT"+s_region+"_CF"+IDsuffix, HT, weight, 0., 1000., 1000);
          FillHist("MET"+s_region+"_CF"+IDsuffix, MET, weight, 0., 1000., 1000);
          FillHist("LT"+s_region+"_CF"+IDsuffix, LT, weight, 0., 1000., 1000);
          FillHist("energy"+s_region+"_CF"+IDsuffix, this_lep.E(), weight, 0., 500., 500);

          if( is_CONV0 ){
            FillHist("Pt_eta_global_CONV0_CF"+IDsuffix, fabs(this_lep.SCEta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
            if(Njets == 0) FillHist("Pt_eta_global_JETS0_CONV0_CF"+IDsuffix, fabs(this_lep.SCEta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
            if(Njets != 0) FillHist("Pt_eta_global_JETS_CONV0_CF"+IDsuffix, fabs(this_lep.SCEta()), this_lep.Pt(), weight, etaarray, 4, ptarray, 6);
  	    FillHist("n_events_global_CONV0_CF"+IDsuffix, 0., weight, 0., 1., 1);
            FillHist("Pt_global_CONV0_CF"+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
            FillHist("invPt_global_CONV0_CF"+IDsuffix, this_lep.Pt(), weight, 0., 0.04, 40);
            FillHist("eta_global_CONV0_CF"+IDsuffix, this_lep.SCEta(), weight, -3., 3., 120);
            FillHist("dXY_global_CONV0_CF"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
            FillHist("HT_global_CONV0_CF"+IDsuffix, HT, weight, 0., 1000., 1000);
            FillHist("MET_global_CONV0_CF"+IDsuffix, MET, weight, 0., 1000., 1000);
            FillHist("LT_global_CONV0_CF"+IDsuffix, LT, weight, 0., 1000., 1000);
            FillHist("energy_global_CONV0_CF"+IDsuffix, this_lep.E(), weight, 0., 500., 500);

            FillHist("n_events"+s_region+"_CONV0_CF"+IDsuffix, 0., weight, 0., 1., 1);
            FillHist("Pt"+s_region+"_CONV0_CF"+IDsuffix, this_lep.Pt(), weight, 0., 500, 500);
	    FillHist("invPt"+s_region+"_CONV0_CF"+IDsuffix, this_lep.Pt(), weight, 0., 0.04, 40);
            FillHist("eta"+s_region+"_CONV0_CF"+IDsuffix, this_lep.SCEta(), weight, -3., 3., 120);
            FillHist("dXY"+s_region+"_CONV0_CF"+IDsuffix, fabs(this_lep.dxy()), weight, 0., 0.02, 100);
            FillHist("HT"+s_region+"_CONV0_CF"+IDsuffix, HT, weight, 0., 1000., 1000);
            FillHist("MET"+s_region+"_CONV0_CF"+IDsuffix, MET, weight, 0., 1000., 1000);
            FillHist("LT"+s_region+"_CONV0_CF"+IDsuffix, LT, weight, 0., 1000., 1000);
            FillHist("energy"+s_region+"_CONV0_CF"+IDsuffix, this_lep.E(), weight, 0., 500., 500);

          }
        }
      }

      if(electronPromptColl.size() == 2){

	CF_MCClosure(electronPromptColl);

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

	    FillHist("RECO_Z_mass_SS"+IDsuffix, (RECOel[0]+RECOel[1]).M(), weight, 70., 110., 40);

/*	    double CFrate_onCFel = Get2DCFRates(false, RECOel[CF_el_index].Pt(), fabs(RECOel[CF_el_index].SCEta()), el_ID, "", "");
	    double cf_onCFel_weight = (CFrate_onCFel/(1-CFrate_onCFel));
	    FillHist("RECO_Z_mass_SS_CFrateToCFEl"+IDsuffix, ((RECOel[0]+RECOel[1]).M()), cf_onCFel_weight, 70., 110., 40);*/
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

            FillHist("RECO_Z_mass_OS"+IDsuffix, (RECOel[0]+RECOel[1]).M(), weight, 70., 110., 40);
	    if(Q_el_index==0){ FillHist("[WARNING]truthmatching_method_needs_check", 0., 1., 0., 1., 1); continue;}
            FillHist("RECO_div_TRUTH_Z_mass_OS"+IDsuffix, ((RECOel[0]+RECOel[1]).M()/(TRUTHel[0]+TRUTHel[1]).M()), weight, 0., 2., 2000);

            double CFrate[2] = {0.,};
            CFrate[0] = GetCFRates(0, RECOel[0].Pt(), RECOel[0].SCEta(), el_ID);
            CFrate[1] = GetCFRates(0, RECOel[1].Pt(), RECOel[1].SCEta(), el_ID);

	    double cf_weight = (CFrate[0] / (1 - CFrate[0])) + (CFrate[1] / (1 - CFrate[1]));

            for(int E_it = 0; E_it < 51; E_it++){
              snu::KParticle SHIFT_RECOel[2];
              TString s_shift = "_SHIFT_"+TString::Itoa(E_it, 10)+"div1000";

              double shift_rate = 0.;
              shift_rate = 1.-0.1*E_it/100.;

              SHIFT_RECOel[0] = ShiftEnergy( RECOel[0], shift_rate );
              SHIFT_RECOel[1] = ShiftEnergy( RECOel[1], shift_rate );

              FillHist("RECO_Z_mass_OS"+IDsuffix+s_shift, (SHIFT_RECOel[0]+SHIFT_RECOel[1]).M(), cf_weight, 70., 110., 40);

            }//E shift
          }//OS charge
        }//!is_CF
      }//prompt size == 2
    }//for different tight id iteration
  }

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


void CFRateCalculator::CF_MCClosure(std::vector<snu::KElectron> electron){

  if(electron.size() != 2 ) return;
  double Z_mass = 91.1876; 
 
  TString Zsuffix[4] = {"_narrowZ", "", "_wideZ", "_verywideZ"};
  double Zwidth[4] = {-2., 0., 5., 10.};
  snu::KElectron lep[2];

  lep[0] = electron.at(0);
  lep[1] = electron.at(1);

  bool is_SS = false;
  if((lep[0].Charge() == lep[1].Charge())) is_SS = true;

  if(!is_SS){
    for(int i=0;i<2;i++){
      lep[i] = ShiftEnergy( lep[i], (1-0.015) );
    }
  }

  snu::KParticle Z_candidate = lep[0]+lep[1]; 

  for(int i=0; i<4; i++){
    double cf_weight = GetCFweight(0, electron, false, Zsuffix[i]);
    bool Z_selection = (((Z_candidate.M() - Z_mass) < (10.+Zwidth[i])) && ((Z_mass - Z_candidate.M()) < (10.+Zwidth[i])));

    if(!is_SS){
      if(Z_selection){

	FillHist("MCClosure_Z_mass_pre"+Zsuffix[i], Z_candidate.M(), cf_weight, 60., 120., 60);
	FillHist("MCClosure_leading_lep_Pt_pre"+Zsuffix[i], lep[0].Pt(), cf_weight, 0., 200., 200);
        FillHist("MCClosure_subleading_lep_Pt_pre"+Zsuffix[i], lep[1].Pt(), cf_weight, 0., 200., 200);
	FillHist("MCClosure_n_events_pre"+Zsuffix[i], 0., cf_weight, 0., 1., 1);

      }
    }
    if(is_SS){
      if(Z_selection){

        FillHist("MCClosure_Z_mass_obs"+Zsuffix[i], Z_candidate.M(), 1, 60., 120., 60);
        FillHist("MCClosure_leading_lep_Pt_obs"+Zsuffix[i], lep[0].Pt(), 1, 0., 200., 200);
        FillHist("MCClosure_subleading_lep_Pt_obs"+Zsuffix[i], lep[1].Pt(), 1, 0., 200., 200);
        FillHist("MCClosure_n_events_obs"+Zsuffix[i], 0., 1, 0., 1., 1);

      }
    }
  }
  return;

}

void CFRateCalculator::CFvalidation(void){

  bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if( !pass_trig ) return;

  double Z_mass = 91.1876;

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  if(METPt > 30) return;

  for(int bbb=0; bbb<1; bbb++){
    TString CFsample="";
    if(bbb==0){
      CFsample="_powheg";
    }
    if(bbb==1){
      CFsample="_madgraph";
    }

    for(int aaa=0; aaa<1; aaa++){
      TString el_ID="";
      TString el_looseID="";
      TString IDsuffix="";
      if(aaa==0){
        el_ID = "HN";
      }
      if(aaa==1){
        el_ID = "MVA";
      }

      IDsuffix = "_"+el_ID+"TIGHT";
      el_looseID = "ELECTRON_"+el_ID+"_FAKELOOSE";
      el_ID = "ELECTRON_"+el_ID+"_TIGHTv4";

      std::vector<snu::KElectron> electronLooseColl = GetElectrons(true, false, el_looseID);
      std::vector<snu::KElectron> electronTightColl = GetElectrons(true, false, el_ID);
      if(electronLooseColl.size() != 2) continue;
      if(electronTightColl.size() != 2) continue;
      std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_VETO", false);
      if( muonLooseColl.size() != 0) continue;

      FillHist("[CHECK]n_of_muons"+IDsuffix+CFsample, muonLooseColl.size(), 1., 0., 5., 5);//check no muons
      FillHist("[CHECK]n_of_electrons"+IDsuffix+CFsample, electronTightColl.size(), 1., 0., 5., 5);//check two electrons

      // define leptons and give Pt, MET cuts
      snu::KElectron lep[2];
      lep[0] = electronTightColl.at(0);
      lep[1] = electronTightColl.at(1);

      bool is_SS = false;
      if( (lep[0].Charge() == lep[1].Charge()) ) is_SS = true;

      if( lep[0].Pt() < 25 || lep[1].Pt() < 25 ) continue;

      if(bbb==0 && (lep[0].Charge() != lep[1].Charge())){
	FillHist("[CHECK]mass_of_OSpair"+IDsuffix, (lep[0]+lep[1]).M(), weight, 80., 100., 20);
	if(((lep[0]+lep[1]).M() < 100) && ((lep[0]+lep[1]).M() > 80)) FillHist("[CHECK]number_of_OSpair"+IDsuffix, 0, weight, 0., 1., 1);
      }

      FillHist("[CHECK]MET"+IDsuffix+CFsample, METPt, 1., 0., 100., 100);

      TString region = "";
      if( (fabs(lep[0].SCEta()) < 0.9) )                                      region += "iB";
      else if( (fabs(lep[0].SCEta()) < 1.4442) )				    region += "oB";
      else if( (fabs(lep[0].SCEta()) > 1.556) && (fabs(lep[0].SCEta()) < 2.5) ) region += "E";
      else continue;

      if( (fabs(lep[1].SCEta()) < 0.9) )                                      region += "iB";
      else if( (fabs(lep[1].SCEta()) < 1.4442) )                              region += "oB";
      else if( (fabs(lep[1].SCEta()) > 1.556) && (fabs(lep[1].SCEta()) < 2.5) ) region += "E";
      else continue;

      if((region == "iBE") || (region == "EiB") || (region == "oBE") || (region == "EoB")) region = "BE";
      if((region == "iBiB") || (region == "iBoB") || (region == "oBiB") || (region == "oBoB")) region = "BB";
      // reduce energy of leptons if OS because of photon radiation
      // increase energy of leptons if SS for better fitting using Gaussian

      double shiftrate=-99.;
      if(CFsample=="_powheg"){
        if(el_ID=="ELECTRON_HN_TIGHTv4") shiftrate = (1.-0.015);
        if(el_ID=="ELECTRON_MVA_TIGHT") shiftrate = (1.-0.014);
      }
      if(CFsample=="_madgraph"){
        if(el_ID=="ELECTRON_HN_TIGHTv4") shiftrate = (1.-0.021);
        if(el_ID=="ELECTRON_MVA_TIGHT") shiftrate = (1.-0.016);
      }
      for(int i=0; i<2; i++){
        if( !is_SS ){
          lep[i] = ShiftEnergy(lep[i], shiftrate);
          FillHist("[CHECK]lepton_E_shift_down"+IDsuffix+CFsample, lep[i].E(), 1., 0., 400., 400);
        }
      }
      snu::KParticle Z_candidate;//define Z candidate after shifting energy (SS : better fitting, energy scale up // OS : photon radiation E loss)
      Z_candidate = (lep[0] + lep[1]);

      TString Zsuffix[4] = {"_narrowZ", "", "_wideZ", "_verywideZ"};
      double Zwidth[4] = {-2., 0., 5., 10.};
   
      std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
      int Njets = jetTightColl.size();
      TString s_njets = "";
      if( Njets == 0 ) s_njets = "JETS0";
      if( Njets != 0 ) s_njets = "JETS";

      double cf_weight = GetCFweight(0, electronTightColl, false, "");

      for(int Z_it = 0; Z_it < 4; Z_it ++){
        double sf_cf_weight = GetCFweight(0, electronTightColl, true, Zsuffix[Z_it]);
        bool Z_selection = (((Z_candidate.M() - Z_mass) < (10.+Zwidth[Z_it])) && ((Z_mass - Z_candidate.M()) < (10.+Zwidth[Z_it])));
        if( is_SS ){
          if(Z_it == 0 && aaa == 0 && bbb == 0){
            snu::KElectron shift_el[2];
	    snu::KParticle Z_candidate_shift;
            for(int shift_it=0; shift_it<2; shift_it++){
              shift_el[shift_it] = ShiftEnergy( lep[shift_it], 1/(1-0.015) );
            }
	    Z_candidate_shift = (shift_el[0]+shift_el[1]);
	    if(((Z_candidate_shift.M() - Z_mass) < 20) && ((Z_mass - Z_candidate_shift.M()) < 20)){
              FillHist("FIT_observed_Z_mass_global"+IDsuffix, Z_candidate_shift.M(), weight, 60., 120., 60);
              FillHist("FIT_observed_n_events_global"+IDsuffix, 0., weight, 0., 1., 1);
              FillHist("FIT_observed_Z_mass_"+region+IDsuffix, Z_candidate_shift.M(), weight, 60., 120., 60);
	      FillHist("FIT_observed_n_events_"+region+IDsuffix, 0., weight, 0., 1., 1);
	      if(Njets == 0){
                FillHist("FIT_observed_Z_mass_global_JETS0"+IDsuffix, Z_candidate_shift.M(), weight, 60., 120., 60);
                FillHist("FIT_observed_n_events_global_JETS0"+IDsuffix, 0., weight, 0., 1., 1);
	      }
              if(Njets == 1){
                FillHist("FIT_observed_Z_mass_global_JETS1"+IDsuffix, Z_candidate_shift.M(), weight, 60., 120., 60);
                FillHist("FIT_observed_n_events_global_JETS1"+IDsuffix, 0., weight, 0., 1., 1);
              }
	    }
	  }
          if( Z_selection ){
            FillHist("observed_Z_mass_global"+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), weight, 60., 120., 60);
            FillHist("observed_n_events_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., weight, 0., 1., 1);
            FillHist("observed_Z_mass_"+region+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), weight, 60., 120., 60);
            FillHist("observed_n_events_"+region+IDsuffix+CFsample+Zsuffix[Z_it], 0., weight, 0., 1., 1);
            FillHist("observed_Z_mass_"+s_njets+"_global"+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), weight, 60., 120., 60);
            FillHist("observed_n_events_"+s_njets+"_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., weight, 0., 1., 1);
            FillHist("observed_Z_mass_"+s_njets+"_"+region+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), weight, 60., 120., 60);
            FillHist("observed_n_events_"+s_njets+"_"+region+IDsuffix+CFsample+Zsuffix[Z_it], 0., weight, 0., 1., 1);
            if(Njets == 0) FillHist("observed_n_jets_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., weight, 0., 2., 2);
            if(Njets != 0) FillHist("observed_n_jets_global"+IDsuffix+CFsample+Zsuffix[Z_it], 1., weight, 0., 2., 2);
            if(Njets == 1){
              FillHist("observed_Z_mass_JETS1_global"+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), weight, 60., 120., 60);
              FillHist("observed_n_events_JETS1_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., weight, 0., 1., 1);
	      FillHist("observed_Z_mass_JETS1_"+region+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), weight, 60., 120., 60);
              FillHist("observed_n_events_JETS1_"+region+IDsuffix+CFsample+Zsuffix[Z_it], 0., weight, 0., 1., 1);
  	    }
          }// Z selection
        }// is_SS
        if( !is_SS ){
          if( Z_selection ){//Z selection after shifting down energy
            FillHist("predicted_Z_mass_global"+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), cf_weight, 60., 120., 60);
            FillHist("predicted_n_events_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., cf_weight, 0., 1., 1);
            FillHist("predicted_Z_mass_"+region+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), cf_weight, 60., 120., 60);
            FillHist("predicted_n_events_"+region+IDsuffix+CFsample+Zsuffix[Z_it], 0., cf_weight, 0., 1., 1);
            FillHist("predicted_Z_mass_"+s_njets+"_global"+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), cf_weight, 60., 120., 60);
            FillHist("predicted_n_events_"+s_njets+"_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., cf_weight, 0., 1., 1);
            FillHist("predicted_Z_mass_"+s_njets+"_"+region+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), cf_weight, 60., 120., 60);
            FillHist("predicted_n_events_"+s_njets+"_"+region+IDsuffix+CFsample+Zsuffix[Z_it], 0., cf_weight, 0., 1., 1);
            if(Njets == 0) FillHist("predicted_n_jets_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., cf_weight, 0., 2., 2);
            if(Njets != 0) FillHist("predicted_n_jets_global"+IDsuffix+CFsample+Zsuffix[Z_it], 1., cf_weight, 0., 2., 2);
            if(Njets == 1){
              FillHist("predicted_Z_mass_JETS1_global"+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), cf_weight, 60., 120., 60);
              FillHist("predicted_n_events_JETS1_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., cf_weight, 0., 1., 1);
              FillHist("predicted_Z_mass_JETS1_"+region+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), cf_weight, 60., 120., 60);
              FillHist("predicted_n_events_JETS1_"+region+IDsuffix+CFsample+Zsuffix[Z_it], 0., cf_weight, 0., 1., 1);
            }

            FillHist("predicted_Z_mass_global"+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], Z_candidate.M(), sf_cf_weight, 60., 120., 60);
            FillHist("predicted_n_events_global"+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], 0., sf_cf_weight, 0., 1., 1);
            FillHist("predicted_Z_mass_"+region+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], Z_candidate.M(), sf_cf_weight, 60., 120., 60);
            FillHist("predicted_n_events_"+region+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], 0., sf_cf_weight, 0., 1., 1);
            FillHist("predicted_Z_mass_"+s_njets+"_global"+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], Z_candidate.M(), sf_cf_weight, 60., 120., 60);
            FillHist("predicted_n_events_"+s_njets+"_global"+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], 0., sf_cf_weight, 0., 1., 1);
            FillHist("predicted_Z_mass_"+s_njets+"_"+region+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], Z_candidate.M(), sf_cf_weight, 60., 120., 60);
            FillHist("predicted_n_events_"+s_njets+"_"+region+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], 0., sf_cf_weight, 0., 1., 1);
            if(Njets == 0) FillHist("predicted_n_jets_global"+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], 0., sf_cf_weight, 0., 2., 2);
            if(Njets != 0) FillHist("predicted_n_jets_global"+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], 1., sf_cf_weight, 0., 2., 2);
            if(Njets == 1){
              FillHist("predicted_Z_mass_JETS1_global"+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], Z_candidate.M(), sf_cf_weight, 60., 120., 60);
              FillHist("predicted_n_events_JETS1_global"+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], 0., sf_cf_weight, 0., 1., 1);
              FillHist("predicted_Z_mass_JETS1_"+region+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], Z_candidate.M(), sf_cf_weight, 60., 120., 60);
              FillHist("predicted_n_events_JETS1_"+region+IDsuffix+CFsample+"_SF"+Zsuffix[Z_it], 0., sf_cf_weight, 0., 1., 1);
            }
        
          } // Z selection
        } // !is_SS
      } // Z width iteration
    }
  }
  return;
}

snu::KElectron CFRateCalculator::ShiftEnergy( snu::KElectron old_lep, double shift_rate ){

  double mass = 0.511e-3;
  snu::KElectron new_lep;
  new_lep.SetPtEtaPhiM((shift_rate*old_lep.Pt()), old_lep.Eta(), old_lep.Phi(), mass) ;
  return new_lep;

}


double CFRateCalculator::GetCFweight(int sys, std::vector<snu::KElectron> electrons, bool apply_sf, TString Zwidth){

  if(electrons.size() != 2) return 0.;

  snu::KElectron lep[2];
  lep[0] = electrons.at(0);
  lep[1] = electrons.at(1);

  TString el_ID = "ELECTRON_HN_TIGHTv4";
  double CFrate[2] = {0.,}, CFweight[2] = {0.,};
  CFrate[0] = GetCFRates(0, lep[0].Pt(), lep[0].SCEta(), el_ID);
  CFrate[1] = GetCFRates(0, lep[1].Pt(), lep[1].SCEta(), el_ID);

  CFweight[0] = CFrate[0] / (1-CFrate[0]);
  CFweight[1] = CFrate[1] / (1-CFrate[1]);
  double sf[2] = {1., 1.};
  if(apply_sf){
    if(Zwidth == "_narrowZ"){
      for(int i=0; i<2; i++){
        if (fabs(lep[i].SCEta()) < 1.4442) sf[i] = 0.765865231;
        else sf[i] = 0.856420327;
      }
    }
    if(Zwidth == ""){
      for(int i=0; i<2; i++){
        if (fabs(lep[i].SCEta()) < 1.4442) sf[i] = 0.75362822;
        else sf[i] = 0.821682654;
      }
    }
    if(Zwidth == "_wideZ"){
      for(int i=0; i<2; i++){
        if (fabs(lep[i].SCEta()) < 1.4442) sf[i] = 0.759713941;
        else sf[i] = 0.784052036;
      }
    }
    if(Zwidth == "_verywideZ"){
      for(int i=0; i<2; i++){
        if (fabs(lep[i].SCEta()) < 1.4442) sf[i] = 0.723099195;
        else sf[i] = 0.757193848;
      }
    }
  }

  return (CFweight[0]*sf[0] + CFweight[1]*sf[1]);
}


double CFRateCalculator::GetCFRates(int sys, double el_pt, double el_eta, TString el_ID){

  el_eta = fabs(el_eta);
  if(el_eta > 1.4442 && el_eta < 1.556) return 0.;

  double invPt = 1./el_pt;
  double a = 999., b= 999.;
  double da = 999., db = 999.;
  if(el_eta < 0.9){
    if(invPt< 0.023){
      a=(-0.00138635);
      b=(4.35054e-05);
    }
    else{
      a=(0.00114356); 
      b=(-1.55941e-05);
    }
  }
  else if(el_eta < 1.4442){
    if(invPt < 0.016){
      a=(-0.0369937);
      b=(0.000797434);
    }
    else if(invPt < 0.024){
      a=(-0.0159017);
      b=(0.00046038);
    }
    else{
      a=(-0.00214657); 
      b=(0.000147245);
    }
  }
  else{
    if(invPt< 0.012){
      a=(-0.4293);
      b=(0.00641511);
    }
    else if(invPt< 0.020){
      a=(-0.104796);
      b=(0.00256146);
    }
    else{
      a=(-0.0161499);
      b=(0.00076872); 
    }
  }

  double rate = (a+sys*da)*invPt + (b+sys*db);
  if(rate < 0) rate = 0.;
  return rate;

}

