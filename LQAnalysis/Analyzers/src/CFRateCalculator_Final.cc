// $Id: CFRateCalculator_Final.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQCFRateCalculator_Final Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "CFRateCalculator_Final.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (CFRateCalculator_Final);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
CFRateCalculator_Final::CFRateCalculator_Final() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("CFRateCalculator_Final");
  
  Message("In CFRateCalculator_Final constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void CFRateCalculator_Final::InitialiseAnalysis() throw( LQError ) {
  
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


void CFRateCalculator_Final::ExecuteEvents()throw( LQError ){

  TString h_sample_suffix = "_sampleA";
  if((eventbase->GetEvent().EventNumber())%2==0) h_sample_suffix = "_sampleB";

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  if(!PassMETFilter()) return;     /// Initial event cuts :
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex

  if(!k_isdata){
    float pileup_reweight=(1.0);
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    weight = pileup_reweight;
  }//Pileup reweights


  std::vector<TString> el_id;
  el_id.push_back("ELECTRON_HN_TIGHTv4");
  el_id.push_back("ELECTRON_HN_FAKELOOSEv7_5");

  double opt_shiftrate[2] = {0.0,};
  opt_shiftrate[0] = (1.-13./1000.);
  opt_shiftrate[1] = (1.-11./1000.);

  //This is for CFRATE and CLOSURE and HALFTEST
  for(int it_id = 0; it_id < el_id.size(); it_id++){
    if(!k_isdata){

      std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
      int Njets = jetTightColl.size();
      double HT = 0.;
      for(int i=0; i<Njets; i++){
        HT += jetTightColl.at(i).Pt();
      }

      double MET = eventbase->GetEvent().MET();
      double LT = 0.;

      std::vector<snu::KElectron> electronTightColl, electronPromptColl;
      electronTightColl = GetElectrons(true, false, el_id.at(it_id));
      electronPromptColl.clear();

      for(int i=0; i<electronTightColl.size(); i++){
        snu::KElectron this_lep=electronTightColl.at(i);
        if(this_lep.Pt() < 25) continue;
        if(fabs(this_lep.SCEta()) > 1.4442 && fabs(this_lep.SCEta()) < 1.5560) continue;
        if(this_lep.MCIsPrompt()){
          electronPromptColl.push_back(this_lep);
          LT += this_lep.Pt();
        }
      }

      FillHist(el_id.at(it_id)+"_[CHECK]size_of_electron_promptcoll", electronPromptColl.size(), weight, 0., 5., 5);

      bool one_is_CF = false;
      for(int i=0; i<electronPromptColl.size(); i++){
        int is_region = 0;
        bool is_CF = false;
        bool is_NOTCONV = false;

        snu::KElectron this_lep=electronPromptColl.at(i);

        if( (fabs(this_lep.SCEta()) < 0.8) )                                          is_region = 1;
        else if( (fabs(this_lep.SCEta()) < 1.4442) )                                  is_region = 2;
        else if( (fabs(this_lep.SCEta()) > 1.5560) )                                  is_region = 3;
        else continue;

        TString s_region = "";
        if( is_region == 1 ) s_region = "Region1";
        else if( is_region == 2 ) s_region = "Region2";
        else s_region = "Region3";
        FillHist(el_id.at(it_id)+"_[CHECK]region_of_electron_eta", is_region, weight, 0., 4., 4);

        if( (MCIsCF(this_lep)) ){ is_CF = true; one_is_CF = true; }
        if( !(this_lep.MCIsFromConversion()) ) is_NOTCONV = true;

        //Draw historgrams for each regions, halfsamples
        FillHist(el_id.at(it_id)+"_CFRATE_Global_N_Events", 0., weight, 0., 1., 1);
        FillHist(el_id.at(it_id)+"_CFRATE_Global_Pt", this_lep.Pt(), weight, 0., 500., 500);
        FillHist(el_id.at(it_id)+"_CFRATE_Global_invPt", 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist(el_id.at(it_id)+"_CFRATE_Global_SCEta", this_lep.SCEta(), weight, -3., 3., 60);
        FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_N_Events", 0., weight, 0., 1., 1);
        FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_Pt", this_lep.Pt(), weight, 0., 500., 500);
        FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_invPt", 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_SCEta", this_lep.SCEta(), weight, -3., 3., 60);
        FillHist(el_id.at(it_id)+"_CFRATE_Global_N_Events"+h_sample_suffix, 0., weight, 0., 1., 1);
        FillHist(el_id.at(it_id)+"_CFRATE_Global_Pt"+h_sample_suffix, this_lep.Pt(), weight, 0., 500., 500);
        FillHist(el_id.at(it_id)+"_CFRATE_Global_invPt"+h_sample_suffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist(el_id.at(it_id)+"_CFRATE_Global_SCEta"+h_sample_suffix, this_lep.SCEta(), weight, -3., 3., 60);
        FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_N_Events"+h_sample_suffix, 0., weight, 0., 1., 1);
        FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_Pt"+h_sample_suffix, this_lep.Pt(), weight, 0., 500., 500);
        FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_invPt"+h_sample_suffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_SCEta"+h_sample_suffix, this_lep.SCEta(), weight, -3., 3., 60);

        //Draw histograms only for half sampletest
        FillHist(el_id.at(it_id)+"_HALFTEST_Global_MET"+h_sample_suffix, MET, 1., 0., 100., 5);
        FillHist(el_id.at(it_id)+"_HALFTEST_Global_METsqdivST"+h_sample_suffix, MET*MET/(MET+LT+HT), 1., 0., 50., 5);
        FillHist(el_id.at(it_id)+"_HALFTEST_Global_N_Jets"+h_sample_suffix, Njets, 1., 0., 5., 5);
        FillHist(el_id.at(it_id)+"_HALFTEST_Global_HT"+h_sample_suffix, HT, 1., 0., 200., 5);
        double halftestrate = GetCFRates((this_lep.Pt()), (this_lep.SCEta()), el_id.at(it_id), true);
        double halftestweight = halftestrate/(1.-halftestrate);
        FillHist(el_id.at(it_id)+"_HALFTEST_Global_MET"+h_sample_suffix+"_CFpredicted", MET, halftestweight, 0., 100., 5);
        FillHist(el_id.at(it_id)+"_HALFTEST_Global_METsqdivST"+h_sample_suffix+"_CFpredicted", MET*MET/(MET+LT+HT), halftestweight, 0., 50., 5);
        FillHist(el_id.at(it_id)+"_HALFTEST_Global_N_Jets"+h_sample_suffix+"_CFpredicted", Njets, halftestweight, 0., 5., 5);
        FillHist(el_id.at(it_id)+"_HALFTEST_Global_HT"+h_sample_suffix+"_CFpredicted", HT, halftestweight, 0., 200., 5);

        if(is_CF){
          FillHist(el_id.at(it_id)+"_CFRATE_Global_N_Events_CF", 0., weight, 0., 1., 1);
          FillHist(el_id.at(it_id)+"_CFRATE_Global_Pt_CF", this_lep.Pt(), weight, 0., 500., 500);
          FillHist(el_id.at(it_id)+"_CFRATE_Global_invPt_CF", 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist(el_id.at(it_id)+"_CFRATE_Global_SCEta_CF", this_lep.SCEta(), weight, -3., 3., 60);
          FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_N_Events_CF", 0., weight, 0., 1., 1);
          FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_Pt_CF", this_lep.Pt(), weight, 0., 500., 500);
          FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_invPt_CF", 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_SCEta_CF", this_lep.SCEta(), weight, -3., 3., 60);
          FillHist(el_id.at(it_id)+"_CFRATE_Global_N_Events"+h_sample_suffix+"_CF", 0., weight, 0., 1., 1);
          FillHist(el_id.at(it_id)+"_CFRATE_Global_Pt"+h_sample_suffix+"_CF", this_lep.Pt(), weight, 0., 500., 500);
          FillHist(el_id.at(it_id)+"_CFRATE_Global_invPt"+h_sample_suffix+"_CF", 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist(el_id.at(it_id)+"_CFRATE_Global_SCEta"+h_sample_suffix+"_CF", this_lep.SCEta(), weight, -3., 3., 60);
          FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_N_Events"+h_sample_suffix+"_CF", 0., weight, 0., 1., 1);
          FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_Pt"+h_sample_suffix+"_CF", this_lep.Pt(), weight, 0., 500., 500);
          FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_invPt"+h_sample_suffix+"_CF", 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_SCEta"+h_sample_suffix+"_CF", this_lep.SCEta(), weight, -3., 3., 60);
          FillHist(el_id.at(it_id)+"_HALFTEST_Global_MET"+h_sample_suffix+"_CFobserved", MET, 1., 0., 100., 5);
          FillHist(el_id.at(it_id)+"_HALFTEST_Global_METsqdivST"+h_sample_suffix+"_CFobserved", MET*MET/(MET+LT+HT), 1., 0., 50., 5);
          FillHist(el_id.at(it_id)+"_HALFTEST_Global_N_Jets"+h_sample_suffix+"_CFobserved", Njets, 1., 0., 5., 5);
          FillHist(el_id.at(it_id)+"_HALFTEST_Global_HT"+h_sample_suffix+"_CFobserved", HT, 1., 0., 200., 5);

          if(is_NOTCONV){
            FillHist(el_id.at(it_id)+"_CFRATE_Global_N_Events_CF_NOTCONV", 0., weight, 0., 1., 1);
            FillHist(el_id.at(it_id)+"_CFRATE_Global_Pt_CF_NOTCONV", this_lep.Pt(), weight, 0., 500., 500);
            FillHist(el_id.at(it_id)+"_CFRATE_Global_invPt_CF_NOTCONV", 1./this_lep.Pt(), weight, 0., 0.04, 40);
            FillHist(el_id.at(it_id)+"_CFRATE_Global_SCEta_CF_NOTCONV", this_lep.SCEta(), weight, -3., 3., 60);
            FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_N_Events_CF_NOTCONV", 0., weight, 0., 1., 1);
            FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_Pt_CF_NOTCONV", this_lep.Pt(), weight, 0., 500., 500);
            FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_invPt_CF_NOTCONV", 1./this_lep.Pt(), weight, 0., 0.04, 40);
            FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_SCEta_CF_NOTCONV", this_lep.SCEta(), weight, -3., 3., 60);
            FillHist(el_id.at(it_id)+"_CFRATE_Global_N_Events"+h_sample_suffix+"_CF_NOTCONV", 0., weight, 0., 1., 1);
            FillHist(el_id.at(it_id)+"_CFRATE_Global_Pt"+h_sample_suffix+"_CF_NOTCONV", this_lep.Pt(), weight, 0., 500., 500);
            FillHist(el_id.at(it_id)+"_CFRATE_Global_invPt"+h_sample_suffix+"_CF_NOTCONV", 1./this_lep.Pt(), weight, 0., 0.04, 40);
            FillHist(el_id.at(it_id)+"_CFRATE_Global_SCEta"+h_sample_suffix+"_CF_NOTCONV", this_lep.SCEta(), weight, -3., 3., 60);
            FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_N_Events"+h_sample_suffix+"_CF_NOTCONV", 0., weight, 0., 1., 1);
            FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_Pt"+h_sample_suffix+"_CF_NOTCONV", this_lep.Pt(), weight, 0., 500., 500);
            FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_invPt"+h_sample_suffix+"_CF_NOTCONV", 1./this_lep.Pt(), weight, 0., 0.04, 40);
            FillHist(el_id.at(it_id)+"_CFRATE_"+s_region+"_SCEta"+h_sample_suffix+"_CF_NOTCONV", this_lep.SCEta(), weight, -3., 3., 60);

          }//Not conversion charge flips
        }//Charge flips
      }//Iterate all prompt electrons

      if(electronPromptColl.size() == 2){
        snu::KElectron this_lep[2];
        this_lep[0]=electronPromptColl.at(0);
        this_lep[1]=electronPromptColl.at(1);

        bool is_Z_loose = (fabs((this_lep[0] + this_lep[1]).M() - 91.1876) < 20.);//Wide range setting (use 15. ultimately) for shiftrate calculation
        bool is_Z_tight = (fabs((this_lep[0] + this_lep[1]).M() - 91.1876) < 15.);
        bool is_SS = (this_lep[0].Charge() == this_lep[1].Charge()) && one_is_CF;

        if((this_lep[0].Charge() == this_lep[1].Charge()) && !one_is_CF) FillHist(el_id.at(it_id)+"_[CHECK]Samesign_electrons_but_not_chargeflip", 0., 1., 0., 1., 1);

        if(is_Z_loose){
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
          if( !(el_index_1.size() == 0) && !(el_index_2.size() == 0) ){
            snu::KTruth truth_lep[2];
            truth_lep[0] = truthColl.at(el_index_1.at(0));
            truth_lep[1] = truthColl.at(el_index_2.at(0));
            std::vector<int> truth_MCmatched_index, reco_MCmatched_index;
            for(int i=0; i<2; i++){
              for(int ii=0; ii<2; ii++){
                if(truth_lep[i].DeltaR(this_lep[ii]) <0.10){
                  truth_MCmatched_index.push_back(i);
                  reco_MCmatched_index.push_back(ii);
                }
              }
            }

            bool well_matched = false;
            if(truth_MCmatched_index.size() == 2){
              if(((truth_MCmatched_index.at(0)+truth_MCmatched_index.at(1) == 1) || (reco_MCmatched_index.at(0)+reco_MCmatched_index.at(1) == 1))){
                well_matched = true;
              }
            }
            if(!well_matched){
              if(truth_MCmatched_index.size() != 2){
                FillHist(el_id.at(it_id)+"_[CHECK]MCmatched_leptons_not_2", truth_MCmatched_index.size(), weight, 0., 4., 4);
              }
              else{
                if((truth_MCmatched_index.at(0)+truth_MCmatched_index.at(1) != 1) || (reco_MCmatched_index.at(0)+reco_MCmatched_index.at(1) != 1)){
                  FillHist(el_id.at(it_id)+"_[CHECK]MCmatched_leptons_doubly_found", 0., weight, 0., 1., 1);
                }
              }
            }//Not well matched
            else{
              if(is_SS){
                if(is_Z_tight){
                  int reco_cf_index = 999, truth_cf_index = 999; bool cf_found = false;
                  if(truth_lep[truth_MCmatched_index.at(0)].PdgId() * this_lep[reco_MCmatched_index.at(0)].Charge() > 0){
                    reco_cf_index = reco_MCmatched_index.at(0);
                    truth_cf_index = truth_MCmatched_index.at(0);
                    cf_found = true;
                  }
                  if(truth_lep[truth_MCmatched_index.at(1)].PdgId() * this_lep[reco_MCmatched_index.at(1)].Charge() > 0){
                    reco_cf_index = reco_MCmatched_index.at(1);
                    truth_cf_index = truth_MCmatched_index.at(1);
                    if(cf_found) FillHist(el_id.at(it_id)+"_[CHECK]Chargeflipped_electron_doubly_found", 0., 1., 0., 1., 1);
                    cf_found = true;
                  }
                  if(!cf_found) FillHist(el_id.at(it_id)+"_[CHECK]Chargeflipped_electron_not_found", 0., 1., 0., 1., 1);
                  else{
                    FillHist(el_id.at(it_id)+"_SHIFTRATE_Electron_Pt_RecoLevel", this_lep[reco_cf_index].Pt(), 1., 0., 500., 500);
                    FillHist(el_id.at(it_id)+"_SHIFTRATE_Electron_Pt_TruthLevel", truth_lep[truth_cf_index].Pt(), 1., 0., 500., 500);
                    FillHist(el_id.at(it_id)+"_SHIFTRATE_Electron_Pt_Ratio", (this_lep[reco_cf_index].Pt()/truth_lep[truth_cf_index].Pt()), 1., 0., 2., 200);
                    if(!MCIsCF(this_lep[reco_cf_index])) FillHist(el_id.at(it_id)+"_[CHECK]Chargeflipped_electron_not_found_GetType", this_lep[reco_cf_index].GetType(), 1., 0., 50., 50);
                  }
                  FillHist(el_id.at(it_id)+"_SHIFTRATE_Zcandidate_Mass_CFobserved", (this_lep[0]+this_lep[1]).M(), 1., (91.1876-40.), (91.1876+40.), 32);
                  FillHist(el_id.at(it_id)+"_CLOSURE_Zcandidate_Mass_CFobserved", (this_lep[0]+this_lep[1]).M(), 1., 70., 110., 40);
                  FillHist(el_id.at(it_id)+"_CLOSURE_LeadingLepton_Pt_CFobserved", this_lep[0].Pt(), 1., 0., 80., 80);
                  FillHist(el_id.at(it_id)+"_CLOSURE_SubLeadingLepton_Pt_CFobserved", this_lep[1].Pt(), 1., 0., 80., 80);
                  FillHist(el_id.at(it_id)+"_CLOSURE_MET_CFobserved", MET, 1., 0., 80., 80);
                  FillHist(el_id.at(it_id)+"_CLOSURE_N_Events_CFobserved", 0., 1., 0., 1., 1);
                }//Pass tight Z requirements for SHIFTRATE and CLOSURE
              }//Is SS dielectron
              if(!is_SS){
                double this_weight=GetCFweight(electronPromptColl, false, el_id.at(it_id), false);
                bool is_shifted_Z_tight = false;
                for(int it_shift=0; it_shift<51; it_shift++){
                  snu::KElectron shifted_lep[2];
                  TString s_shift = "_Shifted_"+TString::Itoa(it_shift, 10)+"div1000";
                  double shift_rate = 1.-0.001*it_shift;
                  shifted_lep[0] = ShiftEnergy(this_lep[0], shift_rate);
                  shifted_lep[1] = ShiftEnergy(this_lep[1], shift_rate);
                  is_shifted_Z_tight = (fabs((shifted_lep[0] + shifted_lep[1]).M() - 91.1876) < 15.);
                  if(is_shifted_Z_tight){
                    FillHist(el_id.at(it_id)+"_SHIFTRATE_Zcandidate_Mass_CFpredicted"+s_shift, (shifted_lep[0]+shifted_lep[1]).M(), this_weight, (91.1876-40.), (91.1876+40.), 32);
                  }//Pass tight Z requirements after shift for SHIFTRATE
                }//Iterate shiftrates for optimal value

                snu::KElectron this_shifted_lep[2];
                this_shifted_lep[0] = ShiftEnergy(this_lep[0], opt_shiftrate[it_id]);
                this_shifted_lep[1] = ShiftEnergy(this_lep[1], opt_shiftrate[it_id]);
 
                is_shifted_Z_tight = ((fabs((this_shifted_lep[0] + this_shifted_lep[1]).M() - 91.1876) < 15.) && (this_shifted_lep[1].Pt() > 25));
                if(is_shifted_Z_tight){
                  FillHist(el_id.at(it_id)+"_CLOSURE_Zcandidate_Mass_CFpredicted", (this_shifted_lep[0]+this_shifted_lep[1]).M(), this_weight, 70., 110., 40);
                  FillHist(el_id.at(it_id)+"_CLOSURE_LeadingLepton_Pt_CFpredicted", this_shifted_lep[0].Pt(), this_weight, 0., 80., 80);
                  FillHist(el_id.at(it_id)+"_CLOSURE_SubLeadingLepton_Pt_CFpredicted", this_shifted_lep[1].Pt(), this_weight, 0., 80., 80);
                  FillHist(el_id.at(it_id)+"_CLOSURE_MET_CFpredicted", MET, this_weight, 0., 80., 80);
                  FillHist(el_id.at(it_id)+"_CLOSURE_N_Events_CFpredicted", 0., this_weight, 0., 1., 1);
                }//Pass tight Z requirements after shift for CLOSURE
              }//Is OS dielectron
            }//Well matched
          }//Find two truth matched electrons
        }//Pass loose Z requirements for SHIFTRATE and CLOSURE without any charge requirements
      }//Has two prompt electrons
    }//Is not data

    //This is for charge flip SCALEFACTOR
    if(k_isdata || !k_isdata){

      bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      if( pass_trig ){

        std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_HN_VETO", false);
        std::vector<snu::KElectron> electronVetoColl_raw=GetElectrons(true,false,"ELECTRON_HN_VETO");
        std::vector<snu::KElectron> electronTightColl=GetElectrons(true, false, el_id.at(it_id));
        std::vector<snu::KElectron> electrons; electrons.clear();
        std::vector<snu::KElectron> electronVetoColl; electronVetoColl.clear();

        for(int i=0; i<electronVetoColl_raw.size();i++){
          if(fabs(electronVetoColl_raw.at(i).SCEta()) < 1.4442 || fabs(electronVetoColl_raw.at(i).SCEta()) > 1.5560){
            electronVetoColl.push_back(electronVetoColl_raw.at(i));
          }
        }

        for(int i=0; i<electronTightColl.size();i++){
          if(fabs(electronTightColl.at(i).SCEta()) < 1.4442 || fabs(electronTightColl.at(i).SCEta()) > 1.5560){
            electrons.push_back(electronTightColl.at(i));
          }
        }

        bool pass_lepton_number = ((muonVetoColl.size() == 0 && electronVetoColl.size() == 2 && electrons.size() == 2));
        if(!k_isdata && pass_lepton_number){

          double this_weight = weight;

          if(!k_isdata){
            this_weight *= mcdata_correction->ElectronScaleFactor(el_id.at(it_id), electrons);
            this_weight *= mcdata_correction->ElectronRecoScaleFactor(electrons);
            double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muonVetoColl, "", 1, 0, 0);
            double trigger_eff_MC   = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "", muonVetoColl, "", 1, 1, 0);
            double trigger_sf = trigger_eff_Data/trigger_eff_MC;
            this_weight *= trigger_sf;
          }
          snu::KElectron this_lep[2];
          this_lep[0] = electrons.at(0);
          this_lep[1] = electrons.at(1);

          snu::KElectron this_shifted_lep[2];
          this_shifted_lep[0] = ShiftEnergy(this_lep[0], 1./opt_shiftrate[it_id]);
          this_shifted_lep[1] = ShiftEnergy(this_lep[1], 1./opt_shiftrate[it_id]);

          bool is_SS = (this_lep[0].Charge() == this_lep[1].Charge());
          bool is_TYPE40 = ((this_lep[0].GetType() == 40) || (this_lep[1].GetType() == 40));

          TString s_region = "";
          if(fabs(this_lep[0].SCEta()) < 1.4442) s_region += "B";
          else s_region += "E";
          if(fabs(this_lep[1].SCEta()) < 1.4442) s_region += "B";
          else s_region += "E";
          if(s_region == "EB") s_region = "BE";

          double Z_Range = 30, Min_Pt = 25; int N_Bins = 24;
          TString s_SCALEFACTOR_syst = "_MassRange30_MinPt25_NBins24_SignalG";
          bool PassCuts_shifted = ((fabs((this_shifted_lep[0]+this_shifted_lep[1]).M() - 91.1876) < Z_Range) && (this_shifted_lep[1].Pt() > Min_Pt));

          if(is_SS && PassCuts_shifted){
            FillHist(el_id.at(it_id)+"_TYPE40CORR_DEN_Global_Zcandidate_Mass_CFobserved", (this_shifted_lep[0]+this_shifted_lep[1]).M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
            FillHist(el_id.at(it_id)+"_TYPE40CORR_DEN_Global_N_Events_CFobserved"+s_SCALEFACTOR_syst, 0., this_weight, 0., 1., 1);
            FillHist(el_id.at(it_id)+"_TYPE40CORR_DEN_"+s_region+"_Zcandidate_Mass_CFobserved"+s_SCALEFACTOR_syst, (this_shifted_lep[0]+this_shifted_lep[1]).M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
            FillHist(el_id.at(it_id)+"_TYPE40CORR_DEN_"+s_region+"_N_Events_CFobserved"+s_SCALEFACTOR_syst, 0., this_weight, 0., 1., 1);
            if(is_TYPE40){
              FillHist(el_id.at(it_id)+"_TYPE40CORR_NUM_Global_Zcandidate_Mass_CFobserved", (this_shifted_lep[0]+this_shifted_lep[1]).M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
              FillHist(el_id.at(it_id)+"_TYPE40CORR_NUM_Global_N_Events_CFobserved"+s_SCALEFACTOR_syst, 0., this_weight, 0., 1., 1);
              FillHist(el_id.at(it_id)+"_TYPE40CORR_NUM_"+s_region+"_Zcandidate_Mass_CFobserved"+s_SCALEFACTOR_syst, (this_shifted_lep[0]+this_shifted_lep[1]).M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
              FillHist(el_id.at(it_id)+"_TYPE40CORR_NUM_"+s_region+"_N_Events_CFobserved"+s_SCALEFACTOR_syst, 0., this_weight, 0., 1., 1);

            }
          }

        }
        if(k_isdata && pass_lepton_number){
          double this_weight = GetCFweight(electrons, false, el_id.at(it_id), false);
          double this_weight_sf = GetCFweight(electrons, true, el_id.at(it_id), false);
          double this_weight_data = 1.;

          snu::KElectron this_lep[2];
          this_lep[0] = electrons.at(0);
          this_lep[1] = electrons.at(1);
          snu::KElectron this_shifted_lep[2];
          this_shifted_lep[0] = ShiftEnergy(this_lep[0], 1./opt_shiftrate[it_id]);
          this_shifted_lep[1] = ShiftEnergy(this_lep[1], 1./opt_shiftrate[it_id]);

          bool is_SS = (this_lep[0].Charge() == this_lep[1].Charge());

          TString s_region = "";
          if(fabs(this_lep[0].SCEta()) < 1.4442) s_region += "B";
          else s_region += "E";
          if(fabs(this_lep[1].SCEta()) < 1.4442) s_region += "B";
          else s_region += "E";
          if(s_region == "EB") s_region = "BE";

          for(int it_SF_syst = 0; it_SF_syst<8; it_SF_syst++){
            double Z_Range = 30, Min_Pt = 25; int N_Bins = 24;
            TString s_SCALEFACTOR_syst = "_MassRange30_MinPt25_NBins24_SignalG";
            if(it_SF_syst == 0){}
            else if(it_SF_syst == 1){ Z_Range = 25; s_SCALEFACTOR_syst = "_MassRange25_MinPt25_NBins24_SignalG"; }
            else if(it_SF_syst == 2){ Z_Range = 35; s_SCALEFACTOR_syst = "_MassRange35_MinPt25_NBins24_SignalG"; }
            else if(it_SF_syst == 3){ Min_Pt = 22;  s_SCALEFACTOR_syst = "_MassRange30_MinPt22_NBins24_SignalG"; }
            else if(it_SF_syst == 4){ Min_Pt = 28;  s_SCALEFACTOR_syst = "_MassRange30_MinPt28_NBins24_SignalG"; }
            else if(it_SF_syst == 5){ N_Bins = 20;  s_SCALEFACTOR_syst = "_MassRange30_MinPt25_NBins20_SignalG"; }
            else if(it_SF_syst == 6){ N_Bins = 30;  s_SCALEFACTOR_syst = "_MassRange30_MinPt25_NBins30_SignalG"; }
            else if(it_SF_syst == 7){               s_SCALEFACTOR_syst = "_MassRange30_MinPt25_NBins24_SignalGG"; }
            else s_SCALEFACTOR_syst = "_WARNING";

            bool PassCuts = ((fabs((this_lep[0]+this_lep[1]).M() - 91.1876) < Z_Range) && (this_lep[1].Pt() > Min_Pt));
            bool PassCuts_shifted = ((fabs((this_shifted_lep[0]+this_shifted_lep[1]).M() - 91.1876) < Z_Range) && (this_shifted_lep[1].Pt() > Min_Pt));

            if(!is_SS && PassCuts){
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_Global_Zcandidate_Mass_CFpredicted"+s_SCALEFACTOR_syst, (this_lep[0]+this_lep[1]).M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_Global_N_Events_CFpredicted"+s_SCALEFACTOR_syst, 0., this_weight, 0., 1., 1);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_Zcandidate_Mass_CFpredicted"+s_SCALEFACTOR_syst, (this_lep[0]+this_lep[1]).M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_N_Events_CFpredicted"+s_SCALEFACTOR_syst, 0., this_weight, 0., 1., 1);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_LeadingLepton_invPt_CFpredicted"+s_SCALEFACTOR_syst, 1./this_lep[0].Pt(), this_weight, 0., 0.04, 40);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_SubLeadingLepton_invPt_CFpredicted"+s_SCALEFACTOR_syst, 1./this_lep[1].Pt(), this_weight, 0., 0.04, 40);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_Lepton_invPt_CFpredicted"+s_SCALEFACTOR_syst, 1./this_lep[0].Pt(), this_weight, 0., 0.04, 40);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_Lepton_invPt_CFpredicted"+s_SCALEFACTOR_syst, 1./this_lep[1].Pt(), this_weight, 0., 0.04, 40);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_Global_Zcandidate_Mass_CFpredictedSF"+s_SCALEFACTOR_syst, (this_lep[0]+this_lep[1]).M(), this_weight_sf, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_Global_N_Events_CFpredictedSF"+s_SCALEFACTOR_syst, 0., this_weight_sf, 0., 1., 1);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_Zcandidate_Mass_CFpredictedSF"+s_SCALEFACTOR_syst, (this_lep[0]+this_lep[1]).M(), this_weight_sf, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_N_Events_CFpredictedSF"+s_SCALEFACTOR_syst, 0., this_weight_sf, 0., 1., 1);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_LeadingLepton_invPt_CFpredictedSF"+s_SCALEFACTOR_syst, 1./this_lep[0].Pt(), this_weight_sf, 0., 0.04, 40);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_SubLeadingLepton_invPt_CFpredictedSF"+s_SCALEFACTOR_syst, 1./this_lep[1].Pt(), this_weight_sf, 0., 0.04, 40);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_Lepton_invPt_CFpredictedSF"+s_SCALEFACTOR_syst, 1./this_lep[0].Pt(), this_weight_sf, 0., 0.04, 40);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_Lepton_invPt_CFpredictedSF"+s_SCALEFACTOR_syst, 1./this_lep[1].Pt(), this_weight_sf, 0., 0.04, 40);
            }
            if(is_SS && PassCuts_shifted){
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_Global_Zcandidate_Mass_CFobserved"+s_SCALEFACTOR_syst, (this_shifted_lep[0]+this_shifted_lep[1]).M(), this_weight_data, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_Global_N_Events_CFobserved"+s_SCALEFACTOR_syst, 0., this_weight_data, 0., 1., 1);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_Zcandidate_Mass_CFobserved"+s_SCALEFACTOR_syst, (this_shifted_lep[0]+this_shifted_lep[1]).M(), this_weight_data, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_N_Events_CFobserved"+s_SCALEFACTOR_syst, 0., this_weight_data, 0., 1., 1);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_LeadingLepton_invPt_CFobserved"+s_SCALEFACTOR_syst, 1./this_lep[0].Pt(), this_weight_data, 0., 0.04, 40);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_SubLeadingLepton_invPt_CFobserved"+s_SCALEFACTOR_syst, 1./this_lep[1].Pt(), this_weight_data, 0., 0.04, 40);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_Lepton_invPt_CFobserved"+s_SCALEFACTOR_syst, 1./this_lep[0].Pt(), this_weight_data, 0., 0.04, 40);
              FillHist(el_id.at(it_id)+"_SCALEFACTOR_"+s_region+"_Lepton_invPt_CFobserved"+s_SCALEFACTOR_syst, 1./this_lep[1].Pt(), this_weight_data, 0., 0.04, 40);
            }
          }//Iterate Fit Systematics
        }//Pass lepton numb
      }//Pass dielectron trigger
    }//Is data
  }

  return;
}// End of execute event loop
  


void CFRateCalculator_Final::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void CFRateCalculator_Final::BeginCycle() throw( LQError ){
  
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

CFRateCalculator_Final::~CFRateCalculator_Final() {
  
  Message("In CFRateCalculator_Final Destructor" , INFO);
  
}


void CFRateCalculator_Final::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void CFRateCalculator_Final::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this CFRateCalculator_FinalCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void CFRateCalculator_Final::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


void CFRateCalculator_Final::GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index ){

  for( int i = it+1 ; i < truthColl.size(); i ++ ){
    if( truthColl.at(i).IndexMother() == it && truthColl.at(i).PdgId() == truthColl.at(it).PdgId() ){
      index.push_back(i);
    }
  }
  return;

}

float CFRateCalculator_Final::GetCFweight(std::vector<snu::KElectron> electrons, bool apply_sf, TString el_ID, bool do_halftest){

  if(electrons.size() > 2) return 0.;

  std::vector<snu::KElectron> lep;
  for(int i=0; i<electrons.size(); i++){
    lep.push_back(electrons.at(i));
  }

  if(lep.size()==2){
    if(lep.at(0).Charge() == lep.at(1).Charge()) return 0.;
  }

  std::vector<double> CFrate, CFweight, sf;
  for(int i=0; i<lep.size(); i++){
    CFrate.push_back(GetCFRates(lep.at(i).Pt(), lep.at(i).SCEta(), el_ID, do_halftest));
    CFweight.push_back( (CFrate.at(i)/(1-CFrate.at(i))) );
  }

  for(int i=0; i<lep.size(); i++){
    if(apply_sf){
      if(el_ID == "ELECTRON_HN_TIGHTv4"){

        if(fabs(lep.at(i).SCEta()) < 1.4442){
          sf.push_back(0.7052);
        }
        else{
          sf.push_back(0.8548);
        }
      }
      if(el_ID == "ELECTRON_HN_FAKELOOSEv7_5"){

        if(fabs(lep.at(i).SCEta()) < 1.4442){
          sf.push_back(0.7233);
        }
        else{
          sf.push_back(0.9120);
        }
      }

    }
    else sf.push_back(1.);
  }

  double cfweight = 0.;
  for(int i=0; i<lep.size(); i++){
    cfweight += sf.at(i) * (CFweight.at(i));
  }
  return cfweight;
}

float CFRateCalculator_Final::GetCFRates(double el_pt, double el_eta, TString el_ID, bool do_halftest){

  el_eta = fabs(el_eta);
  if(el_eta > 1.4442 && el_eta < 1.556) return 0.;

  double invPt = 1./el_pt;
  double a = 999., b= 999.;

  if(!do_halftest){
    if(el_ID == "ELECTRON_HN_FAKELOOSEv7_5"){
      if(el_eta < 0.8){
        if(invPt< 0.022){a=(-2.449e-05); b=(7.006e-05);}
        else{a=(0.004358); b=(-1.979e-05);}
      }
      else if(el_eta < 1.4442){
        if(invPt< 0.020){a=(-0.02186); b=(0.0009426);}
        else{a=(-0.0006124); b=(0.0005317);}
      }
      else{
        if(invPt< 0.011){a=(-0.4791); b=(0.009852);}
        else if(invPt< 0.020){a=(-0.1737); b=(0.006768);}
        else{a=(-0.01247); b=(0.003591);}
      }
    }
    if(el_ID == "ELECTRON_HN_TIGHTv4"){
      if(el_eta < 0.8){
        if(invPt< 0.023){a=(-0.001196); b=(3.806e-05);}
        else{a=(0.0008751); b=(-9.844e-06);}
      }
      else if(el_eta < 1.4442){
        if(invPt< 0.015){a=(-0.03537); b=(0.0007231);}
        else if(invPt< 0.023){a=(-0.01381); b=(0.0004019);}
        else{a=(-0.0007848); b=(0.0000972);}
      }
      else{
        if(invPt< 0.012){a=(-0.4109); b=(0.006406);}
        else if(invPt< 0.021){a=(-0.1107); b=(0.002916);}
        else{a=(-0.01974); b=(0.00103);}
      }
    }
  }
  else{
    if(el_ID == "ELECTRON_HN_TIGHTv4"){
      if(el_eta < 0.8){
        if(invPt< 0.023){a=(-0.00177); b=(4.955e-05);}
        else{a=(0.001059); b=(-1.705e-05);}
      }
      else if(el_eta < 1.4442){
        if(invPt< 0.015){a=(-0.02528); b=(0.0005667);}
        else if(invPt< 0.023){a=(-0.0154); b=(0.0004295);}
        else{a=(-0.0006096); b=(0.00008689);}
      }
      else{
        if(invPt< 0.013){a=(-0.2209); b=(0.004499);}
        else if(invPt< 0.021){a=(-0.1087); b=(0.002858);}
        else{a=(-0.02013); b=(0.001026);}
      }
    }
    else return 0.;
  }

  double rate = (a)*invPt + (b);
  if(rate < 0) rate = 0.;
  return rate;

}

snu::KElectron CFRateCalculator_Final::ShiftEnergy( snu::KElectron old_lep, double shift_rate ){

  double mass = 0.511e-3;
  snu::KElectron new_lep;
  new_lep.SetPtEtaPhiM((shift_rate*old_lep.Pt()), old_lep.Eta(), old_lep.Phi(), mass) ;
  return new_lep;

}
