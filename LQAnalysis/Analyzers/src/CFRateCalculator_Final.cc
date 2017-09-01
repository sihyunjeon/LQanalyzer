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

  double opt_shiftrate = (1.-15./1000.);

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

  //This is for CFRATE and CLOSURE and HALFTEST
  if(!k_isdata){
    if(!k_sample_name.Contains("DYtoEE")) return;

    std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
    int Njets = jetTightColl.size();
    double HT = 0.;
    for(int i=0; i<Njets; i++){
      HT += jetTightColl.at(i).Pt();
    }

    double MET = eventbase->GetEvent().MET();
    double LT = 0.;

    std::vector<snu::KElectron> electronTightColl, electronPromptColl;
    electronTightColl = GetElectrons(true, false, "ELECTRON_HN_TIGHTv4");
    electronPromptColl.clear();

    for(int i=0; i<electronTightColl.size(); i++){
      snu::KElectron this_lep=electronTightColl.at(i);
      if(this_lep.Pt() < 25) continue;
      if((this_lep.SCEta() <1.556) && (this_lep.SCEta() >1.4442)) continue;
      if(this_lep.MCIsPrompt()){
        electronPromptColl.push_back(this_lep);
        LT += this_lep.Pt();
      }
    }
    FillHist("[CHECK]size_of_electron_promptcoll", electronPromptColl.size(), weight, 0., 5., 5);

    bool one_is_CF = false;
    for(int i=0; i<electronPromptColl.size(); i++){
      int is_region = 0;
      bool is_CF = false;
      bool is_NOTCONV = false;

      snu::KElectron this_lep=electronPromptColl.at(i);

      if( (fabs(this_lep.SCEta()) < 0.9) )                                          is_region = 1;
      else if( (fabs(this_lep.SCEta()) < 1.4442) )                                  is_region = 2;
      else if( (fabs(this_lep.SCEta()) > 1.556) && (fabs(this_lep.SCEta()) < 2.5) ) is_region = 3;
      else continue;

      TString s_region = "";
      if( is_region == 1 ) s_region = "Region1";
      else if( is_region == 2 ) s_region = "Region2";
      else if( is_region == 3 ) s_region = "Region3";
      FillHist("[CHECK]region_of_electron_eta", is_region, weight, 0., 4., 4);

      if( (MCIsCF(this_lep)) ){ is_CF = true; one_is_CF = true; }
      if( !(this_lep.MCIsFromConversion()) ) is_NOTCONV = true;

      //Draw historgrams for each regions, halfsamples
      FillHist("CFRATE_Global_N_Events", 0., weight, 0., 1., 1);
      FillHist("CFRATE_Global_Pt", this_lep.Pt(), weight, 0., 500., 500);
      FillHist("CFRATE_Global_invPt", 1./this_lep.Pt(), weight, 0., 0.04, 40);
      FillHist("CFRATE_Global_SCEta", this_lep.SCEta(), weight, -3., 3., 60);
      FillHist("CFRATE_"+s_region+"_N_Events", 0., weight, 0., 1., 1);
      FillHist("CFRATE_"+s_region+"_Pt", this_lep.Pt(), weight, 0., 500., 500);
      FillHist("CFRATE_"+s_region+"_invPt", 1./this_lep.Pt(), weight, 0., 0.04, 40);
      FillHist("CFRATE_"+s_region+"_SCEta", this_lep.SCEta(), weight, -3., 3., 60);
      FillHist("CFRATE_Global_N_Events"+h_sample_suffix, 0., weight, 0., 1., 1);
      FillHist("CFRATE_Global_Pt"+h_sample_suffix, this_lep.Pt(), weight, 0., 500., 500);
      FillHist("CFRATE_Global_invPt"+h_sample_suffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
      FillHist("CFRATE_Global_SCEta"+h_sample_suffix, this_lep.SCEta(), weight, -3., 3., 60);
      FillHist("CFRATE_"+s_region+"_N_Events"+h_sample_suffix, 0., weight, 0., 1., 1);
      FillHist("CFRATE_"+s_region+"_Pt"+h_sample_suffix, this_lep.Pt(), weight, 0., 500., 500);
      FillHist("CFRATE_"+s_region+"_invPt"+h_sample_suffix, 1./this_lep.Pt(), weight, 0., 0.04, 40);
      FillHist("CFRATE_"+s_region+"_SCEta"+h_sample_suffix, this_lep.SCEta(), weight, -3., 3., 60);

      //Draw histograms only for half sampletest
      FillHist("HALFTEST_Global_MET"+h_sample_suffix, MET, 1., 0., 80., 8);
      FillHist("HALFTEST_Global_METsqdivST"+h_sample_suffix, MET*MET/(MET+LT+HT), 1., 0., 40., 8);
      FillHist("HALFTEST_Global_N_Jets"+h_sample_suffix, Njets, 1., 0., 8., 8);
      FillHist("HALFTEST_Global_HT"+h_sample_suffix, HT, 1., 0., 160., 8);

      double halftestrate = GetCFRates((this_lep.Pt()), (this_lep.SCEta()), "ELECTRON_HN_TIGHTv4", true);
      double halftestweight = halftestrate/(1.-halftestrate);
      FillHist("HALFTEST_Global_MET"+h_sample_suffix+"_CFpredicted", MET, halftestweight, 0., 80., 8);
      FillHist("HALFTEST_Global_METsqdivST"+h_sample_suffix+"_CFpredicted", MET*MET/(MET+LT+HT), halftestweight, 0., 40., 8);
      FillHist("HALFTEST_Global_N_Jets"+h_sample_suffix+"_CFpredicted", Njets, halftestweight, 0., 8., 8);
      FillHist("HALFTEST_Global_HT"+h_sample_suffix+"_CFpredicted", HT, halftestweight, 0., 160., 8);

      if(is_CF){
        FillHist("CFRATE_Global_N_Events_CF", 0., weight, 0., 1., 1);
        FillHist("CFRATE_Global_Pt_CF", this_lep.Pt(), weight, 0., 500., 500);
        FillHist("CFRATE_Global_invPt_CF", 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist("CFRATE_Global_SCEta_CF", this_lep.SCEta(), weight, -3., 3., 60);
        FillHist("CFRATE_"+s_region+"_N_Events_CF", 0., weight, 0., 1., 1);
        FillHist("CFRATE_"+s_region+"_Pt_CF", this_lep.Pt(), weight, 0., 500., 500);
        FillHist("CFRATE_"+s_region+"_invPt_CF", 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist("CFRATE_"+s_region+"_SCEta_CF", this_lep.SCEta(), weight, -3., 3., 60);
        FillHist("CFRATE_Global_N_Events"+h_sample_suffix+"_CF", 0., weight, 0., 1., 1);
        FillHist("CFRATE_Global_Pt"+h_sample_suffix+"_CF", this_lep.Pt(), weight, 0., 500., 500);
        FillHist("CFRATE_Global_invPt"+h_sample_suffix+"_CF", 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist("CFRATE_Global_SCEta"+h_sample_suffix+"_CF", this_lep.SCEta(), weight, -3., 3., 60);
        FillHist("CFRATE_"+s_region+"_N_Events"+h_sample_suffix+"_CF", 0., weight, 0., 1., 1);
        FillHist("CFRATE_"+s_region+"_Pt"+h_sample_suffix+"_CF", this_lep.Pt(), weight, 0., 500., 500);
        FillHist("CFRATE_"+s_region+"_invPt"+h_sample_suffix+"_CF", 1./this_lep.Pt(), weight, 0., 0.04, 40);
        FillHist("CFRATE_"+s_region+"_SCEta"+h_sample_suffix+"_CF", this_lep.SCEta(), weight, -3., 3., 60);

        FillHist("HALFTEST_Global_MET"+h_sample_suffix+"_CFobserved", MET, 1., 0., 80., 8);
        FillHist("HALFTEST_Global_METsqdivST"+h_sample_suffix+"_CFobserved", MET*MET/(MET+LT+HT), 1., 0., 40., 8);
        FillHist("HALFTEST_Global_N_Jets"+h_sample_suffix+"_CFobserved", Njets, 1., 0., 8., 8);
        FillHist("HALFTEST_Global_HT"+h_sample_suffix+"_CFobserved", HT, 1., 0., 160., 8);

        if(is_NOTCONV){
          FillHist("CFRATE_Global_N_Events_CF_NOTCONV", 0., weight, 0., 1., 1);
          FillHist("CFRATE_Global_Pt_CF_NOTCONV", this_lep.Pt(), weight, 0., 500., 500);
          FillHist("CFRATE_Global_invPt_CF_NOTCONV", 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist("CFRATE_Global_SCEta_CF_NOTCONV", this_lep.SCEta(), weight, -3., 3., 60);
          FillHist("CFRATE_"+s_region+"_N_Events_CF_NOTCONV", 0., weight, 0., 1., 1);
          FillHist("CFRATE_"+s_region+"_Pt_CF_NOTCONV", this_lep.Pt(), weight, 0., 500., 500);
          FillHist("CFRATE_"+s_region+"_invPt_CF_NOTCONV", 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist("CFRATE_"+s_region+"_SCEta_CF_NOTCONV", this_lep.SCEta(), weight, -3., 3., 60);
          FillHist("CFRATE_Global_N_Events"+h_sample_suffix+"_CF_NOTCONV", 0., weight, 0., 1., 1);
          FillHist("CFRATE_Global_Pt"+h_sample_suffix+"_CF_NOTCONV", this_lep.Pt(), weight, 0., 500., 500);
          FillHist("CFRATE_Global_invPt"+h_sample_suffix+"_CF_NOTCONV", 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist("CFRATE_Global_SCEta"+h_sample_suffix+"_CF_NOTCONV", this_lep.SCEta(), weight, -3., 3., 60);
          FillHist("CFRATE_"+s_region+"_N_Events"+h_sample_suffix+"_CF_NOTCONV", 0., weight, 0., 1., 1);
          FillHist("CFRATE_"+s_region+"_Pt"+h_sample_suffix+"_CF_NOTCONV", this_lep.Pt(), weight, 0., 500., 500);
          FillHist("CFRATE_"+s_region+"_invPt"+h_sample_suffix+"_CF_NOTCONV", 1./this_lep.Pt(), weight, 0., 0.04, 40);
          FillHist("CFRATE_"+s_region+"_SCEta"+h_sample_suffix+"_CF_NOTCONV", this_lep.SCEta(), weight, -3., 3., 60);

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

      if((this_lep[0].Charge() == this_lep[1].Charge()) && !one_is_CF) FillHist("[CHECK]Samesign_electrons_but_not_chargeflip", 0., 1., 0., 1., 1);

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
              FillHist("[CHECK]MCmatched_leptons_not_2", truth_MCmatched_index.size(), weight, 0., 4., 4);
            }
            else{
              if((truth_MCmatched_index.at(0)+truth_MCmatched_index.at(1) != 1) || (reco_MCmatched_index.at(0)+reco_MCmatched_index.at(1) != 1)){
                FillHist("[CHECK]MCmatched_leptons_doubly_found", 0., weight, 0., 1., 1);
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
                  if(cf_found) FillHist("[CHECK]Chargeflipped_electron_doubly_found", 0., 1., 0., 1., 1);
                  cf_found = true;
                }
                if(!cf_found) FillHist("[CHECK]Chargeflipped_electron_not_found", 0., 1., 0., 1., 1);
                else{
                  FillHist("SHIFTRATE_Electron_Pt_RecoLevel", this_lep[reco_cf_index].Pt(), 1., 0., 500., 500);
                  FillHist("SHIFTRATE_Electron_Pt_TruthLevel", truth_lep[truth_cf_index].Pt(), 1., 0., 500., 500);
                  FillHist("SHIFTRATE_Electron_Pt_Ratio", (this_lep[reco_cf_index].Pt()/truth_lep[truth_cf_index].Pt()), 1., 0., 2., 200);
                  if(!MCIsCF(this_lep[reco_cf_index])) FillHist("[CHECK]Chargeflipped_electron_not_found_GetType", this_lep[reco_cf_index].GetType(), 1., 0., 50., 50);
                }
                FillHist("SHIFTRATE_Zcandidate_Mass_CFobserved", (this_lep[0]+this_lep[1]).M(), 1., (91.1876-40.), (91.1876+40.), 32);
                FillHist("CLOSURE_Zcandidate_Mass_CFobserved", (this_lep[0]+this_lep[1]).M(), 1., 70., 110., 16);
                FillHist("CLOSURE_LeadingLepton_Pt_CFobserved", this_lep[0].Pt(), 1., 0., 100., 10);
                FillHist("CLOSURE_SubLeadingLepton_Pt_CFobserved", this_lep[1].Pt(), 1., 0., 100., 10);
                FillHist("CLOSURE_N_Events_CFobserved", 0., 1., 0., 1., 1);

              }//Pass tight Z requirements for SHIFTRATE and CLOSURE
            }//Is SS dielectron
            if(!is_SS){
              double this_weight=GetCFweight(electronPromptColl, false, "ELECTRON_HN_TIGHTv4", false);

              bool is_shifted_Z_tight = false;
              for(int it_shift=0; it_shift<51; it_shift++){
                snu::KElectron shifted_lep[2];
                TString s_shift = "_Shifted_"+TString::Itoa(it_shift, 10)+"div1000";

                double shift_rate = 1.-0.001*it_shift;
                shifted_lep[0] = ShiftEnergy(this_lep[0], shift_rate);
                shifted_lep[1] = ShiftEnergy(this_lep[1], shift_rate);
                is_shifted_Z_tight = (fabs((shifted_lep[0] + shifted_lep[1]).M() - 91.1876) < 15.);
                if(is_shifted_Z_tight){
                  FillHist("SHIFTRATE_Zcandidate_Mass_CFpredicted"+s_shift, (shifted_lep[0]+shifted_lep[1]).M(), this_weight, (91.1876-40.), (91.1876+40.), 32);
                }//Pass tight Z requirements after shift for SHIFTRATE
              }//Iterate shiftrates for optimal value

              snu::KElectron this_shifted_lep[2];
              this_shifted_lep[0] = ShiftEnergy(this_lep[0], opt_shiftrate);
              this_shifted_lep[1] = ShiftEnergy(this_lep[1], opt_shiftrate);

              is_shifted_Z_tight = ((fabs((this_shifted_lep[0] + this_shifted_lep[1]).M() - 91.1876) < 15.) && (this_shifted_lep[1].Pt() > 25));
              if(is_shifted_Z_tight){
                FillHist("CLOSURE_Zcandidate_Mass_CFpredicted", (this_shifted_lep[0]+this_shifted_lep[1]).M(), this_weight, 70., 110., 16);
                FillHist("CLOSURE_LeadingLepton_Pt_CFpredicted", this_shifted_lep[0].Pt(), this_weight, 0., 100., 10);
                FillHist("CLOSURE_SubLeadingLepton_Pt_CFpredicted", this_shifted_lep[1].Pt(), this_weight, 0., 100., 10);
                FillHist("CLOSURE_N_Events_CFpredicted", 0., this_weight, 0., 1., 1);

              }//Pass tight Z requirements after shift for CLOSURE
            }//Is OS dielectron
          }//Well matched
        }//Find two truth matched electrons
      }//Pass loose Z requirements for SHIFTRATE and CLOSURE without any charge requirements
    }//Has two prompt electrons

  }//Is not data


  //This is for charge flip SCALEFACTOR
  else{

    bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
    if( pass_trig ){

      std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_HN_VETO", false);
      std::vector<snu::KElectron> electronVetoColl = GetElectrons(false,false,"ELECTRON_HN_VETO");
      std::vector<snu::KElectron> electronTightColl=GetElectrons(true, false, "ELECTRON_HN_TIGHTv4");
      std::vector<snu::KElectron> electrons=GetElectrons(true, false, "ELECTRON_HN_FAKELOOSEv2");

      bool pass_lepton_number = ((muonVetoColl.size() == 0 && electronVetoColl.size() == 2 && electrons.size() == 2 && electronTightColl.size() == 2));

      if( pass_lepton_number ){
        double this_weight = GetCFweight(electrons, false, "ELECTRON_HN_TIGHTv4", false);
        double this_weight_sf = GetCFweight(electrons, true, "ELECTRON_HN_TIGHTv4", false);

        snu::KElectron this_lep[2];
        this_lep[0] = electrons.at(0);
        this_lep[1] = electrons.at(1);

        snu::KElectron this_shifted_lep[2];
        this_shifted_lep[0] = ShiftEnergy(this_lep[0], 1./opt_shiftrate);
        this_shifted_lep[1] = ShiftEnergy(this_lep[1], 1./opt_shiftrate);

        bool is_SS = (this_lep[0].Charge() == this_lep[1].Charge());

        TString s_region = "";
        if(fabs(this_lep[0].SCEta()) < 1.4442) s_region += "B";
        else if( (fabs(this_lep[0].SCEta()) > 1.556) && (fabs(this_lep[0].SCEta()) < 2.5) ) s_region += "E";
        if(fabs(this_lep[1].SCEta()) < 1.4442) s_region += "B";
        else if( (fabs(this_lep[1].SCEta()) > 1.556) && (fabs(this_lep[1].SCEta()) < 2.5) ) s_region += "E";
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
            FillHist("SCALEFACTOR_Global_Zcandidate_Mass_CFpredicted"+s_SCALEFACTOR_syst, (this_lep[0]+this_lep[1]).M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
            FillHist("SCALEFACTOR_Global_N_Events_CFpredicted"+s_SCALEFACTOR_syst, 0., this_weight, 0., 1., 1);
            FillHist("SCALEFACTOR_"+s_region+"_Zcandidate_Mass_CFpredicted"+s_SCALEFACTOR_syst, (this_lep[0]+this_lep[1]).M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
            FillHist("SCALEFACTOR_"+s_region+"_N_Events_CFpredicted"+s_SCALEFACTOR_syst, 0., this_weight, 0., 1., 1);

            FillHist("SCALEFACTOR_Global_Zcandidate_Mass_CFpredictedSF"+s_SCALEFACTOR_syst, (this_lep[0]+this_lep[1]).M(), this_weight_sf, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
            FillHist("SCALEFACTOR_Global_N_Events_CFpredictedSF"+s_SCALEFACTOR_syst, 0., this_weight_sf, 0., 1., 1);
            FillHist("SCALEFACTOR_"+s_region+"_Zcandidate_Mass_CFpredictedSF"+s_SCALEFACTOR_syst, (this_lep[0]+this_lep[1]).M(), this_weight_sf, 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
            FillHist("SCALEFACTOR_"+s_region+"_N_Events_CFpredictedSF"+s_SCALEFACTOR_syst, 0., this_weight_sf, 0., 1., 1);

          }
          if(is_SS && PassCuts_shifted){
            FillHist("SCALEFACTOR_Global_Zcandidate_Mass_CFobserved"+s_SCALEFACTOR_syst, (this_shifted_lep[0]+this_shifted_lep[1]).M(), 1., 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
            FillHist("SCALEFACTOR_Global_N_Events_CFobserved"+s_SCALEFACTOR_syst, 0., 1., 0., 1., 1);
            FillHist("SCALEFACTOR_"+s_region+"_Zcandidate_Mass_CFobserved"+s_SCALEFACTOR_syst, (this_shifted_lep[0]+this_shifted_lep[1]).M(), 1., 91.1876-Z_Range, 91.1876+Z_Range, N_Bins);
            FillHist("SCALEFACTOR_"+s_region+"_N_Events_CFobserved"+s_SCALEFACTOR_syst, 0., 1., 0., 1., 1);
          }
        }//Iterate Fit Systematics
      }//Pass lepton numb
    }//Pass dielectron trigger
  }//Is data


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
      if(fabs(lep.at(i).SCEta()) < 1.4442){
        sf.push_back(0.691722);
      }
      else{
        sf.push_back(0.68301);
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
    if(el_eta < 0.9){
      if(invPt< 0.023){a=(-0.00138148); b=(4.33442e-05);}
      else{a=(0.00101034); b=(-1.14551e-05);}
    }
    else if(el_eta < 1.4442){
      if(invPt< 0.015){a=(-0.042964); b=(0.000866971);}
      else if(invPt< 0.023){a=(-0.0152852); b=(0.000452217);}
      else{a=(-0.00154575); b=(0.000127211);}
    }
    else{
      if(invPt< 0.012){a=(-0.423831); b=(0.00636555);}
      else if(invPt< 0.020){a=(-0.103982); b=(0.00254955);}
      else{a=(-0.0160296); b=(0.000767227);}
    }
  }
  else{
    if(el_eta < 0.9){
      if(invPt< 0.023){a=(-0.00313374); b=(8.01418e-05);}
      else{a=(0.000961691); b=(-1.42211e-05);}
    }
    else if(el_eta < 1.4442){
      if(invPt< 0.015){a=(-0.0282759); b=(0.000654615);}
      else if(invPt< 0.023){a=(-0.0156543); b=(0.000452466);}
      else{a=(-0.000931418); b=(0.000105082);}
    }
    else{
      if(invPt< 0.012){a=(-0.218108); b=(0.00424753);}
      else if(invPt< 0.020){a=(-0.097057); b=(0.00239728);}
      else{a=(-0.0168018); b=(0.000775736);}
    }
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
