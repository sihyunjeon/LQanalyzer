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



 
  bool trig_pass = (PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v") || PassTrigger("HLT_Ele17_CaloIdL_GsfTrkIdVL_v"));
  TString s_trigger = "single";
  if(!trig_pass) return;

  TString sample_suffix = "";
  if(k_sample_name.Contains("DYtoEE")) sample_suffix = "_powheg";
  else if(k_sample_name.Contains("DYJets_MG")) sample_suffix = "_madgraph";
  else return;

  std::vector<snu::KElectron> electronVLooseColl = GetElectrons(false, false, "ELECTRON_HN_FAKELOOSE");
  double this_reliso_el = 0.08;

  snu::KEvent event = eventbase->GetEvent();

  int N_sys = (2*2+1);
  for(int it_sys = 0; it_sys<N_sys; it_sys++){

    double MET = eventbase->GetEvent().MET();
    TString this_syst;

    if(it_sys==0){
      this_syst = "ElectronE_Up";//
      MET = event.PFMETShifted(snu::KEvent::ElectronEn, snu::KEvent::up);
    }
    else if(it_sys==1){
      this_syst = "ElectronE_Down";//
      MET = event.PFMETShifted(snu::KEvent::ElectronEn, snu::KEvent::down);
    }
    else if(it_sys==2){//
      this_syst = "Unclustered_Up";
      MET = event.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::up);
    }
    else if(it_sys==3){//
      this_syst = "Unclustered_Down";
      MET = event.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::down);
    }
    else if(it_sys==4){//
      this_syst = "Central";
    }
    else{
      Message("it_sys out of range!" , INFO);
      return;
    }

    std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
    int Njets = jetTightColl.size();
    double HT = 0.;
    for(int i=0; i<Njets; i++){
      HT += jetTightColl.at(i).Pt();
    }

    std::vector<snu::KElectron> electronTightColl;
    if(this_syst == "ElectronE_Up"){
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedUp(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        double new_reliso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedUp();
        if( this_electron.Pt() >= 25. && new_reliso < this_reliso_el ) electronTightColl.push_back( this_electron );
      }
    }
    else if(this_syst == "ElectronE_Down"){
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedDown(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        double new_reliso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedDown();
        if( this_electron.Pt() >= 25. && new_reliso < this_reliso_el ) electronTightColl.push_back( this_electron );
      }
    }
    else{
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        if( this_electron.Pt() >= 25. && this_electron.PFRelIso(0.3) < this_reliso_el ) electronTightColl.push_back( this_electron );
      }
    }

    if(electronTightColl.size() == 0) continue;

    for(int aa=0; aa<2; aa++){
      if(aa==0) weight = pileup_reweight;
      if(aa==1) weight = 1.;;

      for(int aaa=0; aaa<1; aaa++){

        std::vector<snu::KElectron> electronPromptColl;
        electronPromptColl.clear();

        TString IDsuffix = "";
        TString el_ID = "";

        if(aaa == 0){
          el_ID = "ELECTRON_HN_TIGHTv4";
          IDsuffix = "_HNTIGHT_"+this_syst;
        }
        else return;

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

      }//for different tight id iteration
    }
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



snu::KElectron CFRateCalculator_syst::ShiftEnergy( snu::KElectron old_lep, double shift_rate ){

  double mass = 0.511e-3;
  snu::KElectron new_lep;
  new_lep.SetPtEtaPhiM((shift_rate*old_lep.Pt()), old_lep.Eta(), old_lep.Phi(), mass) ;
  return new_lep;

}


double CFRateCalculator_syst::GetCFweight(int sys, std::vector<snu::KElectron> electrons, bool apply_sf, TString Zwidth){

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


double CFRateCalculator_syst::GetCFRates(int sys, double el_pt, double el_eta, TString el_ID){

  if(el_ID != "ELECTRON_HN_TIGHTv4") return 0.;

  el_eta = fabs(el_eta);
  if(el_eta > 1.4442 && el_eta < 1.556) return 0.;

  double invPt = 1./el_pt;
  double a = 999., b= 999.;
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

  double rate = (a)*invPt + (b);
  if(rate < 0) rate = 0.;
  return rate;
}

