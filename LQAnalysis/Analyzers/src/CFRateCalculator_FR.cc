// $Id: CFRateCalculator_FR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQCFRateCalculator_FR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "CFRateCalculator_FR.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (CFRateCalculator_FR);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
CFRateCalculator_FR::CFRateCalculator_FR() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("CFRateCalculator_FR");
  
  Message("In CFRateCalculator_FR constructor", INFO);
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


void CFRateCalculator_FR::InitialiseAnalysis() throw( LQError ) {
  
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


void CFRateCalculator_FR::ExecuteEvents()throw( LQError ){

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

  return;
}// End of execute event loop
  


void CFRateCalculator_FR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void CFRateCalculator_FR::BeginCycle() throw( LQError ){
  
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

CFRateCalculator_FR::~CFRateCalculator_FR() {
  
  Message("In CFRateCalculator_FR Destructor" , INFO);
  
}


void CFRateCalculator_FR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void CFRateCalculator_FR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this CFRateCalculator_FRCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void CFRateCalculator_FR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


void CFRateCalculator_FR::GENFindDecayIndex( std::vector<snu::KTruth> truthColl,  int it, std::vector<int>& index ){

  for( int i = it+1 ; i < truthColl.size(); i ++ ){
    if( truthColl.at(i).IndexMother() == it && truthColl.at(i).PdgId() == truthColl.at(it).PdgId() ){
      index.push_back(i);
    }
  }
  return;

}



void CFRateCalculator_FR::CFvalidation(void){

  bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if( !pass_trig ) return;

  double Z_mass = 91.1876;

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  if(METPt > 30) return;

  TString CFsample="";

  for(int aaa=0; aaa<1; aaa++){
    TString el_ID="", el_looseID="";
    TString IDsuffix="";
    TString method_fake = "";

    if(aaa==0){
      el_ID = "HN";
      method_fake = "mva";
    }
    if(aaa==1){
      el_ID = "MVA";
      method_fake = "dijet_ajet40";
    }

    IDsuffix = "_"+el_ID+"TIGHT";
    el_looseID = "ELECTRON_"+el_ID+"_FAKELOOSE";
    el_ID = "ELECTRON_"+el_ID+"_TIGHT";

    std::vector<snu::KElectron> electronLooseColl = GetElectrons(true, false, el_looseID);
    std::vector<snu::KElectron> electronTightColl = GetElectrons(true, false, el_ID);
    if(electronLooseColl.size() != 2) continue;
    if(electronTightColl.size() == 2) continue;
    std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_VETO", false);
    if( muonLooseColl.size() != 0) continue;

    FillHist("[CHECK]n_of_muons"+IDsuffix+CFsample, muonLooseColl.size(), 1., 0., 5., 5);//check no muons
    FillHist("[CHECK]n_of_electrons"+IDsuffix+CFsample, electronTightColl.size(), 1., 0., 5., 5);//check two electrons

    // define leptons and give Pt, MET cuts
    snu::KParticle lep[2];
    lep[0] = electronLooseColl.at(0);
    lep[1] = electronLooseColl.at(1);

    bool is_SS = false;
    if( (lep[0].Charge() != lep[1].Charge()) ) continue;

    if( lep[0].Pt() < 25 || lep[1].Pt() < 25 ) continue;

    FillHist("[CHECK]MET"+IDsuffix+CFsample, METPt, 1., 0., 100., 100);

    TString region = "";
    if( (fabs(lep[0].Eta()) < 0.9) )                                      region += "iB";
    else if( (fabs(lep[0].Eta()) < 1.4442) )                              region += "oB";
    else if( (fabs(lep[0].Eta()) > 1.556) && (fabs(lep[0].Eta()) < 2.5) ) region += "E";
    else continue;

    if( (fabs(lep[1].Eta()) < 0.9) )                                      region += "iB";
    else if( (fabs(lep[1].Eta()) < 1.4442) )                              region += "oB";
    else if( (fabs(lep[1].Eta()) > 1.556) && (fabs(lep[1].Eta()) < 2.5) ) region += "E";
    else continue;

    if((region == "iBE") || (region == "EiB") || (region == "oBE") || (region == "EoB"))  region = "BE";
    if((region == "iBiB") || (region == "iBoB") || (region == "oBiB") || (region == "oBoB")) region = "BB";

    snu::KParticle Z_candidate;//define Z candidate after shifting energy (SS : better fitting, energy scale up // OS : photon radiation E loss)
    Z_candidate = (lep[0] + lep[1]);

    TString Zsuffix[4] = {"_narrowZ", "", "_wideZ", "_verywideZ"};
    double Zwidth[4] = {-2.5, 0., 2.5, 5.,};

    std::vector<snu::KJet> jetTightColl = GetJets("JET_HN", 30., 2.4);
    int Njets = jetTightColl.size();
    TString s_njets = "";
    if( Njets == 0 ) s_njets = "JETS0";
    if( Njets != 0 ) s_njets = "JETS";

    double fake_weight = 0.;
    fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 0, electronLooseColl, el_ID, 2, el_looseID, method_fake);

    for(int Z_it = 0; Z_it < 4; Z_it ++){

      bool Z_selection = (((Z_candidate.M() - Z_mass) < (10.+Zwidth[Z_it])) && ((Z_mass - Z_candidate.M()) < (10.+Zwidth[Z_it])));

      if(Z_it ==0 && aaa==0){
        FillHist("FIT_observed_Z_mass_global"+IDsuffix, Z_candidate.M(), fake_weight, 60., 120., 60);
        FillHist("FIT_observed_n_events_global"+IDsuffix, 0., fake_weight, 0., 1., 1);
        FillHist("FIT_observed_Z_mass_"+region+IDsuffix, Z_candidate.M(), fake_weight, 60., 120., 60);
        FillHist("FIT_observed_n_events_"+region+IDsuffix, 0., fake_weight, 0., 1., 1);
      }

      if( Z_selection ){

        FillHist("observed_Z_mass_global"+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., fake_weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+region+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_"+region+IDsuffix+CFsample+Zsuffix[Z_it], 0., fake_weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+s_njets+"_global"+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_"+s_njets+"_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., fake_weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+s_njets+"_"+region+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_"+s_njets+"_"+region+IDsuffix+CFsample+Zsuffix[Z_it], 0., fake_weight, 0., 1., 1);
        if(Njets == 0) FillHist("observed_n_jets_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., fake_weight, 0., 2., 2);
        if(Njets != 0) FillHist("observed_n_jets_global"+IDsuffix+CFsample+Zsuffix[Z_it], fake_weight, 1., 0., 2., 2);
        if(Njets == 1){
          FillHist("observed_Z_mass_JETS1_global"+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), fake_weight, 70., 110., 40);
          FillHist("observed_n_events_JETS1_global"+IDsuffix+CFsample+Zsuffix[Z_it], 0., fake_weight, 0., 1., 1);
          FillHist("observed_Z_mass_JETS1_"+region+IDsuffix+CFsample+Zsuffix[Z_it], Z_candidate.M(), fake_weight, 70., 110., 40);
          FillHist("observed_n_events_JETS1_"+region+IDsuffix+CFsample+Zsuffix[Z_it], 0., fake_weight, 0., 1., 1);
        }
      }// is_SS
    }
  }
  return;
}

snu::KParticle CFRateCalculator_FR::ShiftEnergy( snu::KParticle old_lep, double shift_rate ){

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

double CFRateCalculator_FR::GetCFRates(int sys, double el_pt, double el_eta, TString el_ID, bool apply_sf){

  el_eta = fabs(el_eta);
  if(el_eta > 1.4442 && el_eta < 1.556) return 0.;

  double sf =999.;
  if(el_eta < 1.4442) sf = 1.0;
  else sf = 1.0;

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
    if(invPt< 0.006){
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
  if(apply_sf) rate = rate*sf;
  return rate;

}

double CFRateCalculator_FR::Get2DCFRates(bool apply_sf, double el_pt, double el_eta, TString el_ID, TString halfsample, TString CFsample){

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

