// $Id: DiEl_SS.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQDiEl_SS Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "DiEl_SS.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (DiEl_SS);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
DiEl_SS::DiEl_SS() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("DiEl_SS");
  
  Message("In DiEl_SS constructor", INFO);
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
  MVATIGHT_CF_sampleB_hist_powheg = (TH2F*)file_powheg ->Get("Pt_eta_global_CF_MVATIGHT_PU_sampleB")->Clone();


}


void DiEl_SS::InitialiseAnalysis() throw( LQError ) {
  
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


void DiEl_SS::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  FillHist("signalEff", 1.5, 1., 0., 15., 15); // NoCut
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);


   if(!PassMETFilter()) return;     /// Initial event cuts : 
   FillCutFlow("EventCut", weight);
   FillHist("signalEff", 2.5, 1., 0., 15., 15); // METFilter

   /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton
   
   if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex  
   FillHist("signalEff", 3.5, 1., 0., 15., 15); // vtx cut

/*   float pileup_reweight=(1.0);
   if (!k_isdata) {   pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());}*/
     
   
   TString diel_trig = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
/*   TString dimuon_trigmuon_trig1="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
   TString dimuon_trigmuon_trig2="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v";
   TString dimuon_trigmuon_trig3="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v";
   TString dimuon_trigmuon_trig4="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v";*/
   // Now you should do an OR of 4 triggers 
 
   vector<TString> trignames;
   trignames.push_back(diel_trig);
/*   trignames.push_back(dimuon_trigmuon_trig1);
   trignames.push_back(dimuon_trigmuon_trig2);
   trignames.push_back(dimuon_trigmuon_trig3);
   trignames.push_back(dimuon_trigmuon_trig4);*/
   bool trig_pass= PassTrigger(diel_trig);
   if(!trig_pass) return;
   FillHist("signalEff", 4.5, 1., 0., 15., 15); // passing trigger

   std::vector<snu::KElectron> electrons;
   if(k_running_nonprompt){
     electrons = GetElectrons(true, false,"ELECTRON_HN_FAKELOOSE");
   }
   else electrons = GetElectrons(false,false, "ELECTRON_HN_TIGHT");

//   std::vector<snu::KElectron> electrons =  GetElectrons(false,false, "ELECTRON_HN_TIGHT");
   std::vector<snu::KElectron> electrons_veto =  GetElectrons(false,false, "ELECTRON_HN_VETO");
   /*
     
   std::vector<snu::KElectron> electrons =  GetElectrons(BaseSelection::ELECTRON_NOCUT);  ... WONT WORK
   std::vector<snu::KElectron> electrons =  GetElectrons("ELECTRON_NOCUT");               ... WILL WORK  
   
   std::vector<snu::KElectron> electrons =  GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);  ... WILL WORK  
   std::vector<snu::KElectron> electrons =  GetElectrons("ELECTRON_POG_TIGHT");                ... WILL WORK  
   
   */

   //   std::vector<snu::KElectron> electrons2 =  GetElectrons(BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0);

   std::vector<snu::KJet> jets = GetJets("JET_HN_TChannel");
//   int nbjet = NBJet(GetJets("JET_HN"));
   int nbjet = 0;
   for(unsigned int ijet=0; ijet<jets.size(); ijet++){
     if(jets.at(ijet).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) nbjet++;
   }   
   int njet = jets.size();
   std::vector<snu::KMuon> muons = GetMuons("MUON_HN_TRI_TIGHT");
   std::vector<snu::KMuon> muons_veto =GetMuons("MUON_HN_VETO"); 


//   mcdata_correction->CorrectMuonMomentum(muons,eventbase->GetTruth()); /// CorrectMuonMomentum(muons);  will also work as Funcion in AnalyzerCore just calls mcdata_correction function
   
   float ev_weight = 1.;

   float trigger_sf = 1.;
   float id_iso_sf = 1.;
   float trigger_ps = 1.;
   float reco_sf = 1.;
   float pu_reweight = 1.;  
 
   if(!isData){
     pu_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());
     trigger_ps = WeightByTrigger(diel_trig, TargetLumi);
     id_iso_sf = mcdata_correction->ElectronScaleFactor("ELECTRON_HN_TIGHT", electrons); 
     reco_sf = mcdata_correction->ElectronRecoScaleFactor(electrons); 
     ev_weight = weight * trigger_ps * id_iso_sf * pu_reweight;
   }
   if(k_running_nonprompt){
     ev_weight = 1.;
     ev_weight *= m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electrons);
   }
   ev_weight *= m_datadriven_bkg->WeightCFEvent(electrons, k_running_chargeflip);

   if(electrons.size() == 2){
     if(electrons_veto.size() > 2) return;
     if(muons_veto.size() > 0) return;
     FillHist("signalEff", 5.5, 1., 0., 15., 15); // 2 tight electrons
     if(electrons.at(0).Charge()==electrons.at(1).Charge()) return;
     FillHist("signalEff", 6.5, 1., 0., 15., 15); // same-sign electron pair
     snu::KParticle X = electrons.at(0)+electrons.at(1); 
     if(electrons.at(0).Pt() > 25. && electrons.at(1).Pt() > 20.){
       FillHist("signalEff", 7.5, 1., 0., 15., 15); // electron pT cuts
       if(X.M() < 10.) return;
       if(fabs(X.M()-10.) < 10.) return;
       FillHist("signalEff", 8.5, 1., 0., 15., 15); // electron mass cut

       if((njet==1) && (X.M()>100.)){
         FillHist("mass_ee_SS1jet", X.M(), ev_weight, 0., 1500., 300);
         FillHist("pT_ee_SS1jet", X.Pt(), ev_weight, 0., 500., 100);
         FillHist("pT_e1_SS1jet", electrons.at(0).Pt(), ev_weight, 0., 500., 100);
         FillHist("pT_e2_SS1jet", electrons.at(1).Pt(), ev_weight, 0., 500., 100);
         FillHist("MET_SS1jet", eventbase->GetEvent().PFMET(), ev_weight, 0., 500., 100);
       }

//       if(eventbase->GetEvent().PFMET() > 35.) return;
//       FillHist("signalEff", 9.5, 1., 0., 15., 15); // MET cut

       if(jets.size() > 1){
         FillHist("signalEff", 9.5, 1., 0., 15., 15); // >= 2jets (s-channel)
         if(jets.size() > 3){
           FillHist("signalEff", 10.5, 1., 0., 15., 15); // >= 4jets (t-channel)
           float wmass = 10000.; int j1 = 0; int j2 = 0;
           bool forward_jet(false), back_jet(false);

           //Select 2 jets closest to m(W)
           for(unsigned int ij=0; ij < jets.size(); ij++){
             for(unsigned int ij2=ij+1; ij2 < jets.size(); ij2++){
               bool btag1 = IsBTagged(jets.at(ij), snu::KJet::CSVv2, snu::KJet::Medium);
               bool btag2 = IsBTagged(jets.at(ij2), snu::KJet::CSVv2, snu::KJet::Medium);
               if((!btag1) && (!btag2)){
                 snu::KParticle W = jets.at(ij) + jets.at(ij2);
                 if(fabs(W.M()-80.4) < wmass){
                    wmass = fabs(W.M()-80.4); j1 = ij; j2 = ij2;
                 }
               }
             }
           }
           
           //Require at least 1 jet in each region that |eta| > 1.5
           for(unsigned int ij3=0; ij3 < jets.size(); ij3++){
             if((jets.at(ij3).Eta()>1.5)) forward_jet=true;
             if((jets.at(ij3).Eta()<-1.5)) back_jet=true;
           }

           snu::KParticle HN1 = electrons.at(0) + jets.at(j1) + jets.at(j2);
           snu::KParticle HN2 = electrons.at(1) + jets.at(j1) + jets.at(j2);

           if(forward_jet && back_jet){

	     //get cfrates and weight
             double sf_CFrate[2] = {-999.};
             snu::KParticle lep[2];
             lep[0] = electrons.at(0);
             lep[1] = electrons.at(1);
             sf_CFrate[0] = Get2DCFRates(true, lep[0].Pt(), fabs(lep[0].Eta()), "ELECTRON_HN_TIGHT", "", "_madgraph"); //sf applied cf rates
             sf_CFrate[1] = Get2DCFRates(true, lep[1].Pt(), fabs(lep[1].Eta()), "ELECTRON_HN_TIGHT", "", "_madgraph");

             double sf_cf_weight = (sf_CFrate[0] / (1 - sf_CFrate[0])) + (sf_CFrate[1] / (1 - sf_CFrate[1]));
	     ev_weight = sf_cf_weight;

	     double old_ptsum = lep[0].Pt() + lep[1].Pt();

             for(int i=0; i<2; i++){
               lep[i] = ShiftEnergy( lep[i] , 1 - 0.02 );
	     }
	     HN1 = lep[0] + jets.at(j1) + jets.at(j2);
	     HN2 = lep[1] + jets.at(j1) + jets.at(j2);

             double new_ptsum = lep[0].Pt() + lep[1].Pt();
	     double ptsum_diff = new_ptsum - old_ptsum;

             FillHist("signalEff", 11.5, 1., 0., 15., 15); // jets exist in forward regions
             FillHist("mass_HN1", HN1.M(), ev_weight, 0., 1500., 300);
             FillHist("mass_HN2", HN2.M(), ev_weight, 0., 1500., 300);
             FillHist("pT_e1", electrons.at(0).Pt(), ev_weight, 0., 500., 100);
             FillHist("pT_e2", electrons.at(1).Pt(), ev_weight, 0., 500., 100);
             FillHist("eta_mu1", electrons.at(0).Eta(), ev_weight, -2.5, 2.5, 25);
             FillHist("eta_mu2", electrons.at(1).Eta(), ev_weight, -2.5, 2.5, 25);
             FillHist("MET_Tch", (eventbase->GetEvent().PFMET()-ptsum_diff), ev_weight, 0., 500., 100);
           }
         }
       }
     }
   } 
	    
   return;
}// End of execute event loop
  


void DiEl_SS::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void DiEl_SS::BeginCycle() throw( LQError ){
  
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

DiEl_SS::~DiEl_SS() {
  
  Message("In DiEl_SS Destructor" , INFO);
  
}


void DiEl_SS::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void DiEl_SS::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this DiEl_SSCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void DiEl_SS::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


double DiEl_SS::Get2DCFRates(bool apply_sf, double el_pt, double el_eta, TString el_ID, TString halfsample, TString CFsample){

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
        if(el_eta < 1.4442) CFrate *= 0.59943;
        else CFrate *= 0.802627;
      }
      if(el_ID == "ELECTRON_MVA_TIGHT"){
        if(el_eta < 1.4442) CFrate *= 0.804595;
        else CFrate *= 1.04709;
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
        if(el_eta < 1.4442) CFrate *= 0.600469;
        else CFrate *= 0.804442;
      }
      if(el_ID == "ELECTRON_MVA_TIGHT"){
        if(el_eta < 1.4442) CFrate *= 0.80544;
        else CFrate *= 1.04851;
      }
    }
  }


  return CFrate;
}


snu::KParticle DiEl_SS::ShiftEnergy( snu::KParticle old_lep, double shift_rate ){

  double mass = 0.511e-3;
  snu::KParticle new_lep;
  new_lep.SetPtEtaPhiM((shift_rate*old_lep.Pt()), old_lep.Eta(), old_lep.Phi(), mass) ;
  return new_lep;

}
