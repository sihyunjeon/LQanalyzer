// $Id: ExampleAnalyzer.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQExampleAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ExampleAnalyzer.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ExampleAnalyzer);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ExampleAnalyzer::ExampleAnalyzer() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ExampleAnalyzer");
  
  Message("In ExampleAnalyzer constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void ExampleAnalyzer::InitialiseAnalysis() throw( LQError ) {
  
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


void ExampleAnalyzer::ExecuteEvents()throw( LQError ){

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
   if (!k_isdata) {   pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());}

  std::vector<TString> triggerlist_emBG1, triggerlist_emBG2, triggerlist_emH1, triggerlist_emH2;

  triggerlist_emBG1.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
  triggerlist_emBG2.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v");

  triggerlist_emH1.push_back("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
  triggerlist_emH2.push_back("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v");
    
  bool Pass_Trigger_em = (PassTriggerOR(triggerlist_emBG1) || PassTriggerOR(triggerlist_emBG2) || PassTriggerOR(triggerlist_emH1) || PassTriggerOR(triggerlist_emH2)); 
  if(!Pass_Trigger_em) return;

  std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_HN_VETO", false);
  std::vector<snu::KMuon> muons;
  std::vector<snu::KMuon> muonTightColl;
  muonTightColl.clear(); muons.clear();
  int muonVetoN = muonVetoColl.size(), muonTightN = 0, muonsN = 0;
  for(unsigned int i=0;i<muonVetoN;i++){
    if(PassID(muonVetoColl.at(i), "MUON_HN_LOOSE")){
      muons.push_back(muonVetoColl.at(i));
      muonsN ++;
      if(PassID(muonVetoColl.at(i), "MUON_HN_TIGHT")){
        muonTightColl.push_back(muonVetoColl.at(i));
        muonTightN++;
      }
    }
  }

  std::vector<snu::KElectron> electronVetoColl = GetElectrons(false,false,"ELECTRON_HN_VETO");
  std::vector<snu::KElectron> electrons;
  std::vector<snu::KElectron> electronTightColl;
  electronTightColl.clear(); electrons.clear();
  int electronVetoN = electronVetoColl.size(), electronTightN = 0, electronsN = 0;
  for(unsigned int i=0;i<electronVetoN;i++){
    if(PassID(electronVetoColl.at(i), "ELECTRON_HN_FAKELOOSE")){
      electrons.push_back(electronVetoColl.at(i));
      electronsN ++;
      if(PassID(electronVetoColl.at(i), "ELECTRON_HN_TIGHTv4")){
        electronTightColl.push_back(electronVetoColl.at(i));
        electronTightN++;
      }
    }
  }

  std::vector<snu::KJet> jets = GetJets("JET_HN");
  std::vector<snu::KFatJet> fatjets = GetFatJets("FATJET_HN");
  std::vector<snu::KJet> jets_for_bjet = GetJets("JET_NOLEPTONVETO");
  std::vector<snu::KJet> Tjets = GetJets("JET_HN_TChannel", 20.);
  std::vector<snu::KJet> bjets, bjetsloose;
  for(unsigned int j=0; j<jets_for_bjet.size(); j++){
    if(jets_for_bjet.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Loose)){
      bjetsloose.push_back(jets_for_bjet.at(j));
      if(jets_for_bjet.at(j).IsBTagged(snu::KJet::CSVv2, snu::KJet::Medium)) bjets.push_back(jets_for_bjet.at(j));
    }
  }

 
  if(fatjets.size() <1 && jets.size() <2) return;
  if(bjets.size() != 0 ) return;
  if(muonVetoN != 1) return;
  if(muonsN != 1) return;
  if(muonTightN != 1) return;
  if(electronVetoN != 1) return;
  if(electronsN != 1) return;
  if(electronTightN != 1) return;
  if(muons.at(0).Charge() != electrons.at(0).Charge()) return;
  if(!PassEMuTriggerPt(electrons, muons)) return;

<<<<<<< HEAD
   mcdata_correction->CorrectMuonMomentum(muons,eventbase->GetTruth()); /// CorrectMuonMomentum(muons);  will also work as Funcion in AnalyzerCore just calls mcdata_correction function
   
   double ev_weight = weight;
   if(!isData){
     //ev_weight = w * trigger_sf * id_iso_sf *  pu_reweight*trigger_ps;
   }

   if(jets.size() > 3){
     if(nbjet > 0){
       if(muons.size() ==2) {
	 if(electrons.size() >= 1){
	   cout << "electrons is tight = " << electrons.at(0).PassTight() << endl;
	   if(!SameCharge(muons)){
	     if(muons.at(0).Pt() > 20. && muons.at(1).Pt() > 10.){
	       if(eventbase->GetEvent().PFMET() > 30){
		 if(trig_pass){
		   FillHist("Massmumu", GetDiLepMass(muons), ev_weight, 0., 200.,400);
		   FillHist("Massmumu_zoomed", GetDiLepMass(muons), ev_weight, 0.,50.,200);
		   FillCLHist(sighist_mm, "DiMuon", eventbase->GetEvent(), muons,electrons,jets, ev_weight);
		 }
	       }
	     }
	   }
	 }
       }
     }
   }

   	    
   float cf_weight = -999.;
   if(CFelectrons.size() == 2){
     if(CFelectrons.at(0).Charge() != CFelectrons.at(1).Charge()){//CF estimation is from OS dielectrons
       cf_weight = GetCFweight(0,CFelectrons, true, "ELECTRON_HN_TIGHTv4");//put in syst=1,0,-1,  electronColl vector, apply sf, electron ID
         //scale factor up downs with syst = 1, -1

       std::vector<snu::KElectron> CFelectrons_shifted = ShiftElectronEnergy(CFelectrons, "ELECTRON_HN_TIGHTv4", true);//apply pt shift considering brem radiation after getting CF weights
       if(CFelectrons_shifted.at(1).Pt()>20){//apply pt cuts again after shifting pt
         FillHist("chargeflipped_Z_mass", (CFelectrons.at(0)+CFelectrons.at(1)).M(),cf_weight, 50., 130., 80);
         FillHist("chargeflipped_leading_lepton", CFelectrons.at(0).Pt(), cf_weight, 0., 200., 200);
   
       }
     }
   }
=======
cout <<eventbase->GetEvent().EventNumber()<<endl;
>>>>>>> 13TeV_v8-0-7.24_AUG2017

  FillHist("Nevents", 0., 1., 0., 1., 1);

  return;
}// End of execute event loop
  


void ExampleAnalyzer::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ExampleAnalyzer::BeginCycle() throw( LQError ){
  
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

ExampleAnalyzer::~ExampleAnalyzer() {
  
  Message("In ExampleAnalyzer Destructor" , INFO);
  
}


void ExampleAnalyzer::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ExampleAnalyzer::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ExampleAnalyzerCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void ExampleAnalyzer::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


bool ExampleAnalyzer::PassEMuTriggerPt(std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons){

  bool pass =false;
  snu::KParticle el,mu;
  el = electrons.at(0);
  mu = muons.at(0);

  if(PassTriggerOR(triggerlist_emBG1)){ pass = ((mu.Pt() >10 && el.Pt() >25)); }
  if(PassTriggerOR(triggerlist_emBG2)){ pass = ((mu.Pt() >25 && el.Pt() >10)); }
  if(PassTriggerOR(triggerlist_emH1)){ pass = ((mu.Pt() >10 && el.Pt() >25)); }
  if(PassTriggerOR(triggerlist_emH2)){ pass = ((mu.Pt() >25 && el.Pt() >10)); }

  return pass;
}
