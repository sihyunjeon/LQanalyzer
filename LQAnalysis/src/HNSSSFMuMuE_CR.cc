// $Id: HNSSSFMuMuE_CR.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNSSSFMuMuE_CR Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNSSSFMuMuE_CR.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNSSSFMuMuE_CR);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNSSSFMuMuE_CR::HNSSSFMuMuE_CR() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNSSSFMuMuE_CR");
  
  Message("In HNSSSFMuMuE_CR constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNSSSFMuMuE_CR::InitialiseAnalysis() throw( LQError ) {
  
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


void HNSSSFMuMuE_CR::ExecuteEvents()throw( LQError ){

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
  if (!k_isdata) {   pileup_reweight = TempPileupWeight();}
    
  TString mumue_trigger="HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v";
  vector<TString> trignames;
  trignames.push_back(mumue_trigger);
  float weight_trigger = WeightByTrigger("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v", TargetLumi);

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_NOCUT");
  std::vector<snu::KJet> jetTightColl = GetJets("JET_HN");
  int nbjet = NBJet(GetJets("JET_HN"));

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_TRI_LOOSE",false);
  std::vector<snu::KMuon> muonTightColl = GetMuons("MUON_HN_TRI_TIGHT",false);

  std::vector<snu::KElectron> electronLooseColl = GetElectrons("ELECTRON_HN_FAKELOOSE", false);
  std::vector<snu::KElectron> electronTightColl = GetElectrons("ELECTRON_HN_TIGHT", false);

  bool trig_pass=PassTrigger("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v");
  if(!trig_pass) return;
  CorrectMuonMomentum(muonLooseColl);
   
  double ev_weight = weight;
  if(!isData){
    weight *= weight_trigger;
    weight *= pileup_reweight;
    //ev_weight = w * trigger_sf * id_iso_sf *  pu_reweight*trigger_ps;
  }

  for(int i=0; i<muonTightColl.size(); i++)
    if(muonTightColl.at(i).Pt() < 10) return;

  for(int i=0; i<electronTightColl.size(); i++)
    if(electronTightColl.at(i).Pt() < 10) return;

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();

  snu::KParticle MET;
  MET.SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), 0., METPt);

  if( ((electronLooseColl.size() == 1) && (electronTightColl.size() == 1))
   && ((muonLooseColl.size() == 2) && (muonTightColl.size() == 2)) ){

    snu::KParticle RAWmu[2], RAWel;
    RAWmu[0] = muonLooseColl.at(0);
    RAWmu[1] = muonLooseColl.at(1);
    RAWel = electronLooseColl.at(0);

    if( RAWmu[0].Charge() != RAWmu[1].Charge() ){

      snu::KParticle RAWnu[2], RECOnu;

      snu::KParticle W_lepton;
      W_lepton = RAWel;

      double nuPz = 999.;
      nuPz = CalculateNuPz(W_lepton, MET, 1);
      RAWnu[0].SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), nuPz, TMath::Sqrt( METPt*METPt + nuPz*nuPz ));
      nuPz = CalculateNuPz(W_lepton, MET, -1);
      RAWnu[1].SetPxPyPzE(METPt*(TMath::Cos(METPhi)), METPt*(TMath::Sin(METPhi)), nuPz, TMath::Sqrt( METPt*METPt + nuPz*nuPz ));
      if( fabs(RAWnu[0].Pz()) < fabs(RAWnu[1].Pz()) ){
        RECOnu = RAWnu[0];
      }
      else RECOnu = RAWnu[1];

      EventSelectionStudy( RAWmu, RAWel, jetTightColl, RECOnu, 1, weight );
    }

    else{
      if( RAWmu[0].Charge() == RAWel.Charge() ){

	FillHist("AllSameSignMuMuE", 0., weight, 0., 1., 1);

      }
    }
  }



 if((muonLooseColl.size() == 2) && (muonTightColl.size() == 2)){

      

  }
   



/*
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
  if( abs(RAWnu[0].Pz()) < abs(RAWnu[1].Pz()) ){
    RECOnu_lowmass = RAWnu[0];
  }
  else RECOnu_lowmass = RAWnu[1];

  RECOW_pri_lowmass = RAWmu[0] + RAWmu[1] + RAWel + RECOnu_lowmass;
  RECOW_sec_lowmass = RAWel + RECOnu_lowmass;

  cout << RECOW_pri_lowmass.M() << endl;

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
  if( abs(RAWnu[0].Pz()) < abs(RAWnu[1].Pz()) ){
    RECOnu_highmass = RAWnu[0];
  }
  else RECOnu_highmass = RAWnu[1];

  RECOW_pri_highmass = RAWmu[0] + RAWmu[1] + RAWel + RECOnu_highmass;
  RECOW_sec_highmass = RAWel + RECOnu_highmass;

//cout << RECOW_sec_highmass.M() << endl;

  // ========== CLASS 3 =====================================
  EventSelectionStudy(RAWmu, RAWel, 3);// RECO particles output
  RECOHN[2] = RECOmu[1] + RECOel + RECOnu_highmass;

  // ========== CLASS 4 =====================================
  EventSelectionStudy(RAWmu, RAWel, 4);// RECO particles output
  RECOHN[3] = RECOmu[1] + RECOel + RECOnu_highmass;


  FillHist("number_of_events_cut0", 0., weight, 0., 1., 1);
  FillHist("W_primary_lowmass_cut0", RECOW_pri_lowmass.M(), weight, 0., 1000., 1000);
  FillHist("W_secondary_lowmass_cut0", RECOW_sec_lowmass.M(), weight, 0., 1000., 1000);
  FillHist("W_primary_highmass_cut0", RECOW_pri_highmass.M(), weight, 0., 1000., 1000);
  FillHist("W_secondary_highmass_cut0", RECOW_sec_highmass.M(), weight, 0., 1000., 1000);
  FillHist("HN_mass_class1_cut0", RECOHN[0].M(), weight, 0., 200., 200);
  FillHist("HN_mass_class2_cut0", RECOHN[1].M(), weight, 0., 200., 200);
  FillHist("HN_mass_class3_cut0", RECOHN[2].M(), weight, 0., 800., 800);
  FillHist("HN_mass_class4_cut0", RECOHN[3].M(), weight, 0., 1500., 1500);
  FillCLHist(sssf_mumue, "cut0", eventbase->GetEvent(), muonLooseColl, electronLooseColl, jetTightColl, weight);

*/
  return;
}// End of execute event loop
  


void HNSSSFMuMuE_CR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNSSSFMuMuE_CR::BeginCycle() throw( LQError ){
  
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

HNSSSFMuMuE_CR::~HNSSSFMuMuE_CR() {
  
  Message("In HNSSSFMuMuE_CR Destructor" , INFO);
  
}


void HNSSSFMuMuE_CR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNSSSFMuMuE_CR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNSSSFMuMuE_CRCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void HNSSSFMuMuE_CR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


void HNSSSFMuMuE_CR::EventSelectionStudy( snu::KParticle RAWmu[], snu::KParticle RAWel, std::vector<snu::KJet> jetColl, snu::KParticle RECOnu, int control_region, double weight ){

  int N_bjet = 0;
  for( unsigned int i = 0; i < jetColl.size() ; i++ ){
    if( jetColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight) ){
      if( (jetColl.at(i).Pt()) > 20 && fabs(jetColl.at(i).Eta()) < 2.4 ){
	N_bjet ++;
      }
    }
  }
 
  StudyControlRegion( "cut0", RAWmu, RAWel, jetColl, RECOnu, 0, weight ) ;
 
  if( control_region == 1 ){

    snu::KParticle RECOmu[2], RECOel;
    snu::KParticle Z_candidate;
    Z_candidate = RAWmu[0] + RAWmu[1];

    if( Z_candidate.M() > 76 && Z_candidate.M() < 106 ){
      StudyControlRegion( "cutZmass", RAWmu, RAWel, jetColl, RECOnu, control_region, weight ) ;

      if( RAWmu[0].Pt() > 20 && RAWmu[1].Pt() > 10 ){
	RECOmu[0] = RAWmu[0];
	RECOmu[1] = RAWmu[1];

	if( RAWel.Pt() > 20 ){
	  RECOel = RAWel;
	  StudyControlRegion( "cutZmass_cutLeptonsPt", RAWmu, RAWel, jetColl, RECOnu, control_region, weight ) ;

	  if( RECOnu.Pt() > 30 ){
            StudyControlRegion( "cutZmass_cutLeptonsPt_cutMET", RAWmu, RAWel, jetColl, RECOnu, control_region, weight ) ;

	    if( (RECOmu[0]+RECOmu[1]).M() > 4 && (RECOmu[0]+RECOel).M() > 4 && (RECOmu[1]+RECOel).M() > 4){
	      if( (RECOmu[0]+RECOmu[1]+RECOel).M() > 100 ){
		StudyControlRegion( "cutZmass_cutLeptonsPt_cutMET_cutLeptonsMass", RAWmu, RAWel, jetColl, RECOnu, control_region, weight ) ;		

		if( N_bjet == 0 ){
		  StudyControlRegion( "cutZmass_cutLeptonsPt_cutMET_cutLeptonsMass_cutBJet", RAWmu, RAWel, jetColl, RECOnu, control_region, weight ) ;

		  if( jetColl.size() == 0 ){
		    StudyControlRegion( "CR1_0Jets", RAWmu, RAWel, jetColl, RECOnu, control_region, weight ) ;}
		  else if( jetColl.size() == 1 ){
                    StudyControlRegion( "CR1_1Jets", RAWmu, RAWel, jetColl, RECOnu, control_region, weight ) ;} 
		  else if( jetColl.size() == 2 ){
                    StudyControlRegion( "CR1_2Jets", RAWmu, RAWel, jetColl, RECOnu, control_region, weight ) ;}
                  else{
                    StudyControlRegion( "CR1_MoreThan3Jets", RAWmu, RAWel, jetColl, RECOnu, control_region, weight ) ;}
		}
	      }
	    }
	  }
	}
      }
    }
  }

  else if( control_region == 2 ){



  }

  return;
}


void HNSSSFMuMuE_CR::StudyControlRegion( TString suffix, snu::KParticle RAWmu[], snu::KParticle RAWel, std::vector<snu::KJet> jetColl, snu::KParticle RECOnu, int control_region, double weight ){

  int N_bjet = 0;
  for( unsigned int i = 0; i < jetColl.size() ; i++ ){
    if( jetColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight) ){
      if( (jetColl.at(i).Pt()) > 20 && fabs(jetColl.at(i).Eta()) < 2.4 ){
        N_bjet ++;
      }
    }
  }

  TString prefix;
  if( control_region == 0 ) prefix = "ControlRegion0";
  else if( control_region == 1 ) prefix = "ControlRegion1";
  else if( control_region == 2 ) prefix = "ControlRegion2";

  FillHist(prefix+"_Z-LeadingMuonPt_"+suffix, RAWmu[0].Pt(), weight, 0., 500., 500);
  FillHist(prefix+"_Z-LeadingMuonEta_"+suffix, RAWmu[0].Eta(), weight, -3., 3., 60);
  FillHist(prefix+"_Z-SecondMuonPt_"+suffix, RAWmu[1].Pt(), weight, 0., 500., 500);
  FillHist(prefix+"_Z-SecondMuonEta_"+suffix, RAWmu[1].Eta(), weight, -3., 3., 60);
  FillHist(prefix+"_Z-CandidateMass_"+suffix, (RAWmu[0]+RAWmu[1]).M(), weight, 0., 200., 200);
  FillHist(prefix+"_W-CandidateTransverseMass_"+suffix, GetTransverseMass( RAWel, RECOnu ), weight, 0., 200., 200);
  FillHist(prefix+"_TrileptonMass_"+suffix, (RAWmu[0]+RAWmu[1]+RAWel).M(), weight, 0., 500., 500);
  FillHist(prefix+"_DileptonMass_"+suffix, (RAWmu[0]+RAWmu[1]).M(), weight, 0., 50., 500);
  FillHist(prefix+"_DileptonMass_"+suffix, (RAWmu[0]+RAWel).M(), weight, 0., 500., 500);
  FillHist(prefix+"_DileptonMass_"+suffix, (RAWmu[1]+RAWel).M(), weight, 0., 500., 500);
  FillHist(prefix+"_NJets_"+suffix, jetColl.size(), weight, 0., 5., 5);
  FillHist(prefix+"_NBJets_"+suffix, N_bjet, weight, 0., 5., 5);
  for(unsigned int i=0; i<jetColl.size(); i++){
    FillHist(prefix+"_JetPt_"+suffix, jetColl.at(i).Pt(), weight, 0., 500., 500);
    FillHist(prefix+"_JetEta_"+suffix, jetColl.at(i).Eta(), weight, -3., 3., 60);
    if( jetColl.at(i).IsBTagged(snu::KJet::CSVv2, snu::KJet::Tight) ){
      FillHist(prefix+"_BJetPt_"+suffix, jetColl.at(i).Pt(), weight, 0., 500., 500);
      FillHist(prefix+"_BJetEta_"+suffix, jetColl.at(i).Eta(), weight, -3., 3., 60);
    }
  }
  FillHist(prefix+"_MET_"+suffix, RECOnu.Pt(), weight, 0., 500., 500);

}
