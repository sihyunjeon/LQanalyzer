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



  if( isData ){
    CFvalidation();
    return;
  }

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


void CFRateCalculator_FR::CFvalidation(void){

  bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if( !pass_trig ) return;

  double Z_mass = 91.1876;

  for(int aaa=0; aaa<2; aaa++){
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
    if(electronLooseColl.size() != 2) return;
    std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_VETO", false);
    if( muonLooseColl.size() != 0) return;

    FillHist("[CHECK]n_of_muons"+IDsuffix, muonLooseColl.size(), 1., 0., 5., 5);//check no muons
    FillHist("[CHECK]n_of_electrons"+IDsuffix, electronTightColl.size(), 1., 0., 5., 5);//check two electrons

    // define leptons and give Pt, MET cuts
    snu::KParticle lep[2];
    lep[0] = electronLooseColl.at(0);
    lep[1] = electronLooseColl.at(1);

    bool is_SS = false;
    if( (lep[0].Charge() == lep[1].Charge()) ) is_SS = true;

    if( lep[0].Pt() < 25 || lep[1].Pt() < 25 ) return;

    double METPt = eventbase->GetEvent().MET();
    double METPhi = eventbase->GetEvent().METPhi();
    if(METPt > 30) return;

    FillHist("[CHECK]MET"+IDsuffix, METPt, 1., 0., 100., 100);

    bool is_region[2][2] = {{false,},};
    if( (fabs(lep[0].Eta()) < 1.4442) )                                   is_region[0][0] = true;
    else if( (fabs(lep[0].Eta()) > 1.556) && (fabs(lep[0].Eta()) < 2.5) ) is_region[1][0] = true;
    if( (fabs(lep[1].Eta()) < 1.4442) )                                   is_region[0][1] = true;
    else if( (fabs(lep[1].Eta()) > 1.556) && (fabs(lep[1].Eta()) < 2.5) ) is_region[1][1] = true;

    TString region = "";
    if((is_region[0][0] && is_region[0][1]))      region = "BB";
    else if((is_region[1][0] && is_region[1][1])) region = "EE";
    else                                          region = "BE";

    // reduce energy of leptons if OS because of photon radiation
    // increase energy of leptons if SS for better fitting using Gaussian

    snu::KParticle Z_candidate;//define Z candidate after shifting energy (SS : better fitting, energy scale up // OS : photon radiation E loss)
    Z_candidate = (lep[0] + lep[1]);
    bool Z_selection = (fabs(Z_candidate.M() - Z_mass) < 10.);

    std::vector<snu::KJet> jetTightColl = GetJets("JET_HN", 30., 2.4);
    int Njets = jetTightColl.size();
    TString s_njets = "", onejet = "";
    if( Njets == 0 ) s_njets = "JETS0";
    if( Njets != 0 ) s_njets = "JETS";

    if( is_SS ){
      if( Z_selection ){

        double fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 0, electronLooseColl, el_ID, 2, el_looseID, method_fake);

        FillHist("observed_Z_mass_global"+IDsuffix, Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_global"+IDsuffix, 0., fake_weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+region+IDsuffix, Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_"+region+IDsuffix, 0., fake_weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+s_njets+"_global"+IDsuffix, Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_"+s_njets+"_global"+IDsuffix, 0., fake_weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+s_njets+"_"+region+IDsuffix, Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_"+s_njets+"_"+region+IDsuffix, 0., fake_weight, 0., 1., 1);
        if(Njets == 0) FillHist("observed_n_jets_global"+IDsuffix, 0., fake_weight, 0., 2., 2);
        if(Njets != 0) FillHist("observed_n_jets_global"+IDsuffix, 1., fake_weight, 0., 2., 2);
        if(Njets == 1){
          FillHist("observed_Z_mass_JETS1_global"+IDsuffix, Z_candidate.M(), fake_weight, 70., 110., 40);
          FillHist("observed_n_events_JETS1_global"+IDsuffix, 0., fake_weight, 0., 1., 1);
          FillHist("observed_Z_mass_JETS1_"+region+IDsuffix, Z_candidate.M(), fake_weight, 70., 110., 40);
          FillHist("observed_n_events_JETS1_"+region+IDsuffix, 0., fake_weight, 0., 1., 1);
        }
      }

      for(int sss=0; sss<2; sss++){
        snu::KParticle Z_candidate_SHIFTED;
        snu::KParticle lep_SHIFTED[2];
	double shiftrate = -999.;
	TString CFsample = "";

        if(sss==0){
          if(el_ID=="ELECTRON_HN_TIGHT") shiftrate = (1.-0.009);
          if(el_ID=="ELECTRON_MVA_TIGHT") shiftrate = (1.-0.014);
	  CFsample = "_powheg";
        }
        if(sss==1){
          if(el_ID=="ELECTRON_HN_TIGHT") shiftrate = (1.-0.019);
          if(el_ID=="ELECTRON_MVA_TIGHT") shiftrate = (1.-0.016);
	  CFsample = "_madgraph";
        }

        for(int i=0; i<2; i++){
          if( !is_SS ){
            lep_SHIFTED[i] = ShiftEnergy(lep[i], 1./shiftrate);
            FillHist("[CHECK]lepton_SHIFTED_E_shift_up"+IDsuffix+CFsample, lep_SHIFTED[i].E(), 1., 0., 400., 400);
          }
        }
        Z_candidate_SHIFTED = lep_SHIFTED[0] + lep_SHIFTED[1];

        if(fabs(Z_candidate_SHIFTED.M() - Z_mass) < 10.){
          double fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 0, electronLooseColl, el_ID, 2, el_looseID, method_fake);
  
          FillHist("SHIFTED_Z_mass_global"+IDsuffix+CFsample, Z_candidate_SHIFTED.M(), fake_weight, 70., 110., 40);
        }
      }

    }
  }
  return;
}


snu::KParticle CFRateCalculator_FR::ShiftEnergy( snu::KParticle old_lep, double shift_rate ){

  double mass = 0.511e-3;
  snu::KParticle new_lep;
  new_lep.SetPtEtaPhiM((shift_rate*old_lep.Pt()), old_lep.Eta(), old_lep.Phi(), mass) ;
  return new_lep;

}
