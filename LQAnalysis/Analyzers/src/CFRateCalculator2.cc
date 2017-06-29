// $Id: CFRateCalculator2.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQCFRateCalculator2 Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "CFRateCalculator2.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (CFRateCalculator2);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
CFRateCalculator2::CFRateCalculator2() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("CFRateCalculator2");
  
  Message("In CFRateCalculator2 constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void CFRateCalculator2::InitialiseAnalysis() throw( LQError ) {
  
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


void CFRateCalculator2::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  if(!PassMETFilter()) return;     /// Initial event cuts :
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex

  bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if( !pass_trig ) return;

  double MET = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  if(MET > 30) return;

  std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_HN_VETO",false);
  if((muonVetoColl.size() != 0)) return;

  double Z_mass = 91.1876;

  float pileup_reweight=(1.0);
  float trigger_reweight=(1.0);

  if(!k_isdata){
    pileup_reweight = mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    trigger_reweight = WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);

    weight *=pileup_reweight;
    weight *=trigger_reweight;
  }



  for(int aaa=0; aaa<2; aaa++){

    std::vector<snu::KElectron> electronTightColl;

    TString IDsuffix = "";
    TString el_ID = "", el_looseID = "";
    TString method_fake = "";

    if(aaa == 0){
      el_ID = "ELECTRON_HN_TIGHT";
      el_looseID = "ELECTRON_HN_FAKELOOSE";
      IDsuffix = "_HNTIGHT";
      method_fake = "mva";
    }
    if(aaa == 1){
      el_ID = "ELECTRON_MVA_TIGHT";
      el_looseID = "ELECTRON_MVA_FAKELOOSE";
      IDsuffix = "_MVATIGHT";
      method_fake = "dijet_ajet40";
    }  

    electronTightColl = GetElectrons(false, false, el_ID);
    if(k_running_nonprompt) electronTightColl = GetElectrons(false, false, el_looseID);
    if((electronTightColl.size() != 2)) continue;
    if((electronTightColl.at(1).Pt() < 25)) continue;
 
    if(!k_isdata){
      std::vector<snu::KElectron> electronLooseColl = GetElectrons(false, false, el_looseID);

      double electron_idsf = mcdata_correction->ElectronScaleFactor(el_ID, electronLooseColl);
      double electron_recosf = mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);

      weight*=electron_idsf;
      weight*=electron_recosf;
    }

    snu::KElectron lep[2];
    for(int i=0; i<2; i++){
      lep[i] = electronTightColl.at(i);
      if((fabs(lep[i].Eta()) < 1.556) && (fabs(lep[i].Eta()) > 1.4442)) continue;
    }
     
    snu::KParticle Z_candidate = lep[0]+lep[1];
    bool is_Z = ((75 < fabs(Z_candidate.M() - Z_mass)) && (fabs(Z_candidate.M() - Z_mass) < 100));
    if(!is_Z) continue;    

    TString charge_suffix = "_OS";
    bool is_SS = false;
    if((lep[0].Charge() == lep[1].Charge())){
      charge_suffix = "_SS";
      is_SS = true;
    }

    double fake_weight=1.0;
    if(k_running_nonprompt){
      if(is_SS){
        fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonVetoColl, "MUON_HN_TRI_TIGHT", 0, electronTightColl, el_ID, 2, el_looseID, method_fake);
        weight = fake_weight;
      }
      else continue;
    }

/*    bool is_region[2][2] = {{false,},};
    if( (fabs(lep[0].Eta()) < 1.4442) )                                   is_region[0][0] = true;
    else if( (fabs(lep[0].Eta()) > 1.556) && (fabs(lep[0].Eta()) < 2.5) ) is_region[1][0] = true;
    if( (fabs(lep[1].Eta()) < 1.4442) )                                   is_region[0][1] = true;
    else if( (fabs(lep[1].Eta()) > 1.556) && (fabs(lep[1].Eta()) < 2.5) ) is_region[1][1] = true;

    TString region = "";
    if((is_region[0][0] && is_region[0][1]))      region = "BB";
    else if((is_region[1][0] && is_region[1][1])) region = "EE";
    else                                          region = "BE";*/

    float etaarray [] = {0.0, 0.9, 1.4442, 1.556, 2.5};
    float ptarray [] = {20., 40., 60., 80., 100., 200., 500.};

    for(int i=0; i<2; i++){
      FillHist("Pt_eta_global"+IDsuffix+charge_suffix, fabs(lep[i].Eta()), lep[i].Pt(), weight, etaarray, 4, ptarray, 6);

      FillHist("n_events_global"+IDsuffix+charge_suffix, 0., weight, 0., 1., 1);
      FillHist("Pt_global"+IDsuffix+charge_suffix, lep[i].Pt(), weight, 0., 500, 500);
      FillHist("invPt_global"+IDsuffix+charge_suffix, 1./lep[i].Pt(), weight, 0., 0.04, 20);
      FillHist("eta_global"+IDsuffix+charge_suffix, lep[i].Eta(), weight, -3., 3., 120);
      FillHist("dXY_global"+IDsuffix+charge_suffix, fabs(lep[i].dxy()), weight, 0., 0.02, 100);

      TString s_region = "";
      if( (fabs(lep[i].Eta()) < 0.9) ) s_region = "IB";
      else if( (fabs(lep[i].Eta()) < 1.4442) ) s_region = "OB";
      else if( (fabs(lep[i].Eta()) < 1.556) ) s_region = "EC";

      FillHist("n_events_"+s_region+IDsuffix+charge_suffix, 0., weight, 0., 1., 1);
      FillHist("Pt_"+s_region+IDsuffix+charge_suffix, lep[i].Pt(), weight, 0., 500, 500);
      FillHist("invPt_"+s_region+IDsuffix+charge_suffix, 1./lep[i].Pt(), weight, 0., 0.04, 20);
      FillHist("eta_"+s_region+IDsuffix+charge_suffix, lep[i].Eta(), weight, -3., 3., 120);
      FillHist("dXY_"+s_region+IDsuffix+charge_suffix, fabs(lep[i].dxy()), weight, 0., 0.02, 100);
    }
    FillHist("MET_global"+IDsuffix+charge_suffix, MET, weight, 0., 1000., 1000);
  }
  return;
}// End of execute event loop
  


void CFRateCalculator2::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void CFRateCalculator2::BeginCycle() throw( LQError ){
  
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

CFRateCalculator2::~CFRateCalculator2() {
  
  Message("In CFRateCalculator2 Destructor" , INFO);
  
}


void CFRateCalculator2::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void CFRateCalculator2::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this CFRateCalculator2Core::MakeHistograms() to make new hists for your analysis
   **/
  
}


void CFRateCalculator2::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

