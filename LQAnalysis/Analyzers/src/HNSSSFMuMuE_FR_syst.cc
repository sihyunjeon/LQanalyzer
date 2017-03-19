// $Id: HNSSSFMuMuE_FR_syst.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNSSSFMuMuE_FR_syst Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNSSSFMuMuE_FR_syst.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNSSSFMuMuE_FR_syst);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNSSSFMuMuE_FR_syst::HNSSSFMuMuE_FR_syst() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNSSSFMuMuE_FR_syst");
  
  Message("In HNSSSFMuMuE_FR_syst constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNSSSFMuMuE_FR_syst::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  Message("Making clever hists for Z ->ll test code", INFO);

  return;
}


void HNSSSFMuMuE_FR_syst::ExecuteEvents()throw( LQError ){

  // ========== Define RelIso ====================
  double this_reliso_mu = 0.4;
  double this_reliso_el = 0.5; 
  // ================================================================================


  // ========== No cut ====================
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  FillCutFlow("NoCut", weight);
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  // ================================================================================


  // ========== MET filter cut ====================
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", weight);
  // ================================================================================


  // ========== Primary vertex cut ====================
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex 
  // ================================================================================


  // ========== Trigger cut ====================
  TString mumu_trigger="HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v";
  vector<TString> trignames;
  trignames.push_back(mumu_trigger);

  bool trig_pass=PassTrigger("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v");
  if(!trig_pass) return;
  // ================================================================================


  // ========== Get Objects (muon, electron, jet) ====================
  std::vector<snu::KMuon> muonVLooseColl = GetMuons("MUON_HN_TRI_VLOOSE",false);
  std::vector<snu::KMuon> muonLooseColl;
  muonLooseColl.clear();

  std::vector<snu::KElectron> electronVLooseColl = GetElectrons(false,false,"ELECTRON16_HN_FAKEVLOOSE");
  std::vector<snu::KElectron> electronLooseColl;
  electronLooseColl.clear();

  std::vector<snu::KJet> jetVLooseColl = GetJets("JET_HN", 30., 2.4);
  // ================================================================================


  /*####################################################################################################
  ##		        Analysis Code 								      ##
  ##				For SameSign MuMuE Channel Analysis				      ##
  ####################################################################################################*/


  std::vector<double> RelIso_array;
  float RelIso_values[7] = {0.2,0.3,0.4,0.6,0.8,1.0,1.1};
  for(int i=0; i<7; i++){
    RelIso_array.push_back(RelIso_values[i]);
  }

  std::vector<double> dXYSig_array;
  float dXYSig_values[4] = {3,4,5,6};
  for(int i=0; i<4; i++){
    dXYSig_array.push_back(dXYSig_values[i]);
  }

  std::vector<double> awayJetPt_array;
  float awayJetPt_values[5] = {20,30,40,60,100};
  for(int i=0;i<5; i++){
    awayJetPt_array.push_back(awayJetPt_values[i]);
  }

  double weight_err = -999.;

  for(int aaa=0; aaa<(RelIso_array.size()-1); aaa++){

    muonLooseColl.clear(); electronLooseColl.clear();

    for(int mu_it=0; mu_it<muonVLooseColl.size(); mu_it++){
      if((muonVLooseColl.at(mu_it).RelIso04() < RelIso_array.at(aaa))){
        muonLooseColl.push_back(muonVLooseColl.at(mu_it));
      }
    }
    for(int el_it=0; el_it<electronVLooseColl.size(); el_it++){
      if((electronVLooseColl.at(el_it).PFRelIso(0.3) < 0.5)){
        electronLooseColl.push_back(electronVLooseColl.at(el_it));
      }
    }
      
    if( !((muonLooseColl.size()==2) && (electronLooseColl.size() == 1)) ) continue;

    snu::KParticle lep[3];
    lep[0] = muonLooseColl.at(0);
    lep[1] = muonLooseColl.at(1);
    lep[2] = electronLooseColl.at(0);
  
    if( (lep[0].Charge() != lep[1].Charge()) ) continue;
    if( (lep[0].Charge() == lep[2].Charge()) ) continue;
    if( (lep[1].Charge() == lep[2].Charge()) ) continue;

    if( (lep[0].Pt() < 20) || (lep[1].Pt() < 10) || (lep[2].Pt() < 10) ) continue;


    for(int bbb=0; bbb<(dXYSig_array.size()-1); bbb++){

      m_datadriven_bkg->GetFakeObj()->SetTrilepWP(dXYSig_array.at(bbb), RelIso_array.at(aaa));
      m_datadriven_bkg->GetFakeObj()->SetTrilepElWP(40);

      weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON16_HN_TIGHT", 1);
      weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON16_HN_TIGHT", 1);

      FillHist("N_of_Fake_Events_Differ_Muon", dXYSig_array.at(bbb), RelIso_array.at(aaa), weight, dXYSig_values, (dXYSig_array.size()-1), RelIso_values, (RelIso_array.size()-1));
      FillHist("N_of_Fake_Events_Differ_Muon_up", dXYSig_array.at(bbb), RelIso_array.at(aaa), weight+weight_err, dXYSig_values, (dXYSig_array.size()-1), RelIso_values, (RelIso_array.size()-1));
      FillHist("N_of_Fake_Events_Differ_Muon_down", dXYSig_array.at(bbb), RelIso_array.at(aaa), weight-weight_err, dXYSig_values, (dXYSig_array.size()-1), RelIso_values, (RelIso_array.size()-1));

      FillHist("N_Yield_Differ_Muon_"+GetSuffix(RelIso_array, aaa, dXYSig_array, bbb), 0., weight, 0., 1., 1);
      FillHist("N_Yield_Differ_Muon_"+GetSuffix(RelIso_array, aaa, dXYSig_array, bbb)+"_up", 0., weight+weight_err, 0., 1., 1);
      FillHist("N_Yield_Differ_Muon_"+GetSuffix(RelIso_array, aaa, dXYSig_array, bbb)+"_down", 0., weight-weight_err, 0., 1., 1);



    }
  }


  for(int aaa=0; aaa<(awayJetPt_array.size()-1); aaa++){

    muonLooseColl.clear(); electronLooseColl.clear();

    for(int mu_it=0; mu_it<muonVLooseColl.size(); mu_it++){
      if((muonVLooseColl.at(mu_it).RelIso04() < 0.4) && (muonVLooseColl.at(mu_it).dXYSig() < 3)){
        muonLooseColl.push_back(muonVLooseColl.at(mu_it));
      }
    }
    for(int el_it=0; el_it<electronVLooseColl.size(); el_it++){
      if((electronVLooseColl.at(el_it).PFRelIso(0.3) < 0.5)){
        electronLooseColl.push_back(electronVLooseColl.at(el_it));
      }
    }

    if( !((muonLooseColl.size()==2) && (electronLooseColl.size() == 1)) ) continue;

    snu::KParticle lep[3];
    lep[0] = muonLooseColl.at(0);
    lep[1] = muonLooseColl.at(1);
    lep[2] = electronLooseColl.at(0);

    if( (lep[0].Charge() != lep[1].Charge()) ) continue;
    if( (lep[0].Charge() == lep[2].Charge()) ) continue;
    if( (lep[1].Charge() == lep[2].Charge()) ) continue;

    if( (lep[0].Pt() < 20) || (lep[1].Pt() < 10) || (lep[2].Pt() < 10) ) continue;

    m_datadriven_bkg->GetFakeObj()->SetTrilepWP(4.0, 0.4);
    m_datadriven_bkg->GetFakeObj()->SetTrilepElWP(awayJetPt_array.at(aaa));

    weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON16_HN_TIGHT", 1);
    weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON16_HN_TIGHT", 1);

    FillHist("N_of_Fake_Events_Differ_Electron", awayJetPt_array.at(aaa), weight, awayJetPt_values, (awayJetPt_array.size()-1));
    FillHist("N_of_Fake_Events_Differ_Electron_up", awayJetPt_array.at(aaa), weight+weight_err, awayJetPt_values, (awayJetPt_array.size()-1));
    FillHist("N_of_Fake_Events_Differ_Electron_down", awayJetPt_array.at(aaa), weight-weight_err, awayJetPt_values, (awayJetPt_array.size()-1));

    FillHist("N_Yield_Differ_Electron_"+GetSuffix(awayJetPt_array, aaa), 0., weight, 0., 1., 1);
    FillHist("N_Yield_Differ_Electron_"+GetSuffix(awayJetPt_array, aaa)+"_up", 0., weight+weight_err, 0., 1., 1);
    FillHist("N_Yield_Differ_Electron_"+GetSuffix(awayJetPt_array, aaa)+"_down", 0., weight-weight_err, 0., 1., 1);

  }


  return;

}// End of execute event loop
  


void HNSSSFMuMuE_FR_syst::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNSSSFMuMuE_FR_syst::BeginCycle() throw( LQError ){
  
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

HNSSSFMuMuE_FR_syst::~HNSSSFMuMuE_FR_syst() {
  
  Message("In HNSSSFMuMuE_FR_syst Destructor" , INFO);
  
}


void HNSSSFMuMuE_FR_syst::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNSSSFMuMuE_FR_syst::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNSSSFMuMuE_FR_systCore::MakeHistograms() to make new hists for your analysis
   **/

}


void HNSSSFMuMuE_FR_syst::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
 
  out_muons.clear();
  out_electrons.clear();
}


TString HNSSSFMuMuE_FR_syst::GetSuffix(std::vector<double> RelIso_array, int RelIso_bin, std::vector<double> dXYSig_array, int dXYSig_bin){
  TString suffix="";

  for(int i=0; i<5; i++){
    if(RelIso_bin == i){
      suffix+="RelIso_0p";
      suffix+=TString::Itoa((int)(RelIso_array.at(i)*10),10);
      suffix+="__";
      break;
    }
  }
  if(RelIso_bin == 5){
    suffix+="RelISo_1p0__";
  }

  for(int i=0; i<3; i++){
    if(dXYSig_bin == i){
      suffix+="dXYSig_";
      suffix+=TString::Itoa((int)(dXYSig_array.at(i)),10);
      break;
    }
  }

  return suffix;
}

TString HNSSSFMuMuE_FR_syst::GetSuffix(std::vector<double> awayJetPt_array, int awayJetPt_bin){

  TString suffix="";

  for(int i=0; i<4; i++){
    if(awayJetPt_bin == i){
      suffix+="awayJetPt_";
      suffix+=TString::Itoa((int)(awayJetPt_array.at(i)),10);
      break;
    }
  }

  return suffix;
}
