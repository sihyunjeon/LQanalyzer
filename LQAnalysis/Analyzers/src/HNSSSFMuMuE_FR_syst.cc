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

  // ========== No cut ====================
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  // ================================================================================


  // ========== MET filter cut ====================
  if(!PassMETFilter()) return;     /// Initial event cuts : 
  // ================================================================================


  // ========== Primary vertex cut ====================
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex 
  // ================================================================================


  // ========== Trigger cut ====================
  std::vector<TString> triggerlist;
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
  triggerlist.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");

  if(!PassTriggerOR(triggerlist)) return;
  // ================================================================================


  // ========== Get Objects (muon, electron, jet) ====================
  std::vector<snu::KMuon> muonVLooseColl = GetMuons("MUON_HN_TRI_VLOOSE",false);
  std::vector<snu::KMuon> muonLooseColl;
  muonLooseColl.clear();

  std::vector<snu::KElectron> electronVLooseColl = GetElectrons(false,false,"ELECTRON_HN_FAKEVLOOSE");
  std::vector<snu::KElectron> electronLooseColl;
  electronLooseColl.clear();

  std::vector<snu::KJet> jetLooseColl = GetJets("JET_HN", 30., 2.4);
  int period_index = 0;
  period_index = GetPeriodIndex();
  int nbjet = NBJet(jetLooseColl, snu::KJet::CSVv2, snu::KJet::Medium, period_index);
  if( nbjet > 0 ) return;
  // ================================================================================


  /*####################################################################################################
  ##		        Analysis Code 								      ##
  ##				For SameSign MuMuE Channel Analysis				      ##
  ####################################################################################################*/

  std::vector<float> RelIso_array_mu;
  std::vector<float> dXYSig_array;
//  std::vector<float> RelIso_array_el;
  std::vector<float> awayJetPt_array;

  std::vector<TString> RelIso_string_mu;
//  std::vector<TString> RelIso_string_el;
  std::vector<TString> dXYSig_string;
  std::vector<TString> awayJetPt_string;

  string RelIso_cstring_mu[6] = {"0p20", "0p30", "0p40", "0p60", "0p80", "1p00"};
//  string RelIso_cstring_el[6] = {"0p45", "0p50", "0p55", "0p60", "0p65", "0p75"};
  string dXYSig_cstring[3] = {"3p0", "4p0", "5p0"};
  string awayJetPt_cstring[4] = {"20", "30", "40", "60"};
  for(int i=0; i<6; i++) RelIso_string_mu.push_back(RelIso_cstring_mu[i]);
//  for(int i=0; i<6; i++) RelIso_string_el.push_back(RelIso_cstring_el[i]);
  for(int i=0; i<3; i++) dXYSig_string.push_back(dXYSig_cstring[i]);
  for(int i=0; i<4; i++) awayJetPt_string.push_back(awayJetPt_cstring[i]);

  float RelIso_values_mu[7] = {0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2};
  for(int i=0; i<6; i++){
    RelIso_array_mu.push_back(RelIso_values_mu[i]);
  }

  float dXYSig_values[4] = {3, 4, 5, 6};
  for(int i=0; i<3; i++){
    dXYSig_array.push_back(dXYSig_values[i]);
  }

/*  float RelIso_values_el[6] = {0.45, 0.50, 0.55, 0.60, 0.65, 0.75};
  for(int i=0; i<6; i++){
    RelIso_array_el.push_back(RelIso_values_el[i]);
  }*/

  float awayJetPt_values[5] = {20, 30, 40, 60, 100};
  for(int i=0;i<4; i++){
    awayJetPt_array.push_back(awayJetPt_values[i]);
  }

  float weight_err = -999.;

  electronLooseColl.clear();
  for(int el_it=0; el_it<electronVLooseColl.size(); el_it++){
    if((electronVLooseColl.at(el_it).PFRelIso(0.3) < 0.5)){
      electronLooseColl.push_back(electronVLooseColl.at(el_it));
    }
  }

  if( (electronLooseColl.size() == 1) ){

    for(int aaa=0; aaa<(RelIso_array_mu.size()); aaa++){

      muonLooseColl.clear();
      for(int mu_it=0; mu_it<muonVLooseColl.size(); mu_it++){
        if((muonVLooseColl.at(mu_it).RelIso04() < RelIso_array_mu.at(aaa))){
          muonLooseColl.push_back(muonVLooseColl.at(mu_it));
        }
      }
      if( !(muonLooseColl.size() == 2) ) continue;

      snu::KParticle lep[3];
      lep[0] = muonLooseColl.at(0);
      lep[1] = muonLooseColl.at(1);
      lep[2] = electronLooseColl.at(0);
  
      if( (lep[0].Charge() != lep[1].Charge()) ) continue;
      if( (lep[0].Charge() == lep[2].Charge()) ) continue;
      if( (lep[1].Charge() == lep[2].Charge()) ) continue;
 
      if( (lep[0].Pt() < 20) || (lep[1].Pt() < 10) || (lep[2].Pt() < 10) ) continue;

      for(int bbb=0; bbb<(dXYSig_array.size()); bbb++){

        for(int ccc=0; ccc<(awayJetPt_array.size()); ccc++){

          m_datadriven_bkg->GetFakeObj()->SetTrilepWP(dXYSig_array.at(bbb), RelIso_array_mu.at(aaa));
          m_datadriven_bkg->GetFakeObj()->SetTrilepElWP(awayJetPt_array.at(ccc));
          weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", 1);
          weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", 1);

      //    TString suffix = "RelIso_mu_"+RelIso_string_mu.at(aaa)+"__RelIso_el_0p5__dXYSig_"+dXYSig_string.at(bbb)+"__awayJetPt_"+awayJetPt_string.at(ccc);

          FillHist("N_MuFake_Syst_awayJetPt"+awayJetPt_string.at(ccc), dXYSig_array.at(bbb), RelIso_array_mu.at(aaa), weight, dXYSig_values, dXYSig_array.size(), RelIso_values_mu, RelIso_array_mu.size());
          FillHist("N_MuFake_Syst_awayJetPt"+awayJetPt_string.at(ccc)+"_up", dXYSig_array.at(bbb), RelIso_array_mu.at(aaa), weight+weight_err, dXYSig_values, dXYSig_array.size(), RelIso_values_mu, RelIso_array_mu.size());
  	  FillHist("N_MuFake_Syst_awayJetPt"+awayJetPt_string.at(ccc)+"_down", dXYSig_array.at(bbb), RelIso_array_mu.at(aaa), weight-weight_err, dXYSig_values, dXYSig_array.size(), RelIso_values_mu, RelIso_array_mu.size());


        }

      }
    }
  }

  muonLooseColl.clear();
  for(int mu_it=0; mu_it<muonVLooseColl.size(); mu_it++){
    if((muonVLooseColl.at(mu_it).RelIso04() < 0.4)){
      muonLooseColl.push_back(muonVLooseColl.at(mu_it));
    }
  }
  if( (muonLooseColl.size() == 2) ){
    electronLooseColl.clear();
    for(int el_it=0; el_it<electronVLooseColl.size(); el_it++){
      if((electronVLooseColl.at(el_it).PFRelIso(0.3) < 0.5)){
        electronLooseColl.push_back(electronVLooseColl.at(el_it));
      }
    }
    if( !(electronLooseColl.size() == 1) ) return;

    snu::KParticle lep[3];
    lep[0] = muonLooseColl.at(0);
    lep[1] = muonLooseColl.at(1);
    lep[2] = electronLooseColl.at(0);

    if( (lep[0].Charge() != lep[1].Charge()) ) return;
    if( (lep[0].Charge() == lep[2].Charge()) ) return;
    if( (lep[1].Charge() == lep[2].Charge()) ) return;

    if( (lep[0].Pt() < 20) || (lep[1].Pt() < 10) || (lep[2].Pt() < 10) ) return;

    for(int aaa=0; aaa<(awayJetPt_array.size()); aaa++){

      m_datadriven_bkg->GetFakeObj()->SetTrilepWP(4.0, 0.4);
      m_datadriven_bkg->GetFakeObj()->SetTrilepElWP(awayJetPt_array.at(aaa));
      weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", 1);
      weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", 1);

      FillHist("N_ElFake_Syst", awayJetPt_array.at(aaa), weight, awayJetPt_values, awayJetPt_array.size());
      FillHist("N_ElFake_Syst_up", awayJetPt_array.at(aaa), weight+weight_err, awayJetPt_values, awayJetPt_array.size());
      FillHist("N_ElFake_Syst_down", awayJetPt_array.at(aaa), weight-weight_err, awayJetPt_values, awayJetPt_array.size());

    }

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


int HNSSSFMuMuE_FR_syst::GetPeriodIndex(void){
  if( isData ){
    if( k_sample_name.Contains("B") ||  k_sample_name.Contains("C") || k_sample_name.Contains("D") || k_sample_name.Contains("E") || k_sample_name.Contains("F") ){
      return 1;
    }
    else return 7;
  }
  else return 0;
}
