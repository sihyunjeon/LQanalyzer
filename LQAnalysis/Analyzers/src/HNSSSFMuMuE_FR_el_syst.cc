// $Id: HNSSSFMuMuE_FR_el_syst.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQHNSSSFMuMuE_FR_el_syst Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "HNSSSFMuMuE_FR_el_syst.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (HNSSSFMuMuE_FR_el_syst);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
HNSSSFMuMuE_FR_el_syst::HNSSSFMuMuE_FR_el_syst() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("HNSSSFMuMuE_FR_el_syst");
  
  Message("In HNSSSFMuMuE_FR_el_syst constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  //MakeCleverHistograms(sighist_mm,"DiMuon");


}


void HNSSSFMuMuE_FR_el_syst::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  Message("Making clever hists for Z ->ll test code", INFO);

  return;
}


void HNSSSFMuMuE_FR_el_syst::ExecuteEvents()throw( LQError ){

  // ========== No cut ====================
  //m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  //m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
  FillCutFlow("NoCut", weight);
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

  std::vector<TString> RelIso_string_el;
  std::vector<TString> awayJetPt_string;
  string RelIso_cstring_el[6] = {"0p45", "0p50", "0p55", "0p60", "0p65", "0p75"};
  string awayJetPt_cstring[4] = {"20", "30", "40", "60"};
  for(int i=0; i<6; i++) RelIso_string_el.push_back(RelIso_cstring_el[i]);
  for(int i=0; i<4; i++) awayJetPt_string.push_back(awayJetPt_cstring[i]);

  std::vector<double> RelIso_array_el;
  std::vector<double> awayJetPt_array;
  float RelIso_values_el[6] = {0.45, 0.50, 0.55, 0.60, 0.65, 0.75};
  for(int i=0; i<6; i++){
    RelIso_array_el.push_back(RelIso_values_el[i]);
  }

  float awayJetPt_values[4] = {20, 30, 40, 60};
  for(int i=0;i<4; i++){
    awayJetPt_array.push_back(awayJetPt_values[i]);
  }

  double weight_err = -999.;


  muonLooseColl.clear();
  for(int mu_it=0; mu_it<muonVLooseColl.size(); mu_it++){
    if((muonVLooseColl.at(mu_it).RelIso04() < 0.4)){
      muonLooseColl.push_back(muonVLooseColl.at(mu_it));
    }
  }
  if( muonLooseColl.size() == 2){

    for(int aaa=0; aaa<(RelIso_array_el.size()); aaa++){

      electronLooseColl.clear();
      for(int el_it=0; el_it<electronVLooseColl.size(); el_it++){
        if((electronVLooseColl.at(el_it).PFRelIso(0.3) < RelIso_array_el.at(aaa))){
          electronLooseColl.push_back(electronVLooseColl.at(el_it));
        }
      }
      if( !(electronLooseColl.size() == 1) ) continue;

      snu::KParticle lep[3];
      lep[0] = muonLooseColl.at(0);
      lep[1] = muonLooseColl.at(1);
      lep[2] = electronLooseColl.at(0);

      if( (lep[0].Charge() != lep[1].Charge()) ) continue;
      if( (lep[0].Charge() == lep[2].Charge()) ) continue;
      if( (lep[1].Charge() == lep[2].Charge()) ) continue;

      if( (lep[0].Pt() < 20) || (lep[1].Pt() < 10) || (lep[2].Pt() < 10) ) continue;

      for(int bbb=0; bbb<(awayJetPt_array.size()); bbb++){
        m_datadriven_bkg->GetFakeObj()->SetTrilepWP(4.0, 0.4);
        m_datadriven_bkg->GetFakeObj()->SetTrilepElWP(awayJetPt_array.at(bbb));
        weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", 1);
        weight_err = m_datadriven_bkg->Get_DataDrivenWeight(true, muonLooseColl, "MUON_HN_TRI_TIGHT", 2, electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", 1);

        TString suffix = "RelIsomu0p4_dXYSig4p0__RelIsoel"+RelIso_string_el.at(aaa)+"_awayJetPt"+awayJetPt_string.at(bbb);

        FillUpDownHist("ElectronSys__"+suffix, 0., weight, weight_err, 0., 1., 1);
      }
    }
    goto stop;
  }

  stop:
  return;

}// End of execute event loop
  


void HNSSSFMuMuE_FR_el_syst::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void HNSSSFMuMuE_FR_el_syst::BeginCycle() throw( LQError ){
  
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

HNSSSFMuMuE_FR_el_syst::~HNSSSFMuMuE_FR_el_syst() {
  
  Message("In HNSSSFMuMuE_FR_el_syst Destructor" , INFO);
  
}


void HNSSSFMuMuE_FR_el_syst::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void HNSSSFMuMuE_FR_el_syst::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this HNSSSFMuMuE_FR_el_systCore::MakeHistograms() to make new hists for your analysis
   **/

}


void HNSSSFMuMuE_FR_el_syst::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
 
  out_muons.clear();
  out_electrons.clear();
}


int HNSSSFMuMuE_FR_el_syst::GetPeriodIndex(void){
  if( isData ){
    if( k_sample_name.Contains("B") ||  k_sample_name.Contains("C") || k_sample_name.Contains("D") || k_sample_name.Contains("E") || k_sample_name.Contains("F") ){
      return 1;
    }
    else return 7;
  }
  else return 0;
}
