// $Id: DataValidation_EMu.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQDataValidation_EMu Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "DataValidation_EMu.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (DataValidation_EMu);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
DataValidation_EMu::DataValidation_EMu() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("DataValidation_EMu");
  
  Message("In DataValidation_EMu constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void DataValidation_EMu::InitialiseAnalysis() throw( LQError ) {
  
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


void DataValidation_EMu::ExecuteEvents()throw( LQError ){

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

   float pu_reweight=(1.0);
   if (!k_isdata) {   pu_reweight = eventbase->GetEvent().PileUpWeight();}
     
   TString emu_trig = "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v";

   // Now you should do an OR of 4 triggers 
 
   vector<TString> triggers;
   triggers.push_back(emu_trig);

   std::vector<snu::KElectron> electronTightColl = GetElectrons(true, true, "ELECTRON_HN_TIGHT");
   std::vector<snu::KElectron> electronLooseColl = GetElectrons(true, true, "ELECTRON_HN_FAKELOOSE");

   std::vector<snu::KJet> jets =   GetJets("JET_HN");
   int nbjet = NBJet(GetJets("JET_HN"));
   std::vector<snu::KMuon> muonTightColl =GetMuons("MUON_HN_TRI_TIGHT",true); 
   std::vector<snu::KMuon> muonLooseColl =GetMuons("MUON_HN_TRI_LOOSE",true);

   bool trig_pass= PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
   if(!trig_pass) return;

   float trigger_sf =  mcdata_correction->TriggerScaleFactor(electronLooseColl, muonLooseColl, "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
   CorrectMuonMomentum(muonLooseColl); /// CorrectMuonMomentum(muons);  will also work as Funcion in AnalyzerCore just calls mcdata_correction function

   float trigger_ps = WeightByTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", TargetLumi);   

   double ev_weight = weight;
   if(!isData){
//     weight = weight * trigger_sf;
     weight = weight * trigger_ps;
     weight = weight * pu_reweight;
   }

  if( !(muonLooseColl.size() == 1 && electronLooseColl.size() == 1) ) return;
  if( (muonLooseColl.at(0).Charge() ==  electronLooseColl.at(0).Charge()) ) return;

  snu::KMuon mu;
  snu::KElectron el;

  mu = muonLooseColl.at(0);
  el = electronLooseColl.at(0);

  if(mu.Pt() < 15 || el.Pt() < 25) return;

  double MET = eventbase->GetEvent().PFMET();
  double PI=3.141592;
  double PIarray[9] = {-PI, -PI*(3./4.), -PI*(2./4.), -PI*(1./4.), 0, PI*(1./4.), PI*(2./4.), PI*(3./4.), PI};

  if( !(muonTightColl.size() == 1 && electronTightColl.size() == 1) ) return;
  FillHist("0", 0, weight, 0., 1., 1);

  FillHist("TT_elecron_charge", el.Charge(), weight, -2., 3., 5);
  FillHist("TT_muon_charge", mu.Charge(), weight, -2., 3., 5);

  FillHist("TT_electron_pt", el.Pt(), weight, 0., 400., 400);
  FillHist("TT_muon_pt", mu.Pt(), weight, 0., 400., 400);

  FillHist("TT_met", MET, weight, 0., 400., 400);

  FillHist("TT_emu_mass", (el+mu).M(), weight, 0., 400., 400);

  for( int i=0; i<8; i++){
    if( (PIarray[i]< el.Phi()) && (el.Phi()<PIarray[i+1]) ){
      FillHist("TT_electron_phi", i, weight, 0., 8., 8);
    }
    if( (PIarray[i]< mu.Phi()) && (mu.Phi()<PIarray[i+1]) ){
      FillHist("TT_muon_phi", i, weight, 0., 8., 8);
    }
  }
 	    
  if(nbjet == 0){
    FillHist("1", 0, weight, 0., 1., 1);

    FillHist("nob_TT_elecron_charge", el.Charge(), weight, -2., 3., 5);
    FillHist("nob_TT_muon_charge", mu.Charge(), weight, -2., 3., 5);

    FillHist("nob_TT_electron_pt", el.Pt(), weight, 0., 400., 400);
    FillHist("nob_TT_muon_pt", mu.Pt(), weight, 0., 400., 400);

    FillHist("nob_TT_met", MET, weight, 0., 400., 400);

    FillHist("nob_TT_emu_mass", (el+mu).M(), weight, 0., 400., 400);

    for( int i=0; i<8; i++){
      if( (PIarray[i]< el.Phi()) && (el.Phi()<PIarray[i+1]) ){
        FillHist("nob_TT_electron_phi", i, weight, 0., 8., 8);
      }
      if( (PIarray[i]< mu.Phi()) && (mu.Phi()<PIarray[i+1]) ){
        FillHist("nob_TT_muon_phi", i, weight, 0., 8., 8);
      }
    }


    if( MET < 40 ){
      FillHist("2", 0, weight, 0., 1., 1);

      FillHist("met40_nob_TT_elecron_charge", el.Charge(), weight, -2., 3., 5);
      FillHist("met40_nob_TT_muon_charge", mu.Charge(), weight, -2., 3., 5);

      FillHist("met40_nob_TT_electron_pt", el.Pt(), weight, 0., 400., 400);
      FillHist("met40_nob_TT_muon_pt", mu.Pt(), weight, 0., 400., 400);

      FillHist("met40_nob_TT_met", MET, weight, 0., 400., 400);

      FillHist("met40_nob_TT_emu_mass", (el+mu).M(), weight, 0., 400., 400);

      for( int i=0; i<8; i++){
        if( (PIarray[i]< el.Phi()) && (el.Phi()<PIarray[i+1]) ){
          FillHist("met40_nob_TT_electron_phi", i, weight, 0., 8., 8);
        }
        if( (PIarray[i]< mu.Phi()) && (mu.Phi()<PIarray[i+1]) ){
          FillHist("met40_nob_TT_muon_phi", i, weight, 0., 8., 8);
        }
      }



    }
  }
  else{
    FillHist("3", 0, weight, 0., 1., 1);

    FillHist("b_TT_elecron_charge", el.Charge(), weight, -2., 3., 5);
    FillHist("b_TT_muon_charge", mu.Charge(), weight, -2., 3., 5);

    FillHist("b_TT_electron_pt", el.Pt(), weight, 0., 400., 400);
    FillHist("b_TT_muon_pt", mu.Pt(), weight, 0., 400., 400);

    FillHist("b_TT_met", MET, weight, 0., 400., 400);

    FillHist("b_TT_emu_mass", (el+mu).M(), weight, 0., 400., 400);

    for( int i=0; i<6; i++){
      if( (PIarray[i]< el.Phi()) && (el.Phi()<PIarray[i+1]) ){
        FillHist("b_TT_electron_phi", i, weight, 0., 6., 6);
      }
      if( (PIarray[i]< mu.Phi()) && (mu.Phi()<PIarray[i+1]) ){
        FillHist("b_TT_muon_phi", i, weight, 0., 6., 6);
      }
    }

  }


   

   return;
}// End of execute event loop
  


void DataValidation_EMu::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void DataValidation_EMu::BeginCycle() throw( LQError ){
  
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

DataValidation_EMu::~DataValidation_EMu() {
  
  Message("In DataValidation_EMu Destructor" , INFO);
  
}


void DataValidation_EMu::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void DataValidation_EMu::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this DataValidation_EMuCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void DataValidation_EMu::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



