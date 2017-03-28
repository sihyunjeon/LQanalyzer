// $Id: CFRateCalculator.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQCFRateCalculator Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "CFRateCalculator.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (CFRateCalculator);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
CFRateCalculator::CFRateCalculator() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("CFRateCalculator");
  
  Message("In CFRateCalculator constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void CFRateCalculator::InitialiseAnalysis() throw( LQError ) {
  
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


void CFRateCalculator::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);


  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", weight);

  bool trig_pass = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if(!trig_pass) return;

  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                                              
  std::vector<snu::KElectron> electronTestColl = GetElectrons(true, false, "ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KElectron> electronLooseColl =  GetElectrons(true,false,"ELECTRON_HN_LOWDXY_FAKELOOSE");
  std::vector<snu::KElectron> electronTightColl =  GetElectrons(true,false,"ELECTRON_HN_LOWDXY_TIGHT");
  FillHist("n_testColl", electronTestColl.size(), 1., 0., 8., 8);
  FillHist("n_lowdxyColl", electronLooseColl.size(), 1., 0., 8., 8);

  //================== weight
  float weight_trigger = WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);

  float pileup_reweight=(1.0);
  if (!k_isdata) {   pileup_reweight = mcdata_correction->PileupWeightByPeriod(eventbase->GetEvent());}
     
  float electron_idsf = mcdata_correction->ElectronScaleFactor("ELECTRON_HN_LOWDXY_TIGHT", electronLooseColl);

  float electron_recosf = mcdata_correction->ElectronRecoScaleFactor(electronLooseColl);
  //====================================================================================================================

  if( electronLooseColl.size() != 2) return;
  if( electronTightColl.size() != 2) return;

  snu::KElectron lep[2];
  lep[0] = electronTightColl.at(0);
  lep[1] = electronTightColl.at(1);
  if( (lep[0].Pt() < 25) || (lep[1].Pt() < 25) ) return;

  double invPt[2];
  invPt[0] = (1./(lep[0].Pt()));
  invPt[1] = (1./(lep[1].Pt()));

  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();

  if(METPt > 30) return;

  snu::KParticle Z_candidate;
  Z_candidate = (lep[0] + lep[1]);

  if(!isData){
    weight *= weight_trigger;
    weight *= pileup_reweight;
    weight *= electron_idsf;
    weight *= electron_recosf;
  }

  FillHist("Z_candidate_mass", Z_candidate.M(), weight, 0., 200., 200);

  bool is_charge_flip = false;
  if( lep[0].MCIsCF() || lep[1].MCIsCF() ) is_charge_flip = true;
  if( isData ) is_charge_flip = true;

  bool is_region[3][2] = {{false,},};
  if( (fabs(lep[0].Eta()) < 0.9) ) 					is_region[0][0] = true;
  else if( (fabs(lep[0].Eta()) < 1.4442) )				is_region[1][0] = true;
  else if( (fabs(lep[0].Eta()) > 1.556) && (fabs(lep[0].Eta()) < 2.5) ) is_region[2][0] = true;

  if( (fabs(lep[1].Eta()) < 0.9) ) 					is_region[0][1] = true;
  else if( (fabs(lep[1].Eta()) < 1.4442) ) 				is_region[1][1] = true;
  else if( (fabs(lep[1].Eta()) > 1.556) && (fabs(lep[1].Eta()) < 2.5) ) is_region[2][1] = true;

  for(int i=0; i<2; i++){
    for(int j=0; j<3; j++){
      if( is_region[j][i] ){
	DrawHistograms( lep[i], j, weight, false, false );
	if( is_charge_flip ) DrawHistograms( lep[i], j, weight, true, false );
      }	
    }
  }
  

  if( (fabs(Z_candidate.M() - 90) < 10) ){
    for(int i=0; i<2; i++){
      for(int j=0; j<3; j++){
        if( is_region[j][i] ){
          DrawHistograms( lep[i], j, weight, false, true );
          if( is_charge_flip ) DrawHistograms( lep[i], j, weight, true, true );
        }
      }
    } 
  } 

  
 
  return;
}// End of execute event loop
  


void CFRateCalculator::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void CFRateCalculator::BeginCycle() throw( LQError ){
  
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

CFRateCalculator::~CFRateCalculator() {
  
  Message("In CFRateCalculator Destructor" , INFO);
  
}


void CFRateCalculator::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void CFRateCalculator::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this CFRateCalculatorCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void CFRateCalculator::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


void CFRateCalculator::DrawHistograms( snu::KElectron lep, int eta_region, double weight, bool is_CF, bool is_Z ){

  if( eta_region == 0 ){
    FillHist("eta_invPt_region1", (1./(lep.Pt())), weight, 0., 0.05, 10);
    FillHist("eta_Pt_region1", lep.Pt(), weight, 0., 200, 20);
  }
  else if( eta_region == 1 ){
    FillHist("eta_invPt_region2", (1./(lep.Pt())), weight, 0., 0.05, 10);
    FillHist("eta_Pt_region2", lep.Pt(), weight, 0., 200, 20);
  }
  else if( eta_region == 2 ){
    FillHist("eta_invPt_region3", (1./(lep.Pt())), weight, 0., 0.05, 10);
    FillHist("eta_Pt_region3", lep.Pt(), weight, 0., 200, 20);
  }

  if(is_CF){
    if( eta_region == 0 ){
      FillHist("eta_invPt_region1_CF", (1./(lep.Pt())), weight, 0., 0.05, 10);
      FillHist("eta_Pt_region1_CF", lep.Pt(), weight, 0., 200, 20);
    }
    else if( eta_region == 1 ){
      FillHist("eta_invPt_region2_CF", (1./(lep.Pt())), weight, 0., 0.05, 10);
      FillHist("eta_Pt_region2_CF", lep.Pt(), weight, 0., 200, 20);
    }
    else if( eta_region == 2 ){
      FillHist("eta_invPt_region3_CF", (1./(lep.Pt())), weight, 0., 0.05, 10);
      FillHist("eta_Pt_region3_CF", lep.Pt(), weight, 0., 200, 20);
    }
  }

  if( is_Z ){
    if( eta_region == 0 ){
      FillHist("eta_invPt_region1_Z", (1./(lep.Pt())), weight, 0., 0.05, 10);
      FillHist("eta_Pt_region1_Z", lep.Pt(), weight, 0., 200, 20);
    }
    else if( eta_region == 1 ){
      FillHist("eta_invPt_region2_Z", (1./(lep.Pt())), weight, 0., 0.05, 10);
      FillHist("eta_Pt_region2_Z", lep.Pt(), weight, 0., 200, 20);
    }
    else if( eta_region == 2 ){
      FillHist("eta_invPt_region3_Z", (1./(lep.Pt())), weight, 0., 0.05, 10);
      FillHist("eta_Pt_region3_Z", lep.Pt(), weight, 0., 200, 20);
    }

    if(is_CF){
      if( eta_region == 0 ){
        FillHist("eta_invPt_region1_Z_CF", (1./(lep.Pt())), weight, 0., 0.05, 10);
        FillHist("eta_Pt_region1_Z_CF", lep.Pt(), weight, 0., 200, 20);
      }
      else if( eta_region == 1 ){
        FillHist("eta_invPt_region2_Z_CF", (1./(lep.Pt())), weight, 0., 0.05, 10);
        FillHist("eta_Pt_region2_Z_CF", lep.Pt(), weight, 0., 200, 20);
      }
      else if( eta_region == 2 ){
        FillHist("eta_invPt_region3_Z_CF", (1./(lep.Pt())), weight, 0., 0.05, 10);
        FillHist("eta_Pt_region3_Z_CF", lep.Pt(), weight, 0., 200, 20);
      }
    }
  }


}


