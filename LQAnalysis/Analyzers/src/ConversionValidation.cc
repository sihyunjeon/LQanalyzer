// $Id: ConversionValidation.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQConversionValidation Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "ConversionValidation.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (ConversionValidation);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
ConversionValidation::ConversionValidation() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("ConversionValidation");
  
  Message("In ConversionValidation constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void ConversionValidation::InitialiseAnalysis() throw( LQError ) {
  
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


void ConversionValidation::ExecuteEvents()throw( LQError ){

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

  bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if(!pass_trig) return;

  // ========== Get Objects (muon, electron, jet) ====================
  std::vector<snu::KMuon> muonVetoColl = GetMuons("MUON_HN_VETO", true);
  std::vector<snu::KMuon> muons;  muons.clear();
  std::vector<snu::KElectron> electronVetoColl = GetElectrons(false, false, "ELECTRON_HN_VETO");
  std::vector<snu::KElectron> electrons;  electrons.clear();

  if(muonVetoColl.size()!=0) return;

  for(unsigned int i=0; i<electronVetoColl.size(); i++){
    if(PassID(electronVetoColl.at(i), "ELECTRON_HN_TIGHTv4")){
      if(fabs(electronVetoColl.at(i).SCEta()) < 1.4442 || fabs(electronVetoColl.at(i).SCEta()) > 1.5560){
        electrons.push_back(electronVetoColl.at(i));
        if(electronVetoColl.at(i).Pt() <25) return;
      }
    }
  }

  if(electronVetoColl.size() != electrons.size()) return;
  if(!(muons.size() == 0 && electrons.size() == 2)) return;
  
  double METPt = eventbase->GetEvent().MET();
  double METPhi = eventbase->GetEvent().METPhi();
  snu::KParticle MET;
  MET.SetPxPyPzE(METPt*TMath::Cos(METPhi), METPt*TMath::Sin(METPhi), 0, METPt);
  // ================================================================================

  // Get event weight
  double this_weight = 1.;
  if(!isData){
    this_weight = weight;
    this_weight *= MCweight;  
    this_weight *= mcdata_correction->ElectronScaleFactor("ELECTRON_HN_TIGHTv4", electrons);
    this_weight *= mcdata_correction->ElectronRecoScaleFactor(electrons);
    this_weight *= mcdata_correction->CatPileupWeight(eventbase->GetEvent(),0);
    double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "ELECTRON_HN_TIGHTv4", muons, "MUON_HN_TIGHT", 0, 0, 0);
    double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, "ELECTRON_HN_TIGHTv4", muons, "MUON_HN_TIGHT", 0, 1, 0);
    double trigger_sf = trigger_eff_Data/trigger_eff_MC;
    this_weight *= trigger_sf;
    this_weight *= GetKFactor();
    this_weight *= 35863.307;
  }
  // ================================================================================

  bool is_SS = false;
  if(electrons.at(0).Charge() == electrons.at(1).Charge()) is_SS = true;
  bool is_TYPE40 = false;
  if((electrons.at(0).GetType() == 40) || (electrons.at(1).GetType() == 40)) is_TYPE40 = true;
  bool is_CONV_JOHN[2] = {false,false};
  bool is_CONV_JIHWAN[2] = {false,false};
  bool is_CONV_OR[2] = {false,false};
  bool is_Barrel[2] = {false,false};
  for(int i=0; i<2; i++){
    if((electrons.at(i).GetType() == 40)) is_CONV_JOHN[i] = true;
    if(JHsConv(electrons.at(i))) is_CONV_JIHWAN[i] = true;
    if(IsExternalConversion(electrons.at(i))) is_CONV_OR[i] = true;

    if(fabs(electrons.at(i).SCEta()) < 1.4442) is_Barrel[i] = true;
    else if(fabs(electrons.at(i).SCEta()) > 1.5560) is_Barrel[i] = false;
    else return;
  }

  snu::KParticle Z_candidate = electrons.at(0) + electrons.at(1);
  double Z_Range = 30., Min_Pt = 25., N_Bins = 40.;
  bool PassCuts = ((fabs((electrons.at(0)+electrons.at(1)).M()-91.1876) < Z_Range) && (electrons.at(1).Pt() > Min_Pt));

  TString s_SIGN= "NULL";
  if(is_SS) s_SIGN = "SS_";  else s_SIGN = "OS_";

  TString s_REGION = "NULL";
  if(is_Barrel[0] && is_Barrel[1]) s_REGION = "BB_";
  else if(!is_Barrel[0] && !is_Barrel[1]) s_REGION = "EE_";
  else s_REGION = "BE_";

  if(PassCuts){
    DrawHistColl(s_SIGN+s_REGION+"DIELECTRON_WOCONV_", electrons.at(0), this_weight);
    DrawHistColl(s_SIGN+s_REGION+"DIELECTRON_WOCONV_", electrons.at(1), this_weight);
    DrawHistColl(s_SIGN+s_REGION+"DIELECTRON_WOCONV_1st_", electrons.at(0), this_weight);
    DrawHistColl(s_SIGN+s_REGION+"DIELECTRON_WOCONV_2nd_", electrons.at(1), this_weight);
    DrawHistColl("ALL_"+s_REGION+"DIELECTRON_WOCONV_", electrons.at(0), this_weight);
    DrawHistColl("ALL_"+s_REGION+"DIELECTRON_WOCONV_", electrons.at(1), this_weight);
    DrawHistColl("ALL_"+s_REGION+"DIELECTRON_WOCONV_1st_", electrons.at(0), this_weight);
    DrawHistColl("ALL_"+s_REGION+"DIELECTRON_WOCONV_2nd_", electrons.at(1), this_weight);

    FillHist(s_SIGN+s_REGION+"DIELECTRON_WOCONV_Zmass", Z_candidate.M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range,N_Bins);
    FillHist(s_SIGN+s_REGION+"DIELECTRON_WOCONV_Zpt", Z_candidate.Pt(), this_weight, 0., 100., 100);
    FillHist(s_SIGN+s_REGION+"DIELECTRON_WOCONV_Nevents", 0., this_weight, 0., 1., 1);

    FillHist("ALL_"+s_REGION+"DIELECTRON_WOCONV_Zmass", Z_candidate.M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range,N_Bins);
    FillHist("ALL_"+s_REGION+"DIELECTRON_WOCONV_Zpt", Z_candidate.Pt(), this_weight, 0., 100., 100);
    FillHist("ALL_"+s_REGION+"DIELECTRON_WOCONV_Nevents", 0., this_weight, 0., 1., 1);


    if(is_CONV_JOHN[0]){
      DrawHistColl(s_SIGN+s_REGION+"DIELECTRON_WCONV_JOHN_", electrons.at(0), this_weight);
      DrawHistColl(s_SIGN+s_REGION+"DIELECTRON_WCONV_JOHN_1st_", electrons.at(0), this_weight);
    }

    if(is_CONV_JOHN[1]){
      DrawHistColl(s_SIGN+s_REGION+"DIELECTRON_WCONV_JOHN_", electrons.at(1), this_weight);
      DrawHistColl(s_SIGN+s_REGION+"DIELECTRON_WCONV_JOHN_2nd_", electrons.at(1), this_weight);
    }

    if(is_CONV_JOHN[0] || is_CONV_JOHN[1]){
      FillHist(s_SIGN+s_REGION+"DIELECTRON_WCONV_JOHN_Zmass_OR", Z_candidate.M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range,N_Bins);
      FillHist(s_SIGN+s_REGION+"DIELECTRON_WCONV_JOHN_Zpt_OR", Z_candidate.Pt(), this_weight, 0., 100., 100);
      FillHist(s_SIGN+s_REGION+"DIELECTRON_WCONV_Nevents", 0., this_weight, 0., 1., 1);

      FillHist("ALL_"+s_REGION+"DIELECTRON_WCONV_JOHN_Zmass_OR", Z_candidate.M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range,N_Bins);
      FillHist("ALL_"+s_REGION+"DIELECTRON_WCONV_JOHN_Zpt_OR", Z_candidate.Pt(), this_weight, 0., 100., 100);
      FillHist("ALL_"+s_REGION+"DIELECTRON_WCONV_Nevents", 0., this_weight, 0., 1., 1);

    }

    if(is_CONV_JIHWAN[0]){
      DrawHistColl(s_SIGN+"DIELECTRON_WCONV_JIHWAN_", electrons.at(0), this_weight);
      DrawHistColl(s_SIGN+"DIELECTRON_WCONV_JIHWAN_1st_", electrons.at(0), this_weight);
    } 

    if(is_CONV_JIHWAN[1]){
      DrawHistColl(s_SIGN+"DIELECTRON_WCONV_JIHWAN_", electrons.at(1), this_weight);
      DrawHistColl(s_SIGN+"DIELECTRON_WCONV_JIHWAN_2nd_", electrons.at(1), this_weight);
    }

    if(is_CONV_JIHWAN[0] || is_CONV_JIHWAN[1]){
      FillHist(s_SIGN+"DIELECTRON_WCONV_JIHWAN_Zmass_OR", Z_candidate.M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range,N_Bins);
      FillHist(s_SIGN+"DIELECTRON_WCONV_JIHWAN_Zpt_OR", Z_candidate.Pt(), this_weight, 0., 100., 100);
    }

    if(is_CONV_OR[0]){
      DrawHistColl(s_SIGN+"DIELECTRON_WCONV_OR_", electrons.at(0), this_weight);
      DrawHistColl(s_SIGN+"DIELECTRON_WCONV_OR_1st_", electrons.at(0), this_weight);
    } 

    if(is_CONV_OR[1]){
      DrawHistColl(s_SIGN+"DIELECTRON_WCONV_OR_", electrons.at(1), this_weight);
      DrawHistColl(s_SIGN+"DIELECTRON_WCONV_OR_2nd_", electrons.at(1), this_weight);
    }

    if(is_CONV_OR[0] || is_CONV_OR[1]){
      FillHist(s_SIGN+"DIELECTRON_WCONV_OR_Zmass_OR", Z_candidate.M(), this_weight, 91.1876-Z_Range, 91.1876+Z_Range,N_Bins);
      FillHist(s_SIGN+"DIELECTRON_WCONV_OR_Zpt_OR", Z_candidate.Pt(), this_weight, 0., 100., 100);
      if(s_SIGN == "SS_"){
        if(Z_candidate.M() >85){
          cout<<electrons.at(0).Pt()<<"\t"<<electrons.at(0).Eta()<<"\t"<<electrons.at(0).Phi()<<"\t"<<electrons.at(0).Charge()<<"\t"<<electrons.at(0).GetType()<<endl;
	  cout<<electrons.at(1).Pt()<<"\t"<<electrons.at(1).Eta()<<"\t"<<electrons.at(1).Phi()<<"\t"<<electrons.at(1).Charge()<<"\t"<<electrons.at(0).GetType()<<endl;
          TruthPrintOut();
        }
      }
    }



  }

  return;
}// End of execute event loop

void ConversionValidation::DrawHistColl( TString this_string, snu::KElectron electron, double this_weight ){

  KLepton lepton;
  lepton = electron;

  TString this_obj_string = this_string+"lepton_";
  FillHist(this_obj_string+"pt", lepton.Pt(), this_weight, 0., 400., 400);
  FillHist(this_obj_string+"eta", lepton.Eta(), this_weight, -5., 5., 200);
  FillHist(this_obj_string+"dxy", lepton.dXY(), this_weight, -1., 1., 200);
  FillHist(this_obj_string+"dz", lepton.dZ(), this_weight, -1., 1., 200);
  FillHist(this_obj_string+"reliso", lepton.RelIso(), this_weight, 0., 1., 200);

  FillHist(this_obj_string+"nevents", 0., this_weight, 0., 1., 1);

}

void ConversionValidation::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void ConversionValidation::BeginCycle() throw( LQError ){
  
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

ConversionValidation::~ConversionValidation() {
  
  Message("In ConversionValidation Destructor" , INFO);
  
}


void ConversionValidation::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void ConversionValidation::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this ConversionValidationCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void ConversionValidation::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


bool ConversionValidation::JHsConv( snu::KElectron lepton ){

  std::vector<snu::KTruth> truthColltemp= eventbase->GetTruth();
//  cout<<"!!!"<<GetLeptonType(lepton, truthColltemp )<<endl;
  FillHist("LEPTON_TYPE_JIHWAN", GetLeptonType(lepton, truthColltemp ), 1., -6., 6., 12);
  FillHist("LEPTON_TYPE_JOHN", lepton.GetType(), 1., 0., 41., 41);

  if(GetLeptonType(lepton, truthColltemp ) == -5 || GetLeptonType(lepton, truthColltemp ) == -6 || GetLeptonType(lepton, truthColltemp ) == -6){
//    cout <<"!"<<endl;
    return true;
  }
  else return false;

}

