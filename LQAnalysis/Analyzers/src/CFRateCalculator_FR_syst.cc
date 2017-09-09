// $Id: CFRateCalculator_FR_syst.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQCFRateCalculator_FR_syst Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "CFRateCalculator_FR_syst.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (CFRateCalculator_FR_syst);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
CFRateCalculator_FR_syst::CFRateCalculator_FR_syst() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("CFRateCalculator_FR_syst");
  
  Message("In CFRateCalculator_FR_syst constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();


}


void CFRateCalculator_FR_syst::InitialiseAnalysis() throw( LQError ) {
  
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


void CFRateCalculator_FR_syst::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  if(!PassMETFilter()) return;     /// Initial event cuts :
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex

  CFvalidation();
  if(isData)  return;
 
  return;
}// End of execute event loop
  


void CFRateCalculator_FR_syst::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void CFRateCalculator_FR_syst::BeginCycle() throw( LQError ){
  
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

CFRateCalculator_FR_syst::~CFRateCalculator_FR_syst() {
  
  Message("In CFRateCalculator_FR_syst Destructor" , INFO);
  
}


void CFRateCalculator_FR_syst::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void CFRateCalculator_FR_syst::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this CFRateCalculator_FR_systCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void CFRateCalculator_FR_syst::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}


void CFRateCalculator_FR_syst::CFvalidation(void){

  bool pass_trig = PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  if( !pass_trig ) return;

  std::vector<snu::KMuon> muonLooseColl = GetMuons("MUON_HN_VETO", false);
  if( muonLooseColl.size() != 0) return;

  double Z_mass = 91.1876;
  double this_reliso_el = 0.6;

  TString el_ID = "ELECTRON_HN_TIGHT";
  TString IDsuffix = "_HNTIGHT";
  std::vector<snu::KElectron> electronVLooseColl = GetElectrons(false, false, "ELECTRON_HN_FAKELOOSE");
  std::vector<snu::KJet> jetVLooseColl = GetJets("JET_HN");

  snu::KEvent event = eventbase->GetEvent();

  int N_sys = (2*2+1);
  for(int it_sys = 0; it_sys<N_sys; it_sys++){

    double METPt = event.MET();
    double METPhi = event.METPhi();

    TString this_syst;

    if(it_sys==0){
      this_syst = "ElectronE_Up";//
      METPt = event.PFMETShifted(snu::KEvent::ElectronEn, snu::KEvent::up);
    }
    else if(it_sys==1){
      this_syst = "ElectronE_Down";//
      METPt = event.PFMETShifted(snu::KEvent::ElectronEn, snu::KEvent::down);
    }
    else if(it_sys==2){//
      this_syst = "Unclustered_Up";
      METPt = event.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::up);
    }
    else if(it_sys==3){//
      this_syst = "Unclustered_Down";
      METPt = event.PFMETShifted(snu::KEvent::Unclustered, snu::KEvent::down);
    }
    else if(it_sys==4){//
      this_syst = "Central";
    }
    else{
      Message("it_sys out of range!" , INFO);
      return;
    }

    std::vector<snu::KElectron> electronLooseColl, electronTightColl;
    if(this_syst == "ElectronE_Up"){
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedUp(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        double new_reliso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedUp();
        if( this_electron.Pt() >= 25. && new_reliso < this_reliso_el ) electronLooseColl.push_back( this_electron );
      }
    }
    else if(this_syst == "ElectronE_Down"){
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        this_electron.SetPtEtaPhiM( this_electron.Pt()*this_electron.PtShiftedDown(), this_electron.Eta(), this_electron.Phi(), this_electron.M() );
        double new_reliso = this_electron.PFRelIso(0.3)/this_electron.PtShiftedDown();
        if( this_electron.Pt() >= 25. && new_reliso < this_reliso_el ) electronLooseColl.push_back( this_electron );
      }
    }
    else{
      for(unsigned int i=0; i<electronVLooseColl.size(); i++){
        snu::KElectron this_electron = electronVLooseColl.at(i);
        if( this_electron.Pt() >= 25. && this_electron.PFRelIso(0.3) < this_reliso_el ) electronLooseColl.push_back( this_electron );
      }
    }
    for(unsigned int i=0; i<electronLooseColl.size(); i++){
      if(PassID(electronLooseColl.at(i), "ELECTRON_HN_TIGHT")) electronTightColl.push_back( electronLooseColl.at(i) );
    }

    if(METPt > 30) continue;

    TString CFsample="";

    if(electronLooseColl.size() != 2) continue;
    if(electronTightColl.size() == 2) continue;

    FillHist("[CHECK]n_of_muons"+IDsuffix+CFsample, muonLooseColl.size(), 1., 0., 5., 5);//check no muons
    FillHist("[CHECK]n_of_electrons"+IDsuffix+CFsample, electronTightColl.size(), 1., 0., 5., 5);//check two electrons

    // define leptons and give Pt, MET cuts
    snu::KParticle lep[2];
    lep[0] = electronLooseColl.at(0);
    lep[1] = electronLooseColl.at(1);

    bool is_SS = true;
    if( (lep[0].Charge() != lep[1].Charge()) ){ continue; }

    FillHist("[CHECK]MET"+IDsuffix+CFsample+"_"+this_syst, METPt, 1., 0., 100., 100);

    TString region = "";
    if( (fabs(lep[0].Eta()) < 0.9) )                                      region += "iB";
    else if( (fabs(lep[0].Eta()) < 1.4442) )				    region += "oB";
    else if( (fabs(lep[0].Eta()) > 1.556) && (fabs(lep[0].Eta()) < 2.5) ) region += "E";
    else continue;

    if( (fabs(lep[1].Eta()) < 0.9) )                                      region += "iB";
    else if( (fabs(lep[1].Eta()) < 1.4442) )                              region += "oB";
    else if( (fabs(lep[1].Eta()) > 1.556) && (fabs(lep[1].Eta()) < 2.5) ) region += "E";
    else continue;

    if((region == "iBE") || (region == "EiB") || (region == "oBE") || (region == "EoB"))  region = "BE";

    double fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", 0, electronLooseColl, "ELECTRON_HN_TIGHT", 2, "ELECTRON_HN_FAKELOOSE", "mva");

    snu::KParticle Z_candidate;//define Z candidate after shifting energy (SS : better fitting, energy scale up // OS : photon radiation E loss)
    Z_candidate = (lep[0] + lep[1]);
    bool Z_selection = (((Z_candidate.M() - Z_mass) < 10.) && ((Z_mass - Z_candidate.M()) < 10.));

    int Njets = jetVLooseColl.size();
    TString s_njets = "";
    if( Njets == 0 ) s_njets = "JETS0";
    if( Njets != 0 ) s_njets = "JETS";

    if( is_SS ){
      if( Z_selection ){
        FillHist("observed_Z_mass_global"+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_global"+IDsuffix+CFsample+"_"+this_syst, 0., fake_weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+region+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_"+region+IDsuffix+CFsample+"_"+this_syst, 0., fake_weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+s_njets+"_global"+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_"+s_njets+"_global"+IDsuffix+CFsample+"_"+this_syst, 0., fake_weight, 0., 1., 1);
        FillHist("observed_Z_mass_"+s_njets+"_"+region+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), fake_weight, 70., 110., 40);
        FillHist("observed_n_events_"+s_njets+"_"+region+IDsuffix+CFsample+"_"+this_syst, 0., fake_weight, 0., 1., 1);
        if(Njets == 0) FillHist("observed_n_jets_global"+IDsuffix+CFsample+"_"+this_syst, 0., fake_weight, 0., 2., 2);
        if(Njets != 0) FillHist("observed_n_jets_global"+IDsuffix+CFsample+"_"+this_syst, 1., fake_weight, 0., 2., 2);
        if(Njets == 1){
          FillHist("observed_Z_mass_JETS1_global"+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), fake_weight, 70., 110., 40);
          FillHist("observed_n_events_JETS1_global"+IDsuffix+CFsample+"_"+this_syst, 0., fake_weight, 0., 1., 1);	 
          FillHist("observed_Z_mass_JETS1_"+region+IDsuffix+CFsample+"_"+this_syst, Z_candidate.M(), fake_weight, 70., 110., 40);
          FillHist("observed_n_events_JETS1_"+region+IDsuffix+CFsample+"_"+this_syst, 0., fake_weight, 0., 1., 1);
        }
      }// Z selection
    }// is_SS
  }
  
  return;
}


