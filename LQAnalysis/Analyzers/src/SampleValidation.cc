// $Id: SampleValidation.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQSampleValidation Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "SampleValidation.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (SampleValidation);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
SampleValidation::SampleValidation() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("SampleValidation");
  
  Message("In SampleValidation constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void SampleValidation::InitialiseAnalysis() throw( LQError ) {
  
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


void SampleValidation::ExecuteEvents()throw( LQError ){

  if(isData) return;
  TruthPrintOut();

  if( (k_sample_name.Contains("WR")) ){

    std::vector<snu::KTruth> truthColl;
    eventbase->GetTruthSel()->Selection(truthColl);
    int max = truthColl.size();

    for( int i=2; i<max; i++){
      if( (truthColl.at(i).PdgId() == 9900012) || (truthColl.at(i).PdgId() == 9900016) ) return;
    }

    int hn_index = -999;
    for( int i=2; i<max; i++){
      if( truthColl.at(i).PdgId() == 9900014 ){
        hn_index = i;
        break;
      }
    }
    if( hn_index < 0 ) return;

    int lep_index[2] = {-999, -999};
    for( int i=2; i<max; i++){
      if( abs(truthColl.at(i).PdgId()) == 13 ){
        if( (truthColl.at(i).IndexMother() == truthColl.at(hn_index).IndexMother()) ){
          lep_index[0] = i;
          break;
        }
      }
    }
    if( lep_index[0] < 0 ) return;
    for( int i=2; i<max; i++){
      if( abs(truthColl.at(i).PdgId()) == 13 ){
        if( truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 9900014 ){
          lep_index[1] = i;
          break;
        }
      }
    }
    if( lep_index[1] < 0 ) return;

    std::vector<snu::KParticle> jets_temp, jets; jets_temp.clear(); jets.clear();
    int jet_index[2] = {-999,-999};
    for( int i=2; i<max; i++){
      if( abs(truthColl.at(i).PdgId()) < 7 ){
        if( truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 9900014 ){
          jet_index[0] = i;
          jets_temp.push_back(truthColl.at(i));
          break;
        }
      }
    }
    if( jet_index[0] < 0 ) return;
    for( int i=2; i<max; i++){
      if( abs(truthColl.at(i).PdgId()) < 7 ){
        if( truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 9900014 ){
          if( jet_index[0] != i ){
            jet_index[1] = i;
            jets_temp.push_back(truthColl.at(i));
            break;
          }
        }
      }
    }
    if( jet_index[1] < 0 ) return;

    snu::KParticle N, lep[2], jet[2];
    N = truthColl.at(hn_index);
    lep[0] = truthColl.at(lep_index[0]);
    lep[1] = truthColl.at(lep_index[1]);

    std::vector<snu::KParticle> leptons;
    if(lep[0].Pt() > lep[1].Pt()){
      leptons.push_back(lep[0]);
      leptons.push_back(lep[1]);
    }
    else{
      leptons.push_back(lep[1]);
      leptons.push_back(lep[0]);
    }
    jet[0] = truthColl.at(jet_index[0]);
    jet[1] = truthColl.at(jet_index[1]);
    jets = SortByPtOrder(jets_temp);

    FillHist("Mass_WRoffshell", (jet[0]+jet[1]).M(), 1., 0., 5000., 500);
    FillHist("Pt_WRoffshell", (jet[0]+jet[1]).Pt(), 1., 0., 1000., 100);
    FillHist("Eta_WRoffshell", (jet[0]+jet[1]).Eta(), 1., -5., 5., 50);
    FillHist("Mass_WRonshell", (lep[0]+lep[1]+jet[0]+jet[1]).M(), 1., 0., 5000., 500);
    FillHist("Pt_WRonshell", (lep[0]+lep[1]+jet[0]+jet[1]).Pt(), 1., 0., 1000., 100);
    FillHist("Eta_WRonshell", (lep[0]+lep[1]+jet[0]+jet[1]).Eta(), 1., -5., 5., 50);

    FillHist("Mass_HN", N.M(), 1., 0., 2500., 250);
    FillHist("Pt_HN", N.Pt(), 1., 0., 2000., 200);
    FillHist("Eta_HN", N.Eta(), 1., -5., 5., 50);

    FillHist("Mass_Reco_HN", (lep[1] + jet[0] + jet[1]).M(), 1., 0., 2500., 250);
    FillHist("Pt_Reco_HN", (lep[1] + jet[0] + jet[1]).Pt(), 1., 0., 1000., 100);
    FillHist("Eta_Reco_HN", (lep[1] + jet[0] + jet[1]).Eta(), 1., -5., 5., 50);

    FillHist("Pt_lepton", lep[0].Pt(), 1., 0., 2000., 200);
    FillHist("Pt_lepton", lep[1].Pt(), 1., 0., 2000., 200);
    FillHist("Eta_lepton", lep[0].Eta(), 1., -5., 5., 50);
    FillHist("Eta_lepton", lep[1].Eta(), 1., -5., 5., 50);
    FillHist("Pt_1stlepton", leptons.at(0).Pt(), 1., 0., 2000., 200);
    FillHist("Pt_2ndlepton", leptons.at(1).Pt(), 1., 0., 2000., 200);
    FillHist("Eta_1stlepton", leptons.at(0).Eta(), 1., -5., 5., 50);
    FillHist("Eta_2ndlepton", leptons.at(1).Eta(), 1., -5., 5., 50);

    FillHist("Pt_jet", jet[0].Pt(), 1., 0., 2000., 200);
    FillHist("Pt_jet", jet[1].Pt(), 1., 0., 2000., 200);
    FillHist("Eta_jet", jet[0].Eta(), 1., -5., 5., 50);
    FillHist("Eta_jet", jet[1].Eta(), 1., -5., 5., 50);

    FillHist("Pt_1stjet", jets.at(0).Pt(), 1., 0., 2000., 200);
    FillHist("Pt_2ndjet", jets.at(1).Pt(), 1., 0., 2000., 200);
    FillHist("Eta_1stjet", jets.at(0).Eta(), 1., -5., 5., 50);
    FillHist("Eta_2ndjet", jets.at(1).Eta(), 1., -5., 5., 50);



  }

  if( (k_sample_name.Contains("HNpair")) || (k_sample_name.Contains("Zprime")) ){

    std::vector<snu::KTruth> truthColl;
    eventbase->GetTruthSel()->Selection(truthColl);
    int max = truthColl.size();

    for( int i=2; i<max; i++){
      if( (truthColl.at(i).PdgId() == 9900012) || (truthColl.at(i).PdgId() == 9900016) ) return;
    }

    int hn_index[2] = {-999, -999};
    for( int i=2; i<max; i++){
      if( truthColl.at(i).PdgId() == 9900014 ){
        hn_index[0] = i;
        break;
      }
    }
    if( hn_index[0] < 0 ) return;
    for( int i=2; i<max; i++){
      if( truthColl.at(i).PdgId() == 9900014 ){
        if( hn_index[0] != i ){
          hn_index[1] = i;
          break;
        }
      }
    }
    if( hn_index[1] < 0 ) return;

    int lep_index[2] = {-999, -999};
    for( int i=2; i<max; i++){
      if( abs(truthColl.at(i).PdgId()) == 13 ){
        if( truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 9900014 ){
          lep_index[0] = i;
          break;
        }
      }
    }
    if( lep_index[0] < 0 ) return;
    for( int i=2; i<max; i++){
      if( abs(truthColl.at(i).PdgId()) == 13 ){
        if( truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 9900014 ){
          if( lep_index[0] != i ){
            lep_index[1] = i;
            break;
          }
        }
      }
    }
    if( lep_index[1] < 0 ) return;

    std::vector<snu::KParticle> jets_temp, jets; jets_temp.clear(); jets.clear();
    int jet_index[2][2] = {{-999,-999},{-999,-999}};
    for( int i=2; i<max; i++){
      if( abs(truthColl.at(i).PdgId()) < 7 ){
        if( (truthColl.at(lep_index[0]).IndexMother() == truthColl.at(i).IndexMother()) ){
          jet_index[0][0] = i;
          jets_temp.push_back(truthColl.at(i));
          break;
        }
      }
    }
    if( jet_index[0][0] < 0 ) return;
    for( int i=2; i<max; i++){
      if( abs(truthColl.at(i).PdgId()) < 7 ){
        if( (truthColl.at(lep_index[0]).IndexMother() == truthColl.at(i).IndexMother()) ){
          if( jet_index[0][0] != i ){
            jet_index[0][1] = i;
            jets_temp.push_back(truthColl.at(i));
            break;
          }
        }
      }
    }
    if( jet_index[0][1] < 0 ) return;
    for( int i=2; i<max; i++){
      if( abs(truthColl.at(i).PdgId()) < 7 ){
        if( (truthColl.at(lep_index[1]).IndexMother() == truthColl.at(i).IndexMother()) ){
          jet_index[1][0] = i;
          jets_temp.push_back(truthColl.at(i));
          break;
        }
      }
    }
    if( jet_index[1][0] < 0 ) return;
    for( int i=2; i<max; i++){
      if( abs(truthColl.at(i).PdgId()) < 7 ){
        if( (truthColl.at(lep_index[1]).IndexMother() == truthColl.at(i).IndexMother()) ){
          if( jet_index[1][0] != i ){
            jet_index[1][1] = i;
            jets_temp.push_back(truthColl.at(i));
            break;
          }
        }
      }
    }
    if( jet_index[1][1] < 0 ) return;

    snu::KParticle N[2], lep[2], jet[2][2];
    N[0] = truthColl.at(hn_index[0]);
    N[1] = truthColl.at(hn_index[1]);
    lep[0] = truthColl.at(lep_index[0]);
    lep[1] = truthColl.at(lep_index[1]);

    std::vector<snu::KParticle> leptons;
    if(lep[0].Pt() > lep[1].Pt()){
      leptons.push_back(lep[0]);
      leptons.push_back(lep[1]);
    }
    else{
      leptons.push_back(lep[1]);
      leptons.push_back(lep[0]);
    }
    jet[0][0] = truthColl.at(jet_index[0][0]);
    jet[0][1] = truthColl.at(jet_index[0][1]);
    jet[1][0] = truthColl.at(jet_index[1][0]);
    jet[1][1] = truthColl.at(jet_index[1][1]);
    jets = SortByPtOrder(jets_temp);

    FillHist("Mass_Zprime", (N[0]+N[1]).M(), 1., 0., 5000., 500);
    FillHist("Mass_Reco_Zprime", (lep[0]+jet[0][0]+jet[0][1]+lep[1]+jet[1][0]+jet[1][1]).M(), 1., 0., 5000., 500);

    FillHist("Mass_HN", N[0].M(), 1., 0., 2500., 250);
    FillHist("Mass_HN", N[1].M(), 1., 0., 2500., 250);
    FillHist("Pt_HN", N[0].Pt(), 1., 0., 2000., 200);
    FillHist("Pt_HN", N[1].Pt(), 1., 0., 2000., 200);
    FillHist("Eta_HN", N[0].Eta(), 1., -5., 5., 50);
    FillHist("Eta_HN", N[1].Eta(), 1., -5., 5., 50);

    FillHist("Mass_Reco_HN", (lep[0] + jet[0][0] + jet[0][1]).M(), 1., 0., 2500., 250);
    FillHist("Mass_Reco_HN", (lep[1] + jet[1][0] + jet[1][1]).M(), 1., 0., 2500., 250);
    FillHist("Pt_Reco_HN", (lep[0] + jet[0][0] + jet[0][1]).Pt(), 1., 0., 2000., 200);
    FillHist("Pt_Reco_HN", (lep[1] + jet[1][0] + jet[1][1]).Pt(), 1., 0., 2000., 200);
    FillHist("Eta_Reco_HN", (lep[0] + jet[0][0] + jet[0][1]).Eta(), 1., -5., 5., 50);
    FillHist("Eta_Reco_HN", (lep[1] + jet[1][0] + jet[1][1]).Eta(), 1., -5., 5., 50);

    FillHist("Pt_lepton", lep[0].Pt(), 1., 0., 2000., 200);
    FillHist("Pt_lepton", lep[1].Pt(), 1., 0., 2000., 200);
    FillHist("Eta_lepton", lep[0].Eta(), 1., -5., 5., 50);
    FillHist("Eta_lepton", lep[1].Eta(), 1., -5., 5., 50);
    FillHist("Pt_1stlepton", leptons.at(0).Pt(), 1., 0., 2000., 200);
    FillHist("Pt_2ndlepton", leptons.at(1).Pt(), 1., 0., 2000., 200);
    FillHist("Eta_1stlepton", leptons.at(0).Eta(), 1., -5., 5., 50);
    FillHist("Eta_2ndlepton", leptons.at(1).Eta(), 1., -5., 5., 50);

    FillHist("Pt_jet", jet[0][0].Pt(), 1., 0., 2000., 200);
    FillHist("Pt_jet", jet[0][1].Pt(), 1., 0., 2000., 200);
    FillHist("Pt_jet", jet[1][0].Pt(), 1., 0., 2000., 200);
    FillHist("Pt_jet", jet[1][1].Pt(), 1., 0., 2000., 200);
    FillHist("Eta_jet", jet[0][0].Eta(), 1., -5., 5., 50);
    FillHist("Eta_jet", jet[0][1].Eta(), 1., -5., 5., 50);
    FillHist("Eta_jet", jet[1][0].Eta(), 1., -5., 5., 50);
    FillHist("Eta_jet", jet[1][1].Eta(), 1., -5., 5., 50);

    FillHist("Pt_1stjet", jets.at(0).Pt(), 1., 0., 2000., 200);
    FillHist("Pt_2ndjet", jets.at(1).Pt(), 1., 0., 2000., 200);
    FillHist("Pt_3rdjet", jets.at(2).Pt(), 1., 0., 2000., 200);
    FillHist("Pt_4thjet", jets.at(3).Pt(), 1., 0., 2000., 200);
    FillHist("Eta_1stjet", jets.at(0).Eta(), 1., -5., 5., 50);
    FillHist("Eta_2ndjet", jets.at(1).Eta(), 1., -5., 5., 50);
    FillHist("Eta_3rdjet", jets.at(2).Eta(), 1., -5., 5., 50);
    FillHist("Eta_4thjet", jets.at(3).Eta(), 1., -5., 5., 50);

  }
/*
  TruthPrintOut();
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);

  int max = truthColl.size();
  std::vector<int> hn_index, lep1_index, lep2_index, q1_index, q2_index;
  for( int i=2 ; i<max ; i++){
    if( ((truthColl.at(i).PdgId()) == 9900012) || ((truthColl.at(i).PdgId()) == 90) ){
      hn_index.push_back(i);
      break;
    }
  }
  if(hn_index.size() ==0) return;

  for( int i=2; i<max ; i++){
    if( (fabs(truthColl.at(i).PdgId()) == 13 ) && truthColl.at(hn_index.at(0)).IndexMother()==truthColl.at(i).IndexMother()){
      lep1_index.push_back(i);
      break;
    }
  }
  for( int i=2; i<max ; i++){
    if( (fabs(truthColl.at(i).PdgId()) == 13 )  && ( truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 9900012 ||truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 90 )){
      lep2_index.push_back(i);
      break;
    }
  }

  for( int i=2; i<max ; i++){
    if( (0 < truthColl.at(i).PdgId() && truthColl.at(i).PdgId() < 6)  && ( truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 9900012 ||truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 24 ||truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 90)){
      q1_index.push_back(i);
      break;
    }
  }
  for( int i=2; i<max ; i++){
    if( (-6 < truthColl.at(i).PdgId() && truthColl.at(i).PdgId() < 0)  && ( truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 9900012 ||truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 24 ||truthColl.at(truthColl.at(i).IndexMother()).PdgId() == 90) ){
      q2_index.push_back(i);
      break;
    }
  }
  if(lep1_index.size() * lep2_index.size() * q1_index.size() * q2_index.size() == 0) return;

  snu::KTruth lep[2], jet[2];
  lep[0] = truthColl.at(lep1_index.at(0));
  lep[1] = truthColl.at(lep2_index.at(0));
  jet[0] = truthColl.at(q1_index.at(0));
  jet[1] = truthColl.at(q2_index.at(0));


  if(lep[0].Pt() < lep[1].Pt()){
    snu::KTruth temp;
    temp = lep[0];
    lep[0] = lep[1];
    lep[1] = temp;
  } 
  if(jet[0].Pt() < jet[1].Pt()){
    snu::KTruth temp;
    temp = jet[0];
    jet[0] = jet[1];
    jet[1] = temp;
  }

  std::vector<snu::KJet> jets = GetJets("JET_NOLEPTONVETO", 20, 5.0);
  FillHist("JETS", jets.size(), 1., 0., 10., 10);

  std::vector<snu::KMuon> muons = GetMuons("MUON_HN_TIGHT", false);
  FillHist("NEVENTS", 0., 1., 0., 1., 1);
  FillHist("NVERTICES", eventbase->GetEvent().nVertices(), 1., 0., 50., 50);
*/
  return;
}// End of execute event loop
  


void SampleValidation::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void SampleValidation::BeginCycle() throw( LQError ){
  
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

SampleValidation::~SampleValidation() {
  
  Message("In SampleValidation Destructor" , INFO);
  
}


void SampleValidation::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void SampleValidation::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this SampleValidationCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void SampleValidation::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}

std::vector<snu::KParticle> SampleValidation::SortByPtOrder( std::vector<snu::KParticle> jets ){
  int const N=jets.size();
  int order[N];
  for(int i=0; i<N; i++){
    order[i] = i;
  }

  std::vector<snu::KParticle> jets_temp;
  bool can_stop = true;
  for(int i=0; i<jets.size()-1; i++){
    snu::KParticle this_jet = jets.at(i);
    snu::KParticle next_jet = jets.at(i+1);

    if(this_jet.Pt()<next_jet.Pt()){
      order[i] = i+1;
      order[i+1] = i;
      can_stop = false;
    }
  }
  for(int i=0; i<jets.size(); i++){
    jets_temp.push_back(jets.at(order[i]));
  }
  if(!can_stop) SortByPtOrder( jets_temp );
  else return jets_temp;
}

