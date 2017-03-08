#include "SSSFMuMuEPlots.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

SSSFMuMuEPlots::SSSFMuMuEPlots(TString name): StdPlots(name){

  TH1::SetDefaultSumw2(true);

  map_sig["h_NMuons"]				= new TH1D("h_NMuons_"			+ name,"number of muon",5,0,5);
  map_sig["h_NElectrons"]	 		= new TH1D("h_NElectrons_"		+ name,"number of electron",5,0,5);

  map_sig["h_MuonPt"] 				= new TH1D("h_MuonPt_"			+ name,"Pt of muon",500,0,500);
  map_sig["h_MuonEta"] 				= new TH1D("h_MuonEta_"			+ name,"Eta of muon",60,-3,3);
  map_sig["h_MuondZ"] 				= new TH1D("h_MuondZ_"  		+ name,"dZ of muon",100,-0.5,0.5);
  map_sig["h_MuondXY"] 				= new TH1D("h_MuondXY_"			+ name,"dXY of muon",200,-1,1);
  map_sig["h_MuondXYSig"] 			= new TH1D("h_MuondXYSig_"			+ name,"dXYSig of muon",200,-10,10);
  map_sig["h_MuonRelIso04"] 			= new TH1D("h_MuonRelIso04_"		+ name,"RelIso04 of muon",100,0,1);
  map_sig["h_leadingMuonPt"]                   	= new TH1D("h_leadingMuonPt_"  	        + name,"Pt of leading muon",500,0,500);
  map_sig["h_leadingMuonEta"]                  	= new TH1D("h_leadingMuonEta_" 		+ name,"Eta of leading muon",60,-3,3);
  map_sig["h_leadingMuondZ"]                   	= new TH1D("h_leadingMuondZ_"           + name,"dZ of leading muon",100,-0.5,0.5);
  map_sig["h_leadingMuondXY"]                  	= new TH1D("h_leadingMuondXY_"          + name,"dXY of leading muon",200,-1,1);
  map_sig["h_leadingMuondXYSig"]               	= new TH1D("h_leadingMuondXYSig_"           + name,"dXYSig of leading muon",200,-10,10);
  map_sig["h_leadingMuonRelIso04"]             	= new TH1D("h_leadingMuonRelIso04_"     + name,"RelIso04 of leading muon",100,0,1);
  map_sig["h_secondMuonPt"]            	        = new TH1D("h_secondMuonPt_"            + name,"Pt of second muon",500,0,500);
  map_sig["h_secondMuonEta"]           	        = new TH1D("h_secondMuonEta_"           + name,"Eta of second muon",60,-3,3);
  map_sig["h_secondMuondZ"]            	        = new TH1D("h_secondMuondZ_"            + name,"dZ of second muon",100,-0.5,0.5);
  map_sig["h_secondMuondXY"]           	        = new TH1D("h_secondMuondXY_"           + name,"dXY of second muon",200,-1,1);
  map_sig["h_secondMuondXYSig"]        	        = new TH1D("h_secondMuondXYSig_"            + name,"dXYSig of second muon",200,-10,10);
  map_sig["h_secondMuonRelIso04"]      	        = new TH1D("h_secondMuonRelIso04_"      + name,"RelIso04 of second muon",100,0,1);
  map_sig["h_ElectronPt"]                       = new TH1D("h_ElectronPt_"            	+ name,"Pt of electron",500,0,500);
  map_sig["h_ElectronEta"]                      = new TH1D("h_ElectronEta_"           	+ name,"Eta of electron",60,-3,3);
  map_sig["h_ElectrondZ_b"]                    	= new TH1D("h_ElectrondZ_b_"            + name,"dZ_b of electron",200,-1,1);
  map_sig["h_ElectrondZ_e"]                    	= new TH1D("h_ElectrondZ_e_"          	+ name,"dZ_e of electron",200,-1,1);
  map_sig["h_ElectrondXY_b"]                    = new TH1D("h_ElectrondXY_b_"           + name,"dXY_b of electron",100,-0.5,0.5);
  map_sig["h_ElectrondXY_e"]                    = new TH1D("h_ElectrondXY_e_"           + name,"dXY_e of electron",100,-0.5,0.5);
  map_sig["h_ElectronRelIso03_b"]		= new TH1D("h_ElectronRelIso03_b_"	+ name,"RelIso03_b of electron",100,0,1);
  map_sig["h_ElectronRelIso03_e"]               = new TH1D("h_ElectronRelIso03_e_"      + name,"RelIso03_e of electron",100,0,1);

  map_sig["h_sumcharge"] 	 	        = new TH1D("h_sumcharge_" 		+ name,"chargesum of all lepton",5,-2,3);
  map_sig["h_NJets"]       		        = new TH1D("h_NJets_"		        + name,"number of jet",10,0,10);
  map_sig["h_JetPt"]     		        = new TH1D("h_JetPt_"           	+ name,"Pt of jet",500,0,500);
  map_sig["h_JetEta"]               		= new TH1D("h_JetEta_"        		+ name,"Eta of jet",60,-3,3);
  map_sig["h_PFMET"]  		                = new TH1D("h_PFMET_"               	+ name,"missing ET",500,0.0,500.0);
  map_sig["h_PFMETPhi"]         		= new TH1D("h_PFMETPhi_"           	+ name,"Phi of missing ET",80,-4,4);
  map_sig["h_nVertices"]	                = new TH1D("h_nVertices_" 	        + name,"number of even vertices",60,0,60);
  map_sig["h_HT"]				= new TH1D("h_HT_"			+ name,"Pt sum of jet",500,0,500);

}

void SSSFMuMuEPlots::Fill(snu::KEvent ev, std::vector<snu::KMuon>& muons, std::vector<snu::KElectron>& electrons, std::vector<snu::KJet>& jets, Double_t weight) {
  
  if(! ((electrons.size() == 1 && muons.size() ==2) || (muons.size() ==2) || (electrons.size() == 1 && muons.size() ==1) )) return;
  Fill("h_NMuons" ,muons.size(), weight);
  Fill("h_NElectrons" ,electrons.size(), weight);

  int imu=0;
  for(std::vector<snu::KMuon>::iterator muit = muons.begin(); muit != muons.end(); muit++, imu++){
    Fill("h_MuonPt", muit->Pt(),weight);
    Fill("h_MuonEta",muit->Eta(),weight);
    Fill("h_MuondZ", muit->dZ(),weight);
    Fill("h_MuondXY", muit->dXY(),weight);
    Fill("h_MuondXYSig", muit->dXYSig(),weight);
    Fill("h_MuonRelIso04", muit->RelIso04(), weight);
    if(imu ==0) {
      Fill("h_leadingMuonPt", muit->Pt(),weight);
      Fill("h_leadingMuonEta",muit->Eta(),weight);
      Fill("h_leadingMuondZ", muit->dZ(),weight);
      Fill("h_leadingMuondXY", muit->dXY(),weight);
      Fill("h_leadingMuondXYSig", muit->dXYSig(),weight);
      Fill("h_leadingMuonRelIso04", muit->RelIso04(), weight);
    }
    if(imu ==1) {
      Fill("h_secondMuonPt", muit->Pt(),weight);
      Fill("h_secondMuonEta",muit->Eta(),weight);
      Fill("h_secondMuondZ", muit->dZ(),weight);
      Fill("h_secondMuondXY", muit->dXY(),weight);
      Fill("h_secondMuondXYSig", muit->dXYSig(),weight);
      Fill("h_secondMuonRelIso04", muit->RelIso04(), weight);
    }
  }

  int iel=0, n_b=0, n_e=0;
  for(std::vector<snu::KElectron>::iterator elit = electrons.begin(); elit != electrons.end(); elit++, iel++){
    if(iel ==0){
      Fill("h_ElectronPt", elit->Pt(),weight);
      Fill("h_ElectronEta",elit->Eta(),weight);
      if(((elit->Eta()) < 1.479)){
        Fill("h_ElectrondZ_b", elit->dz(),weight);
        Fill("h_ElectrondXY_b", elit->dxy(),weight);
        Fill("h_ElectronRelIso03_b", elit->PFRelIso(0.3), weight);
      }
      else{
        Fill("h_ElectrondZ_e", elit->dz(),weight);
        Fill("h_ElectrondXY_e", elit->dxy(),weight);
        Fill("h_ElectronRelIso03_e", elit->PFRelIso(0.3), weight);
      }
    }
  }

  int sum_charge=0;
  for(std::vector<snu::KElectron>::iterator elit = electrons.begin(); elit != electrons.end(); elit++){
    sum_charge+= elit->Charge();
  }
  for(std::vector<snu::KMuon>::iterator muit = muons.begin(); muit != muons.end(); muit++){
    sum_charge+= muit->Charge();
  }
  Fill("h_sumcharge",sum_charge, weight);
   
  Fill("h_NJets",jets.size(), weight);
  
  Fill("h_PFMET",ev.PFMET(), weight);
  Fill("h_PFMETPhi",ev.METPhi(snu::KEvent::pfmet), weight);

  Fill("h_nVertices", ev.nVertices(), weight);
 
  double HT=0.; 
  for(UInt_t j=0; j < jets.size(); j++){
    Fill("h_JetPt", jets[j].Pt(),weight);
    Fill("h_JetEta",jets[j].Eta(),weight);
    HT += jets[j].Pt();
  }
  Fill("h_HT",HT, weight);

  return;
}/// End of Fill



void SSSFMuMuEPlots::Write() {
 
  for(map<TString, TH1*>::iterator it = map_sig.begin(); it != map_sig.end(); it++){
    it->second->Write();
  }

  for(map<TString, TH2*>::iterator it = map_sig2.begin(); it != map_sig2.end(); it++){
    it->second->Write();
  }


  for(map<TString, TH3*>::iterator it = map_sig3.begin(); it != map_sig3.end(); it++){
    it->second->Write();
  }

}


SSSFMuMuEPlots::SSSFMuMuEPlots(): StdPlots(){
}


/**
 * Copy constructor.
 */
SSSFMuMuEPlots::SSSFMuMuEPlots(const SSSFMuMuEPlots& sp): StdPlots(sp)
{
  for(std::map<TString, TH1*>::iterator mit = map_sig.begin(); mit != map_sig.end() ; mit++){
    std::map<TString, TH1*>::iterator mit2 = sp.GetMap().find(mit->first);
    mit->second = mit2->second;
  }

  for(std::map<TString, TH2*>::iterator mit = map_sig2.begin(); mit != map_sig2.end() ; mit++){
    std::map<TString, TH2*>::iterator mit2 = sp.GetMap2().find(mit->first);
    mit->second = mit2->second;
  }

  for(std::map<TString, TH3*>::iterator mit = map_sig3.begin(); mit != map_sig3.end() ; mit++){
    std::map<TString, TH3*>::iterator mit2 = sp.GetMap3().find(mit->first);
    mit->second = mit2->second;
  }
  
}


SSSFMuMuEPlots& SSSFMuMuEPlots::operator= (const SSSFMuMuEPlots& sp)
{
  if (this != &sp) {

    for(std::map<TString, TH1*>::iterator mit = map_sig.begin(); mit != map_sig.end() ; mit++){
      std::map<TString, TH1*>::iterator mit2 = sp.GetMap().find(mit->first);
      mit->second = mit2->second;
    }
    
    for(std::map<TString, TH2*>::iterator mit = map_sig2.begin(); mit != map_sig2.end() ; mit++){
      std::map<TString, TH2*>::iterator mit2 = sp.GetMap2().find(mit->first);
      mit->second = mit2->second;
    }

    for(std::map<TString, TH3*>::iterator mit = map_sig3.begin(); mit != map_sig3.end() ; mit++){
      std::map<TString, TH3*>::iterator mit2 = sp.GetMap3().find(mit->first);
      mit->second = mit2->second;
    }
  }
  return *this;
}

SSSFMuMuEPlots::~SSSFMuMuEPlots() {
   for(std::map<TString, TH1*>::iterator mit = map_sig.begin(); mit != map_sig.end() ; mit++){
     delete mit->second ;
  }

   for(std::map<TString, TH2*>::iterator mit = map_sig2.begin(); mit != map_sig2.end() ; mit++){
     delete mit->second ;
   }

   for(std::map<TString, TH3*>::iterator mit = map_sig3.begin(); mit != map_sig3.end() ; mit++){
     delete mit->second ;
   }
   
}

void SSSFMuMuEPlots::Fill(TString name, double value, double w){
  std::map<TString, TH1*>::iterator it = map_sig.find(name);
  if(it!= map_sig.end())   it->second->Fill(value, w);

  else cout << name << " not found in map_sig" << endl;
  return;
}
 
void SSSFMuMuEPlots::Fill(TString name, double value1, double value2, double w){
   std::map<TString, TH2*>::iterator it = map_sig2.find(name);
   if(it!= map_sig2.end()) it->second->Fill(value1, value2, w);
   else cout << name << " not found in map_sig" << endl;
   return;
 }



void SSSFMuMuEPlots::Fill(TString name, double value1, double value2, double value3, double w){
  std::map<TString, TH3*>::iterator it = map_sig3.find(name);
  if(it!= map_sig3.end()) it->second->Fill(value1, value2, value3, w);
  else cout << name << " not found in map_sig" << endl;
  return;
}

