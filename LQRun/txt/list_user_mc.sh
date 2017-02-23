#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################
declare -a cr=('WZTo3LNu_powheg' 'ZZTo4L_powheg')

declare -a hn_muemu=(
 'HN_SSSF_MuEMu_30'
 'HN_SSSF_MuEMu_40'
 'HN_SSSF_MuEMu_50'
 'HN_SSSF_MuEMu_60'
 'HN_SSSF_MuEMu_70'
 'HN_SSSF_MuEMu_100'
 'HN_SSSF_MuEMu_150'
 'HN_SSSF_MuEMu_300'
 'HN_SSSF_MuEMu_500'
 'HN_SSSF_MuEMu_700'
)

declare -a hn_mumue=(
 'HN_SSSF_MuMuE_5'
 'HN_SSSF_MuMuE_10'
 'HN_SSSF_MuMuE_20'
 'HN_SSSF_MuMuE_30'
 'HN_SSSF_MuMuE_40'
 'HN_SSSF_MuMuE_50'
 'HN_SSSF_MuMuE_60'
 'HN_SSSF_MuMuE_70'
 'HN_SSSF_MuMuE_90'
 'HN_SSSF_MuMuE_100'
 'HN_SSSF_MuMuE_150'
 'HN_SSSF_MuMuE_200'
 'HN_SSSF_MuMuE_300'
 'HN_SSSF_MuMuE_400'
 'HN_SSSF_MuMuE_500'
 'HN_SSSF_MuMuE_700'
 'HN_SSSF_MuMuE_1000'
)

declare -a background_prompt=(
 'WWW' 'WWZ' 'WZZ' 'ZZZ'
 'ttW' 'ttZ'
)

declare -a test=('WGtoLNuG' 'WgstarToLNuEE' 'ZGto2LG' 'DYJets_10to50' 'DYJets' 'TTJets_aMC' 'WJets' 'WW' 'WZ' 'ZZ' 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW' 'SingleTop_tW')

declare -a tmplist=('WpWp_qcd_madgraph' 'ZG_llG_MCatNLO' 'ZZ_llnunu_powheg' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llll_powheg' 'ZZ_pythia8' 'ttHnobb_Powheg' 'ttHtobb_Powheg')

declare -a tmpall_mc=('TT_powheg' 'WJets' 'WW'  'WZ' 'ZZ' 'DYJets')

declare -a hn=('DYJets_10to50'  'DYJets' 'WW' 'ZZ' 'WZ' 'TTJets_MG')

declare -a fake_muon=(
 'DYJets_10to50' 'DYJets'
 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW_noHadron' 'SingleTop_tW_noHadron'
 'TTJets_aMC'
 'WGtoLNuG' 'ZGto2LG'
 'WJets'
 'WWW' 'WWZ' 'WZZ' 'ZZZ'
 'WW' 'WZ' 'ZZ'
 'ttH_nonbb' 'ttH_bb' 'ttW' 'ttZ'
)

