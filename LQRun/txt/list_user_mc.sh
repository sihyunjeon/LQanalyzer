#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

declare -a test=(
 'WWW' 'WWZ' 'WZZ' 'ZZZ'
)

declare -a hn_mumumu=(
 'HN_MuMuMu_5'
 'HN_MuMuMu_20'
 'HN_MuMuMu_30'
 'HN_MuMuMu_40'
 'HN_MuMuMu_50'
 'HN_MuMuMu_60'
 'HN_MuMuMu_70'
 'HN_MuMuMu_90'
 'HN_MuMuMu_100'
 'HN_MuMuMu_150'
 'HN_MuMuMu_200'
 'HN_MuMuMu_300'
 'HN_MuMuMu_400'
 'HN_MuMuMu_500'
 'HN_MuMuMu_700'
 'HN_MuMuMu_1000'
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

declare -a background=(
 'DYJets_10to50' 'DYJets' 
 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW_noHadron' 'SingleTbar_tW' 'SingleTop_tW_noHadron' 'SingleTop_tW'
 'TTJets_aMC' 'TTLL_powheg' 'TTLJ_powheg'
 'VBF_HToMuMu' 'GG_HToMuMu'
 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG'
 'WJets'
 'WWW' 'WWZ' 'WZZ' 'ZZZ'
 'WW' 'WZ' 'ZZ' 'WpWpEWK' 'WpWpQCD'
 'ttH_nonbb' 'ttH_bb' 'ttW' 'ttZ'
)

declare -a background_cr=(
 'WGtoLNuMM'
 'WZTo3LNu_powheg' 'ZZTo4L_powheg'
 'WWW' 'WWZ' 'WZZ' 'ZZZ'
 'ttW' 'ttZ'
)

declare -a fake_muon=(
 'DYJets_10to50' 'DYJets'
 'GG_HToMuMu' 'VBF_HToMuMu'
 'SingleTop_s' 'SingleTbar_t' 'SingleTop_t' 'SingleTbar_tW' 'SingleTop_tW'
 'TTLL_powheg' 'TTLJ_powheg'
 'WGtoLNuG' 'WGtoLNuEE' 'WGtoLNuMM' 'ZGto2LG'
 'WJets'
 'WWW' 'WWZ' 'WZZ' 'ZZZ'
 'WW' 'WZ' 'ZZ'
 'ttH_nonbb' 'ttH_bb' 'ttW' 'ttZ'
)

declare -a fake_qcd=(
 'QCD_Pt-15to20_MuEnriched'
 'QCD_Pt-20to30_MuEnriched'
 'QCD_Pt-30to50_MuEnriched'
 'QCD_Pt-50to80_MuEnriched'
 'QCD_Pt-80to120_MuEnriched'
 'QCD_Pt-120to170_MuEnriched'
 'QCD_Pt-170to300_MuEnriched'
 'QCD_Pt-300to470_MuEnriched'
 'QCD_Pt-470to600_MuEnriched'
 'QCD_Pt-600to800_MuEnriched'
 'QCD_Pt-800to1000_MuEnriched'
 'QCD_Pt-1000toInf_MuEnriched'
)
