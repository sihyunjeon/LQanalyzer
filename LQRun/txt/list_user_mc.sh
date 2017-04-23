#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

declare -a signal_prompt_os=('ggHtoZZ' 'ggZZto2e2mu' 'TTG' 'ttWToLNu' 'ttZToLL_M-10' 'ttZToLL_M-1to10' 'vbhHtoZZ' 'WgstarToLNuEE' 'WgstarToLNuMuMu' 'WWG' 'WWW' 'WWZ' 'WZG' 'WZZ' 'WZ' 'ZGto2LG' 'ZZ') 

declare -a example=("WW" "WZ" "HNMupMup_100" "DYJets")

declare -a tmplist=('WpWp_qcd_madgraph' 'ZG_llG_MCatNLO' 'ZZ_llnunu_powheg' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llll_powheg' 'ZZ_pythia8' 'ttHnobb_Powheg' 'ttHtobb_Powheg')

declare -a tmpall_mc=('TTJets_aMC' 'WJets' 'WW'  'WZ' 'ZZ' 'DYJets')

declare -a closure=('TTJets_aMC' 'WJets' 'DYJets')

declare -a chargeflip=(
 'DYtoEE'
 'DYJets_10to50'
 'DYJets'
 'TTJets_aMC'
 'TTLL_powheg'
 'TTLJ_powheg'
 'WJets'
 'WJets_MG'
)

declare -a chargeflip_powheg=(
 'DYtoEE'
 'TTLL_powheg'
 'TTLJ_powheg'
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

declare -a signal_prompt=(

 'WWW' 'WWZ' 'WZZ' 'ZZZ'

 'ttW' 'ttZ'

 'WZ' 'ZZ'

)

declare -a control_prompt=(

 'ttWToLNu'
 'ttZToLL_M-1to10' 'ttZ'
 'ttH_nonbb'

 'WZTo3LNu_amcatnlo'
 'ZZTo4L_powheg'

 'WWW' 'WWZ' 'WZZ' 'ZZZ'

 'WGtoLNuG' 'ZGto2LG'

 'ggZZto2e2mu'
)

declare -a pu_dilepton_list=('DYJets_10to50' 'DYJets' 'WJets' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg' 'DYJets_MG_10to50' 'DYJets_MG' 'TTJets_aMC' )

declare -a hntmp=('TTTT' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds'  )

declare -a hn_eetmp=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'WWTo2L2Nu' 'WWToLNuQQ' 'QCD_DoubleEMEnriched_30-40_mgg80toinf' 'QCD_DoubleEMEnriched_30-inf_mgg40to80' 'QCD_DoubleEMEnriched_40-inf_mgg80toinf' 'TTTT' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds'  )

declare -a sktmp=('DYJets' 'WJets' 'TT_powheg')

declare -a hn_ee_sig=('WpWpEWK' 'WpWpQCD'  'ZZZ' 'WZZ'  'ww_ds'  'ggZZto4e' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'ggHtoZZ'  'WWG' 'WZG'    'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' )
declare -a hn_ee_sigcf=('DYJets' 'TT_powheg' 'WWTo2L2Nu' )

declare -a hn_ee_type=('DYJets' 'TT_powheg' 'WJets' 'ZGto2LG' 'QCD_Pt-30to50_EMEnriched')
