#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################
declare -a signal_prompt=(

 'WWW' 'WWZ' 'WZZ' 'ZZZ'

 'ttW' 'ttZ'

 'WZ' 'ZZ'

)

declare -a dilep_emu=(
 'HNMoriondLLEmMum_100' 'HNMoriondLLEmMum_1100' 'HNMoriondLLEmMum_200' 'HNMoriondLLEmMum_50' 'HNMoriondLLEmMum_500' 
 'HNMoriondLLEpMup_100' 'HNMoriondLLEpMup_1100' 'HNMoriondLLEpMup_200' 'HNMoriondLLEpMup_50' 'HNMoriondLLEpMup_500'
 'HNMoriondLLMumEm_100' 'HNMoriondLLMumEm_1100' 'HNMoriondLLMumEm_200' 'HNMoriondLLMumEm_50' 'HNMoriondLLMumEm_500'
 'HNMoriondLLMupEp_100' 'HNMoriondLLMupEp_1100' 'HNMoriondLLMupEp_200' 'HNMoriondLLMupEp_50' 'HNMoriondLLMupEp_500'
 'HNMoriondLL_Tchannel_EmMum_100' 'HNMoriondLL_Tchannel_EmMum_1100' 'HNMoriondLL_Tchannel_EmMum_200' 'HNMoriondLL_Tchannel_EmMum_500'
 'HNMoriondLL_Tchannel_EpMup_100' 'HNMoriondLL_Tchannel_EpMup_1100' 'HNMoriondLL_Tchannel_EpMup_200' 'HNMoriondLL_Tchannel_EpMup_500'
 'HNMoriondLL_Tchannel_MumEm_100' 'HNMoriondLL_Tchannel_MumEm_1100' 'HNMoriondLL_Tchannel_MumEm_200' 'HNMoriondLL_Tchannel_MumEm_500'
 'HNMoriondLL_Tchannel_MupEp_100' 'HNMoriondLL_Tchannel_MupEp_1100' 'HNMoriondLL_Tchannel_MupEp_200' 'HNMoriondLL_Tchannel_MupEp_500'
)

declare -a diel=('DYJets' 'WGtoLNuG' 'ZGto2LG' 'TTJets_aMC' 'WW' 'WZ' 'ZZ')

declare -a hntmp=('TTTT' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds'  )

declare -a hn_eetmp=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'WWTo2L2Nu' 'WWToLNuQQ' 'QCD_DoubleEMEnriched_30-40_mgg80toinf' 'QCD_DoubleEMEnriched_30-inf_mgg40to80' 'QCD_DoubleEMEnriched_40-inf_mgg80toinf' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds'  )

declare -a sktmp=('DYJets' 'WJets' 'TT_powheg')
declare -a vv=('WZTo3LNu_powheg' 'ZZTo4L_powheg' 'WGtoLNuG'  'ZGto2LG' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds' 'TG' 'TTG' 'ttWToLNu')

declare -a hn_ee_sig=('WpWpEWK' 'WpWpQCD'  'ZZZ' 'WZZ'  'ww_ds'  'ggZZto4e' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'ggHtoZZ'  'WWG' 'WZG'    'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' )
declare -a hn_ee_sigcf=('DYJets' 'TT_powheg' 'WWTo2L2Nu' )

declare -a hn_ee_type=('DYJets' 'TT_powheg' 'WJets' 'ZGto2LG' 'QCD_Pt-30to50_EMEnriched')


declare -a tmp=(
'DYJets_MG'
'GG_HToMuMu' 
'SingleTop_t'
'TTJets_aMC'
'TTLL_powheg'
'ttZToLL_M-10'
'WWTo2L2Nu_DS'
'ZGto2LG' )

declare -a tmp2=(
'HNEpMup_100' 
'HNEpMup_1100' 
'HNEpMup_1500' 
'HNEpMup_200' 
'HNEpMup_40' 
'HNEpMup_500' 
'HNEpMup_50' 
'HNEpMup_60' 
'HNMoriondLL_Tchannel_EpMup_100' 
'HNMoriondLL_Tchannel_EpMup_1100' 
'HNMoriondLL_Tchannel_EpMup_200' 
'HNMoriondLL_Tchannel_EpMup_500' 
'HN_MuMuMu_1000'
'HN_MuMuMu_10'
'HN_MuMuMu_150'
'HN_MuMuMu_20'
'HN_MuMuMu_200'
 )
