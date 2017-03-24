#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

declare -a example=("WW" "WZ" "HNMupMup_100" "DYJets")

declare -a tmplist=('WpWp_qcd_madgraph' 'ZG_llG_MCatNLO' 'ZZ_llnunu_powheg' 'ZZ_llqq_MCatNLO' 'ZZ_llll_MCatNLO' 'ZZ_llll_powheg' 'ZZ_pythia8' 'ttHnobb_Powheg' 'ttHtobb_Powheg')

declare -a tmpall_mc=('TTJets_aMC' 'WJets' 'WW'  'WZ' 'ZZ' 'DYJets')

declare -a hn=('DYJets_10to50'  'DYJets' 'WW' 'ZZ' 'WZ' 'TTJets_MG')

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

)

