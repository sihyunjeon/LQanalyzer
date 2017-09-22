#!/bin/sh

########################
### SAMPLE LIST ########## 
#######################

declare -a  qcd_mu_prompt=('QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 
'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched' 'TT_powheg' 'WJets'  'DYJets')
declare -a hnfail=('HNMumEp_50' 'HNMumEp_500' 'HNMumEp_60' 'HNMumMum_100' 'HNMumMum_1100' 'HNMumMum_1500' 'HNMumMum_200' 'HNMumMum_40' 'HNMumMum_50' 'HNMumMum_500' 'HNMumMum_60' 'HNMoriondLL_Tchannel_MumMum_100' 'HNMoriondLL_Tchannel_MumMum_1100' 'HNMoriondLL_Tchannel_MumMum_200' 'HNMoriondLL_Tchannel_MumMum_500' 'HNMumMup_100' 'HNMumMup_1100' 'HNMumMup_1500' 'HNMumMup_200' 'HNMumMup_40' 'HNMumMup_50' 'HNMumMup_500' 'HNMumMup_60' 'HNMupEm_100' 'HNMupEm_1100' 'HNMupEm_1500' 'HNMupEm_200' 'HNMupEm_40' 'HNMupEm_50' 'HNMupEm_500' 'HNMupEm_60' 'HNMupEp_100')

declare -a hn_ss_mm=(
 'HNMumMum_40' 'HNMumMum_60' 'HNMumMum_100' 'HNMumMum_200' 'HNMumMum_600' 'HNMumMum_1000' 'HNMumMum_1400'

)

declare -a hn_ss=(
 'HNMumEm_40' 'HNMumEm_50' 'HNMumEm_60' 'HNMumEm_70' 'HNMumEm_90' 'HNMumEm_100' 'HNMumEm_125' 'HNMumEm_150' 'HNMumEm_200' 'HNMumEm_300' 'HNMumEm_400' 'HNMumEm_600' 'HNMumEm_1000' 'HNMumEm_1400'
 'HNMupEp_40' 'HNMupEp_50' 'HNMupEp_60' 'HNMupEp_70' 'HNMupEp_90' 'HNMupEp_100' 'HNMupEp_125' 'HNMupEp_150' 'HNMupEp_200' 'HNMupEp_300' 'HNMupEp_400' 'HNMupEp_600' 'HNMupEp_1000' 'HNMupEp_1400'
 'HNEmMum_40' 'HNEmMum_50' 'HNEmMum_60' 'HNEmMum_70' 'HNEmMum_90' 'HNEmMum_100' 'HNEmMum_125' 'HNEmMum_150' 'HNEmMum_200' 'HNEmMum_300' 'HNEmMum_400' 'HNEmMum_600' 'HNEmMum_1000' 'HNEmMum_1400'
 'HNEpMup_40' 'HNEpMup_50' 'HNEpMup_60' 'HNEpMup_70' 'HNEpMup_90' 'HNEpMup_100' 'HNEpMup_125' 'HNEpMup_150' 'HNEpMup_200' 'HNEpMup_300' 'HNEpMup_400' 'HNEpMup_600' 'HNEpMup_1000' 'HNEpMup_1400'

 'HNMumMum_40' 'HNMumMum_60' 'HNMumMum_100' 'HNMumMum_200' 'HNMumMum_600' 'HNMumMum_1000' 'HNMumMum_1400'
 'HNMupMup_40' 'HNMupMup_60' 'HNMupMup_100' 'HNMupMup_200' 'HNMupMup_600' 'HNMupMup_1000' 'HNMupMup_1400'
)

declare -a hn_bkg=(
 'SingleTbar_tW' 'SingleTop_tW' 'ttW' 'ttH_nonbb' 'ttH_bb' 'TG' 'TTG' 'tZq' 'ttZToLL_M-1to10' 'ttZToLL_M-10'
 'WZTo3LNu_amcatnlo' 'ZZTo4L_powheg'
 'ZZZ' 'WZZ' 'WWZ' 'WWW'
 'WGtoLNuG' 'ZGto2LG'
 'ggHtoZZ' 
 'ww_ds' 'WZG' 'WWG'
 'ggZZto2e2mu' 'ggZZto2mu2tau' 'ggZZto2e2tau' 'ggZZto4tau'
 'WpWpEWK' 'WpWpQCD'
)

declare -a dilepton_list2=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ZZZ' 'WZZ' 'WWZ' 'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg' 'WWTo2L2Nu' 'WWToLNuQQ' 'TG' 'TTG' 'ggHtoWW' 'ggHtoZZ' 'vbfHtoWW' 'tZq' 'ww_ds' 'ttZToLL_M-1to10' 'WZG' 'WWG' 'ggZZto4e' 'ttWToLNu' 'ttZToLL_M-10' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-300toInf_EMEnriched' 'QCD_DoubleEMEnriched_30-40_mgg80toinf' 'QCD_Pt-30to50_EMEnriched' 'QCD_DoubleEMEnriched_30-inf_mgg40to80' 'QCD_DoubleEMEnriched_40-inf_mgg80toinf' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched' 'qcd_15to20_bctoe' 'qcd_170to250_bctoe' 'qcd_20to30_bctoe' 'qcd_250toinf_bctoe' 'qcd_30to80_bctoe' 'qcd_80to170_bctoe')

declare -a hn_mme=('HN_SSSF_MuMuE_1000' 'HN_SSSF_MuMuE_100' 'HN_SSSF_MuMuE_10' 'HN_SSSF_MuMuE_150' 'HN_SSSF_MuMuE_200' 'HN_SSSF_MuMuE_20' 'HN_SSSF_MuMuE_300' 'HN_SSSF_MuMuE_30' 'HN_SSSF_MuMuE_400' 'HN_SSSF_MuMuE_40' 'HN_SSSF_MuMuE_500' 'HN_SSSF_MuMuE_50' 'HN_SSSF_MuMuE_5' 'HN_SSSF_MuMuE_60' 'HN_SSSF_MuMuE_700' 'HN_SSSF_MuMuE_70' 'HN_SSSF_MuMuE_90')

declare -a qcd_ee_hn=('qcd_15to20_bctoe' 'qcd_170to250_bctoe' 'qcd_20to30_bctoe' 'qcd_250toinf_bctoe' 'qcd_30to80_bctoe' 'qcd_80to170_bctoe' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-300toInf_EMEnriched'  'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched' 'WJets' 'TTJets_aMC')
declare -a qcd_ee=('qcd_15to20_bctoe' 'qcd_170to250_bctoe' 'qcd_20to30_bctoe' 'qcd_250toinf_bctoe' 'qcd_30to80_bctoe' 'qcd_80to170_bctoe' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-300toInf_EMEnriched'  'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched')
declare -a qcd_eemm=('qcd_15to20_bctoe' 'qcd_170to250_bctoe' 'qcd_20to30_bctoe' 'qcd_250toinf_bctoe' 'qcd_30to80_bctoe' 'qcd_80to170_bctoe' 'QCD_Pt-120to170_EMEnriched' 'QCD_Pt-170to300_EMEnriched' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-300toInf_EMEnriched'  'QCD_Pt-30to50_EMEnriched' 'QCD_Pt-50to80_EMEnriched' 'QCD_Pt-80to120_EMEnriched' 'WJets' 'TTJets_aMC' 'QCD_Pt-1000toInf_MuEnriched' 'QCD_Pt-120to170_MuEnriched' 'QCD_Pt-15to20_MuEnriched' 'QCD_Pt-170to300_MuEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_Pt-300to470_MuEnriched' 'QCD_Pt-30to50_MuEnriched' 'QCD_Pt-470to600_MuEnriched' 'QCD_Pt-50to80_MuEnriched' 'QCD_Pt-600to800_MuEnriched' 'QCD_Pt-800to1000_MuEnriched' 'QCD_Pt-80to120_MuEnriched' 'TT_powheg')

declare -a hn_ll_ee2=('HN_t_ch_EmEm_100_official' 'HN_t_ch_EmEm_1100_official' 'HN_t_ch_EmEm_200_official' 'HN_t_ch_EmEm_500_official' 'HN_t_ch_EpEp_100_official' 'HN_t_ch_EpEp_1100_official' 'HN_t_ch_EpEp_200_official' 'HN_t_ch_EpEp_500_official' 'HNpair_ElEl_WR5000_Zp1000_HN100_official' 'HNpair_ElEl_WR5000_Zp1000_HN200_official' 'HNpair_ElEl_WR5000_Zp1000_HN300_official' 'HNpair_ElEl_WR5000_Zp1000_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN100_official' 'HNpair_ElEl_WR5000_Zp1500_HN200_official' 'HNpair_ElEl_WR5000_Zp1500_HN300_official' 'HNpair_ElEl_WR5000_Zp1500_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN500_official' 'HNpair_ElEl_WR5000_Zp1500_HN600_official' 'HNpair_ElEl_WR5000_Zp1500_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN100_official' 'HNpair_ElEl_WR5000_Zp2000_HN200_official' 'HNpair_ElEl_WR5000_Zp2000_HN300_official' 'HNpair_ElEl_WR5000_Zp2000_HN400_official' 'HNpair_ElEl_WR5000_Zp2000_HN500_official' 'HNpair_ElEl_WR5000_Zp2000_HN600_official' 'HNpair_ElEl_WR5000_Zp2000_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN800_official' 'HNpair_ElEl_WR5000_Zp2000_HN900_official' 'HNpair_ElEl_WR5000_Zp2500_HN1000_official' 'HNpair_ElEl_WR5000_Zp2500_HN100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1200_official' 'HNpair_ElEl_WR5000_Zp2500_HN200_official' 'HNpair_ElEl_WR5000_Zp2500_HN300_official' 'HNpair_ElEl_WR5000_Zp2500_HN400_official' 'HNpair_ElEl_WR5000_Zp2500_HN500_official' 'HNpair_ElEl_WR5000_Zp2500_HN600_official' 'HNpair_ElEl_WR5000_Zp2500_HN700_official' 'HNpair_ElEl_WR5000_Zp2500_HN800_official' 'HNpair_ElEl_WR5000_Zp2500_HN900_official' 'HNpair_ElEl_WR5000_Zp3000_HN1000_official' 'HNpair_ElEl_WR5000_Zp3000_HN100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1200_official' 'HNpair_ElEl_WR5000_Zp3000_HN1300_official' 'HNpair_ElEl_WR5000_Zp3000_HN1400_official' 'HNpair_ElEl_WR5000_Zp3000_HN200_official' 'HNpair_ElEl_WR5000_Zp3000_HN300_official' 'HNpair_ElEl_WR5000_Zp3000_HN400_official' 'HNpair_ElEl_WR5000_Zp3000_HN500_official' 'HNpair_ElEl_WR5000_Zp3000_HN600_official' 'HNpair_ElEl_WR5000_Zp3000_HN700_official' 'HNpair_ElEl_WR5000_Zp3000_HN800_official' 'HNpair_ElEl_WR5000_Zp3000_HN900_official' 'HNpair_ElEl_WR5000_Zp4000_HN1000_official' 'HNpair_ElEl_WR5000_Zp4000_HN100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1200_official' 'HNpair_ElEl_WR5000_Zp4000_HN1300_official' 'HNpair_ElEl_WR5000_Zp4000_HN1400_official' 'HNpair_ElEl_WR5000_Zp4000_HN1500_official' 'HNpair_ElEl_WR5000_Zp4000_HN1700_official' 'HNpair_ElEl_WR5000_Zp4000_HN1800_official' 'HNpair_ElEl_WR5000_Zp4000_HN1900_official' 'HNpair_ElEl_WR5000_Zp4000_HN200_official' 'HNpair_ElEl_WR5000_Zp4000_HN300_official' 'HNpair_ElEl_WR5000_Zp4000_HN400_official' 'HNpair_ElEl_WR5000_Zp4000_HN500_official' 'HNpair_ElEl_WR5000_Zp4000_HN600_official' 'HNpair_ElEl_WR5000_Zp4000_HN700_official' 'HNpair_ElEl_WR5000_Zp4000_HN800_official' 'HNpair_ElEl_WR5000_Zp4000_HN900_official' 'HNpair_ElEl_WR5000_Zp500_HN100_official' 'HNpair_ElEl_WR5000_Zp500_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN100_official' 'HNpair_ElEl_WR5000_Zp750_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN300_official'  'HNEmEm_100' 'HNEmEm_1100' 'HNEmEm_1500' 'HNEmEm_200' 'HNEmEm_40' 'HNEmEm_50' 'HNEmEm_500' 'HNEmEm_60' 'HNEmEp_100' 'HNEmEp_1100' 'HNEmEp_1500' 'HNEmEp_200' 'HNEmEp_40' 'HNEmEp_50' 'HNEmEp_500' 'HNEmEp_60' 'HNEpEm_100' 'HNEpEm_1100' 'HNEpEm_1500' 'HNEpEm_200' 'HNEpEm_40' 'HNEpEm_50' 'HNEpEm_500' 'HNEpEm_60' 'HNEpEp_100' 'HNEpEp_1100' 'HNEpEp_1500' 'HNEpEp_200' 'HNEpEp_40' 'HNEpEp_50' 'HNEpEp_500' 'HNEpEp_60' )

declare -a hn_ll_ee1=( 'HNEmEm_100' 'HNEmEm_1100' 'HNEmEm_1500' 'HNEmEm_200' 'HNEmEm_40' 'HNEmEm_50' 'HNEmEm_500' 'HNEmEm_60' 'HNEmEm_70' 'HNEmEm_80' 'HNEmEm_90' 'HNEmEm_125' 'HNEmEm_150' 'HNEmEm_300' 'HNEmEm_400' 'HNEmEm_600')

declare -a hn_ll_mm_1=( 'HNMumMum_100' 'HNMumMum_1000' 'HNMumMum_1100' 'HNMumMum_1200' 'HNMumMum_125' 'HNMumMum_1300' 'HNMumMum_1400' 'HNMumMum_150' 'HNMumMum_1500' 'HNMumMum_200' 'HNMumMum_250' 'HNMumMum_300' 'HNMumMum_40' 'HNMumMum_400' 'HNMumMum_50' 'HNMumMum_500' 'HNMumMum_60' 'HNMumMum_600' 'HNMumMum_70' 'HNMumMum_700' 'HNMumMum_80' 'HNMumMum_800' 'HNMumMum_90' 'HNMumMum_900' 'HNMumMum_Tchannel_1000' 'HNMoriondLL_Tchannel_MumMum_100' 'HNMoriondLL_Tchannel_MumMum_1100' 'HNMumMum_Tchannel_1200' 'HNMumMum_Tchannel_1300' 'HNMumMum_Tchannel_1400' 'HNMumMum_Tchannel_1500' 'HNMoriondLL_Tchannel_MumMum_200' 'HNMumMum_Tchannel_300' 'HNMumMum_Tchannel_400' 'HNMoriondLL_Tchannel_MumMum_500' 'HNMumMum_Tchannel_600' 'HNMumMum_Tchannel_700' 'HNMumMum_Tchannel_800' 'HNMumMum_Tchannel_900'  'HNMupMup_100' 'HNMupMup_1000' 'HNMupMup_1100' 'HNMupMup_1200' 'HNMupMup_125' 'HNMupMup_1300' 'HNMupMup_1400' 'HNMupMup_150' 'HNMupMup_1500' 'HNMupMup_200' 'HNMupMup_250' 'HNMupMup_300' 'HNMupMup_40' 'HNMupMup_400' 'HNMupMup_50' 'HNMupMup_500' 'HNMupMup_60' 'HNMupMup_600' 'HNMupMup_70' 'HNMupMup_700' 'HNMupMup_80' 'HNMupMup_800' 'HNMupMup_90' 'HNMupMup_900' 'HNMupMup_Tchannel_1000' 'HNMoriondLL_Tchannel_MupMup_100' 'HNMoriondLL_Tchannel_MupMup_1100' 'HNMupMup_Tchannel_1200' 'HNMupMup_Tchannel_1300' 'HNMupMup_Tchannel_1400' 'HNMupMup_Tchannel_1500' 'HNMoriondLL_Tchannel_MupMup_200' 'HNMupMup_Tchannel_300' 'HNMupMup_Tchannel_400' 'HNMoriondLL_Tchannel_MupMup_500' 'HNMupMup_Tchannel_600' 'HNMupMup_Tchannel_700' 'HNMupMup_Tchannel_800' 'HNMupMup_Tchannel_900' )




declare -a hn_pair_all=('HNpair_ElEl_WR5000_Zp1000_HN100_official' 'HNpair_ElEl_WR5000_Zp1000_HN200_official' 'HNpair_ElEl_WR5000_Zp1000_HN300_official' 'HNpair_ElEl_WR5000_Zp1000_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN100_official' 'HNpair_ElEl_WR5000_Zp1500_HN200_official' 'HNpair_ElEl_WR5000_Zp1500_HN300_official' 'HNpair_ElEl_WR5000_Zp1500_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN500_official' 'HNpair_ElEl_WR5000_Zp1500_HN600_official' 'HNpair_ElEl_WR5000_Zp1500_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN100_official' 'HNpair_ElEl_WR5000_Zp2000_HN200_official' 'HNpair_ElEl_WR5000_Zp2000_HN300_official' 'HNpair_ElEl_WR5000_Zp2000_HN400_official' 'HNpair_ElEl_WR5000_Zp2000_HN500_official' 'HNpair_ElEl_WR5000_Zp2000_HN600_official' 'HNpair_ElEl_WR5000_Zp2000_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN800_official' 'HNpair_ElEl_WR5000_Zp2000_HN900_official' 'HNpair_ElEl_WR5000_Zp2500_HN1000_official' 'HNpair_ElEl_WR5000_Zp2500_HN100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1200_official' 'HNpair_ElEl_WR5000_Zp2500_HN200_official' 'HNpair_ElEl_WR5000_Zp2500_HN300_official' 'HNpair_ElEl_WR5000_Zp2500_HN400_official' 'HNpair_ElEl_WR5000_Zp2500_HN500_official' 'HNpair_ElEl_WR5000_Zp2500_HN600_official' 'HNpair_ElEl_WR5000_Zp2500_HN700_official' 'HNpair_ElEl_WR5000_Zp2500_HN800_official' 'HNpair_ElEl_WR5000_Zp2500_HN900_official' 'HNpair_ElEl_WR5000_Zp3000_HN1000_official' 'HNpair_ElEl_WR5000_Zp3000_HN100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1200_official' 'HNpair_ElEl_WR5000_Zp3000_HN1300_official' 'HNpair_ElEl_WR5000_Zp3000_HN1400_official' 'HNpair_ElEl_WR5000_Zp3000_HN200_official' 'HNpair_ElEl_WR5000_Zp3000_HN300_official' 'HNpair_ElEl_WR5000_Zp3000_HN400_official' 'HNpair_ElEl_WR5000_Zp3000_HN500_official' 'HNpair_ElEl_WR5000_Zp3000_HN600_official' 'HNpair_ElEl_WR5000_Zp3000_HN700_official' 'HNpair_ElEl_WR5000_Zp3000_HN800_official' 'HNpair_ElEl_WR5000_Zp3000_HN900_official' 'HNpair_ElEl_WR5000_Zp4000_HN1000_official' 'HNpair_ElEl_WR5000_Zp4000_HN100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1200_official' 'HNpair_ElEl_WR5000_Zp4000_HN1300_official' 'HNpair_ElEl_WR5000_Zp4000_HN1400_official' 'HNpair_ElEl_WR5000_Zp4000_HN1500_official' 'HNpair_ElEl_WR5000_Zp4000_HN1700_official' 'HNpair_ElEl_WR5000_Zp4000_HN1800_official' 'HNpair_ElEl_WR5000_Zp4000_HN1900_official' 'HNpair_ElEl_WR5000_Zp4000_HN200_official' 'HNpair_ElEl_WR5000_Zp4000_HN300_official' 'HNpair_ElEl_WR5000_Zp4000_HN400_official' 'HNpair_ElEl_WR5000_Zp4000_HN500_official' 'HNpair_ElEl_WR5000_Zp4000_HN600_official' 'HNpair_ElEl_WR5000_Zp4000_HN700_official' 'HNpair_ElEl_WR5000_Zp4000_HN800_official' 'HNpair_ElEl_WR5000_Zp4000_HN900_official' 'HNpair_ElEl_WR5000_Zp500_HN100_official' 'HNpair_ElEl_WR5000_Zp500_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN100_official' 'HNpair_ElEl_WR5000_Zp750_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN300_official' 'HNpair_MuMu_WR5000_Zp1000_HN100_official' 'HNpair_MuMu_WR5000_Zp1000_HN200_official' 'HNpair_MuMu_WR5000_Zp1000_HN300_official' 'HNpair_MuMu_WR5000_Zp1000_HN400_official' 'HNpair_MuMu_WR5000_Zp1500_HN100_official' 'HNpair_MuMu_WR5000_Zp1500_HN200_official' 'HNpair_MuMu_WR5000_Zp1500_HN300_official' 'HNpair_MuMu_WR5000_Zp1500_HN400_official' 'HNpair_MuMu_WR5000_Zp1500_HN500_official' 'HNpair_MuMu_WR5000_Zp1500_HN600_official' 'HNpair_MuMu_WR5000_Zp1500_HN700_official' 'HNpair_MuMu_WR5000_Zp2000_HN100_official' 'HNpair_MuMu_WR5000_Zp2000_HN200_official' 'HNpair_MuMu_WR5000_Zp2000_HN300_official' 'HNpair_MuMu_WR5000_Zp2000_HN400_official' 'HNpair_MuMu_WR5000_Zp2000_HN500_official' 'HNpair_MuMu_WR5000_Zp2000_HN600_official' 'HNpair_MuMu_WR5000_Zp2000_HN700_official' 'HNpair_MuMu_WR5000_Zp2000_HN800_official' 'HNpair_MuMu_WR5000_Zp2000_HN900_official' 'HNpair_MuMu_WR5000_Zp2500_HN1000_official' 'HNpair_MuMu_WR5000_Zp2500_HN100_official' 'HNpair_MuMu_WR5000_Zp2500_HN1100_official' 'HNpair_MuMu_WR5000_Zp2500_HN1200_official' 'HNpair_MuMu_WR5000_Zp2500_HN200_official' 'HNpair_MuMu_WR5000_Zp2500_HN300_official' 'HNpair_MuMu_WR5000_Zp2500_HN400_official' 'HNpair_MuMu_WR5000_Zp2500_HN500_official' 'HNpair_MuMu_WR5000_Zp2500_HN600_official' 'HNpair_MuMu_WR5000_Zp2500_HN700_official' 'HNpair_MuMu_WR5000_Zp2500_HN800_official' 'HNpair_MuMu_WR5000_Zp2500_HN900_official' 'HNpair_MuMu_WR5000_Zp3000_HN1000_official' 'HNpair_MuMu_WR5000_Zp3000_HN100_official' 'HNpair_MuMu_WR5000_Zp3000_HN1100_official' 'HNpair_MuMu_WR5000_Zp3000_HN1200_official' 'HNpair_MuMu_WR5000_Zp3000_HN1300_official' 'HNpair_MuMu_WR5000_Zp3000_HN1400_official' 'HNpair_MuMu_WR5000_Zp3000_HN200_official' 'HNpair_MuMu_WR5000_Zp3000_HN300_official' 'HNpair_MuMu_WR5000_Zp3000_HN400_official' 'HNpair_MuMu_WR5000_Zp3000_HN500_official' 'HNpair_MuMu_WR5000_Zp3000_HN600_official' 'HNpair_MuMu_WR5000_Zp3000_HN700_official' 'HNpair_MuMu_WR5000_Zp3000_HN800_official' 'HNpair_MuMu_WR5000_Zp3000_HN900_official' 'HNpair_MuMu_WR5000_Zp4000_HN1000_official' 'HNpair_MuMu_WR5000_Zp4000_HN100_official' 'HNpair_MuMu_WR5000_Zp4000_HN1100_official' 'HNpair_MuMu_WR5000_Zp4000_HN1200_official' 'HNpair_MuMu_WR5000_Zp4000_HN1300_official' 'HNpair_MuMu_WR5000_Zp4000_HN1400_official' 'HNpair_MuMu_WR5000_Zp4000_HN1500_official' 'HNpair_MuMu_WR5000_Zp4000_HN1600_official' 'HNpair_MuMu_WR5000_Zp4000_HN1700_official' 'HNpair_MuMu_WR5000_Zp4000_HN1800_official' 'HNpair_MuMu_WR5000_Zp4000_HN1900_official' 'HNpair_MuMu_WR5000_Zp4000_HN200_official' 'HNpair_MuMu_WR5000_Zp4000_HN300_official' 'HNpair_MuMu_WR5000_Zp4000_HN400_official' 'HNpair_MuMu_WR5000_Zp4000_HN500_official' 'HNpair_MuMu_WR5000_Zp4000_HN600_official' 'HNpair_MuMu_WR5000_Zp4000_HN700_official' 'HNpair_MuMu_WR5000_Zp4000_HN800_official' 'HNpair_MuMu_WR5000_Zp4000_HN900_official' 'HNpair_MuMu_WR5000_Zp500_HN100_official' 'HNpair_MuMu_WR5000_Zp500_HN200_official' 'HNpair_MuMu_WR5000_Zp750_HN100_official' 'HNpair_MuMu_WR5000_Zp750_HN200_official' 'HNpair_MuMu_WR5000_Zp750_HN300_official')


declare -a hn_ll_ee_john=('HNEmEm_100' 'HNEmEm_1000' 'HNEmEm_1100' 'HNEmEm_1200' 'HNEmEm_125' 'HNEmEm_1300' 'HNEmEm_1400' 'HNEmEm_150' 'HNEmEm_1500' 'HNEmEm_200' 'HNEmEm_300' 'HNEmEm_40' 'HNEmEm_400' 'HNEmEm_50' 'HNEmEm_500' 'HNEmEm_60'  'HNEmEm_600' 'HNEmEm_70' 'HNEmEm_700' 'HNEmEm_80' 'HNEmEm_800' 'HNEmEm_90' 'HNEmEm_900' 'HNEmEm_Tchannel_1000' 'HNEmEm_Tchannel_1200' 'HNEmEm_Tchannel_1300' 'HNEmEm_Tchannel_1400' 'HNEmEm_Tchannel_1500' 'HNEmEm_Tchannel_300' 'HNEmEm_Tchannel_400' 'HNEmEm_Tchannel_600' 'HNEmEm_Tchannel_700' 'HNEmEm_Tchannel_800' 'HNEmEm_Tchannel_900' 'HNEmEp_100' 'HNEmEp_1000' 'HNEmEp_1100' 'HNEmEp_1200' 'HNEmEp_125' 'HNEmEp_1300' 'HNEmEp_1400' 'HNEmEp_150' 'HNEmEp_1500' 'HNEmEp_200' 'HNEmEp_250' 'HNEmEp_300' 'HNEmEp_40' 'HNEmEp_400' 'HNEmEp_50' 'HNEmEp_500' 'HNEmEp_60' 'HNEmEp_600' 'HNEmEp_70' 'HNEmEp_700' 'HNEmEp_80' 'HNEmEp_800' 'HNEmEp_90' 'HNEmEp_900' 'HNEpEm_100' 'HNEpEm_1000' 'HNEpEm_1100' 'HNEpEm_1200' 'HNEpEm_125' 'HNEpEm_1300' 'HNEpEm_1400' 'HNEpEm_150' 'HNEpEm_1500' 'HNEpEm_200' 'HNEpEm_250' 'HNEpEm_300' 'HNEpEm_40' 'HNEpEm_400' 'HNEpEm_50' 'HNEpEm_500' 'HNEpEm_60' 'HNEpEm_600' 'HNEpEm_70' 'HNEpEm_700' 'HNEpEm_80' 'HNEpEm_800' 'HNEpEm_90' 'HNEpEm_900' 'HNEpEp_100' 'HNEpEp_1000' 'HNEpEp_1100' 'HNEpEp_1200' 'HNEpEp_125' 'HNEpEp_1300' 'HNEpEp_1400' 'HNEpEp_150' 'HNEpEp_1500' 'HNEpEp_200' 'HNEpEp_250' 'HNEpEp_300' 'HNEpEp_40' 'HNEpEp_400' 'HNEpEp_50' 'HNEpEp_500' 'HNEpEp_60' 'HNEpEp_600' 'HNEpEp_70' 'HNEpEp_700' 'HNEpEp_80' 'HNEpEp_800' 'HNEpEp_90' 'HNEpEp_900' 'HNEpEp_Tchannel_1000' 'HNEpEp_Tchannel_1200' 'HNEpEp_Tchannel_1300' 'HNEpEp_Tchannel_1400' 'HNEpEp_Tchannel_1500' 'HNEpEp_Tchannel_300' 'HNEpEp_Tchannel_400' 'HNEpEp_Tchannel_600' 'HNEpEp_Tchannel_700' 'HNEpEp_Tchannel_800' 'HNEpEp_Tchannel_900' )

declare -a hn_pair_ee=('HNpair_ElEl_WR5000_Zp1000_HN100_official' 'HNpair_ElEl_WR5000_Zp1000_HN200_official' 'HNpair_ElEl_WR5000_Zp1000_HN300_official' 'HNpair_ElEl_WR5000_Zp1000_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN100_official' 'HNpair_ElEl_WR5000_Zp1500_HN200_official' 'HNpair_ElEl_WR5000_Zp1500_HN300_official' 'HNpair_ElEl_WR5000_Zp1500_HN400_official' 'HNpair_ElEl_WR5000_Zp1500_HN500_official' 'HNpair_ElEl_WR5000_Zp1500_HN600_official' 'HNpair_ElEl_WR5000_Zp1500_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN100_official' 'HNpair_ElEl_WR5000_Zp2000_HN200_official' 'HNpair_ElEl_WR5000_Zp2000_HN300_official' 'HNpair_ElEl_WR5000_Zp2000_HN400_official' 'HNpair_ElEl_WR5000_Zp2000_HN500_official' 'HNpair_ElEl_WR5000_Zp2000_HN600_official' 'HNpair_ElEl_WR5000_Zp2000_HN700_official' 'HNpair_ElEl_WR5000_Zp2000_HN800_official' 'HNpair_ElEl_WR5000_Zp2000_HN900_official' 'HNpair_ElEl_WR5000_Zp2500_HN1000_official' 'HNpair_ElEl_WR5000_Zp2500_HN100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1100_official' 'HNpair_ElEl_WR5000_Zp2500_HN1200_official' 'HNpair_ElEl_WR5000_Zp2500_HN200_official' 'HNpair_ElEl_WR5000_Zp2500_HN300_official' 'HNpair_ElEl_WR5000_Zp2500_HN400_official' 'HNpair_ElEl_WR5000_Zp2500_HN500_official' 'HNpair_ElEl_WR5000_Zp2500_HN600_official' 'HNpair_ElEl_WR5000_Zp2500_HN700_official' 'HNpair_ElEl_WR5000_Zp2500_HN800_official' 'HNpair_ElEl_WR5000_Zp2500_HN900_official' 'HNpair_ElEl_WR5000_Zp3000_HN1000_official' 'HNpair_ElEl_WR5000_Zp3000_HN100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1100_official' 'HNpair_ElEl_WR5000_Zp3000_HN1200_official' 'HNpair_ElEl_WR5000_Zp3000_HN1300_official' 'HNpair_ElEl_WR5000_Zp3000_HN1400_official' 'HNpair_ElEl_WR5000_Zp3000_HN200_official' 'HNpair_ElEl_WR5000_Zp3000_HN300_official' 'HNpair_ElEl_WR5000_Zp3000_HN400_official' 'HNpair_ElEl_WR5000_Zp3000_HN500_official' 'HNpair_ElEl_WR5000_Zp3000_HN600_official' 'HNpair_ElEl_WR5000_Zp3000_HN700_official' 'HNpair_ElEl_WR5000_Zp3000_HN800_official' 'HNpair_ElEl_WR5000_Zp3000_HN900_official' 'HNpair_ElEl_WR5000_Zp4000_HN1000_official' 'HNpair_ElEl_WR5000_Zp4000_HN100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1100_official' 'HNpair_ElEl_WR5000_Zp4000_HN1200_official' 'HNpair_ElEl_WR5000_Zp4000_HN1300_official' 'HNpair_ElEl_WR5000_Zp4000_HN1400_official' 'HNpair_ElEl_WR5000_Zp4000_HN1500_official' 'HNpair_ElEl_WR5000_Zp4000_HN1700_official' 'HNpair_ElEl_WR5000_Zp4000_HN1800_official' 'HNpair_ElEl_WR5000_Zp4000_HN1900_official' 'HNpair_ElEl_WR5000_Zp4000_HN200_official' 'HNpair_ElEl_WR5000_Zp4000_HN300_official' 'HNpair_ElEl_WR5000_Zp4000_HN400_official' 'HNpair_ElEl_WR5000_Zp4000_HN500_official' 'HNpair_ElEl_WR5000_Zp4000_HN600_official' 'HNpair_ElEl_WR5000_Zp4000_HN700_official' 'HNpair_ElEl_WR5000_Zp4000_HN800_official' 'HNpair_ElEl_WR5000_Zp4000_HN900_official' 'HNpair_ElEl_WR5000_Zp500_HN100_official' 'HNpair_ElEl_WR5000_Zp500_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN100_official' 'HNpair_ElEl_WR5000_Zp750_HN200_official' 'HNpair_ElEl_WR5000_Zp750_HN300_official')
declare -a hn_pair_mm=('HNpair_MuMu_WR5000_Zp1000_HN100_official' 'HNpair_MuMu_WR5000_Zp1000_HN200_official' 'HNpair_MuMu_WR5000_Zp1000_HN300_official' 'HNpair_MuMu_WR5000_Zp1000_HN400_official' 'HNpair_MuMu_WR5000_Zp1500_HN100_official' 'HNpair_MuMu_WR5000_Zp1500_HN200_official' 'HNpair_MuMu_WR5000_Zp1500_HN300_official' 'HNpair_MuMu_WR5000_Zp1500_HN400_official' 'HNpair_MuMu_WR5000_Zp1500_HN500_official' 'HNpair_MuMu_WR5000_Zp1500_HN600_official' 'HNpair_MuMu_WR5000_Zp1500_HN700_official' 'HNpair_MuMu_WR5000_Zp2000_HN100_official' 'HNpair_MuMu_WR5000_Zp2000_HN200_official' 'HNpair_MuMu_WR5000_Zp2000_HN300_official' 'HNpair_MuMu_WR5000_Zp2000_HN400_official' 'HNpair_MuMu_WR5000_Zp2000_HN500_official' 'HNpair_MuMu_WR5000_Zp2000_HN600_official' 'HNpair_MuMu_WR5000_Zp2000_HN700_official' 'HNpair_MuMu_WR5000_Zp2000_HN800_official' 'HNpair_MuMu_WR5000_Zp2000_HN900_official' 'HNpair_MuMu_WR5000_Zp2500_HN1000_official' 'HNpair_MuMu_WR5000_Zp2500_HN100_official' 'HNpair_MuMu_WR5000_Zp2500_HN1100_official' 'HNpair_MuMu_WR5000_Zp2500_HN1200_official' 'HNpair_MuMu_WR5000_Zp2500_HN200_official' 'HNpair_MuMu_WR5000_Zp2500_HN300_official' 'HNpair_MuMu_WR5000_Zp2500_HN400_official' 'HNpair_MuMu_WR5000_Zp2500_HN500_official' 'HNpair_MuMu_WR5000_Zp2500_HN600_official' 'HNpair_MuMu_WR5000_Zp2500_HN700_official' 'HNpair_MuMu_WR5000_Zp2500_HN800_official' 'HNpair_MuMu_WR5000_Zp2500_HN900_official' 'HNpair_MuMu_WR5000_Zp3000_HN1000_official' 'HNpair_MuMu_WR5000_Zp3000_HN100_official' 'HNpair_MuMu_WR5000_Zp3000_HN1100_official' 'HNpair_MuMu_WR5000_Zp3000_HN1200_official' 'HNpair_MuMu_WR5000_Zp3000_HN1300_official' 'HNpair_MuMu_WR5000_Zp3000_HN1400_official' 'HNpair_MuMu_WR5000_Zp3000_HN200_official' 'HNpair_MuMu_WR5000_Zp3000_HN300_official' 'HNpair_MuMu_WR5000_Zp3000_HN400_official' 'HNpair_MuMu_WR5000_Zp3000_HN500_official' 'HNpair_MuMu_WR5000_Zp3000_HN600_official' 'HNpair_MuMu_WR5000_Zp3000_HN700_official' 'HNpair_MuMu_WR5000_Zp3000_HN800_official' 'HNpair_MuMu_WR5000_Zp3000_HN900_official' 'HNpair_MuMu_WR5000_Zp4000_HN1000_official' 'HNpair_MuMu_WR5000_Zp4000_HN100_official' 'HNpair_MuMu_WR5000_Zp4000_HN1100_official' 'HNpair_MuMu_WR5000_Zp4000_HN1200_official' 'HNpair_MuMu_WR5000_Zp4000_HN1300_official' 'HNpair_MuMu_WR5000_Zp4000_HN1400_official' 'HNpair_MuMu_WR5000_Zp4000_HN1500_official' 'HNpair_MuMu_WR5000_Zp4000_HN1600_official' 'HNpair_MuMu_WR5000_Zp4000_HN1700_official' 'HNpair_MuMu_WR5000_Zp4000_HN1800_official' 'HNpair_MuMu_WR5000_Zp4000_HN1900_official' 'HNpair_MuMu_WR5000_Zp4000_HN200_official' 'HNpair_MuMu_WR5000_Zp4000_HN300_official' 'HNpair_MuMu_WR5000_Zp4000_HN400_official' 'HNpair_MuMu_WR5000_Zp4000_HN500_official' 'HNpair_MuMu_WR5000_Zp4000_HN600_official' 'HNpair_MuMu_WR5000_Zp4000_HN700_official' 'HNpair_MuMu_WR5000_Zp4000_HN800_official' 'HNpair_MuMu_WR5000_Zp4000_HN900_official' 'HNpair_MuMu_WR5000_Zp500_HN100_official' 'HNpair_MuMu_WR5000_Zp500_HN200_official' 'HNpair_MuMu_WR5000_Zp750_HN100_official' 'HNpair_MuMu_WR5000_Zp750_HN200_official' 'HNpair_MuMu_WR5000_Zp750_HN300_official')
declare -a example=("WW" "WZ" "HNMupMup_100" "DYJets")

declare -a tmplist=('ggZZto4e' 'ttZToLL_M-1to10' 'ttWToLNu' 'ttZToLL_M-10' 'WZG' 'WWG')

declare -a tmpall_mc=('TTJets_aMC' 'LowStat_WJets' 'WW'  'WZ' 'ZZ' 'LowStat_DYJets' 'DYJets_10to50' )

declare -a hn=('DYJets_10to50'  'DYJets' 'TT_powheg' 'WJets' 'WGtoLNuG'  'ZGto2LG' 'qcd_15to20_bctoe' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-20to30_MuEnriched')
declare -a hn_fake=('DYJets' 'TT_powheg' 'WJets' 'qcd_15to20_bctoe' 'QCD_Pt-20to30_EMEnriched' 'QCD_Pt-20to30_MuEnriched' 'QCD_DoubleEMEnriched_30-40_mgg80toinf' 'QCD_DoubleEMEnriched_30-inf_mgg40to80' 'QCD_DoubleEMEnriched_40-inf_mgg80toinf')

declare -a pu_dilepton_list=('DYJets_10to50' 'DYJets' 'WJets' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg' 'DYJets_MG_10to50' 'DYJets_MG' 'TTJets_aMC' )

declare -a hntmp=('TTTT' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds'  )

declare -a hn_eetmp=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'WWTo2L2Nu' 'WWToLNuQQ' 'QCD_DoubleEMEnriched_30-40_mgg80toinf' 'QCD_DoubleEMEnriched_30-inf_mgg40to80' 'QCD_DoubleEMEnriched_40-inf_mgg80toinf' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2e2nu' 'ggZZto2e2tau'  'ggZZto4e' 'ggWWto2L2Nu' 'ww_ds'  )
declare -a hn_mm_all=('DYJets_10to50' 'DYJets' 'WJets' 'WpWpEWK' 'WpWpQCD' 'TT_powheg'  'SingleTop_s' 'SingleTbar_t' 'SingleTop_t'  'SingleTbar_tW' 'SingleTop_tW' 'WWW' 'ttW' 'ttZ' 'ttH_nonbb' 'ttH_bb' 'ZZZ' 'WZZ'  'VBF_HToMuMu' 'WGtoLNuG'  'ZGto2LG' 'WZTo3LNu_powheg' 'ZZTo4L_powheg'  'WWTo2L2Nu' 'WWToLNuQQ' 'TG' 'TTG' 'ttWToLNu' 'ttZToLL_M-1to10' 'ttZToLL_M-10'  'tZq' 'ggHtoWW' 'ggHtoZZ' 'WWG' 'WZG' 'WZto2L2Q_amcatnlo' 'ZZTo2L2Nu_Powheg' 'ZZTo2L2Q_Powheg' 'ggZZto2e2mu' 'ggZZto2mu2nu' 'ggZZto2mu2tau'  'ggZZto4mu' 'ggWWto2L2Nu' 'ww_ds'  )


declare -a sktmp=('WJets' 'TT_powheg')
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
declare -a new_list=('DYJets' 'LowStatDYJets' )
