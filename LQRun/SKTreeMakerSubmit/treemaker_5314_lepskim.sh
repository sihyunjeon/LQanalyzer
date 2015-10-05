#!/bin/sh
### sets all configurable variables to defaul values

###### SET WHAT JOBS TO RUN
runMC=false

runDoubleMuon=true
runDoubleElectron=false
runSingleMuon=false
runElectronMuon=false

if [[ $1  == "ALL" ]]; 
then
    runMC=true
    runDoubleMuon=true
    runDoubleElectron=true
    runSingleMuon=true
    runElectronMuon=true
fi


if [[ $1  == "MC" ]];
then
    runMC=true
    runDoubleMuon=false
    runDoubleElectron=false
    runSingleMuon=false
fi

if [[ $1  == "DATA" ]];
then
    runMC=false
    runDoubleMuon=true
    runDoubleElectron=true
    runSingleMuon=true
    runElectronMuon=true
fi    



if [[ $runMC  == "true" ]]; 
then
    source functions.sh
    cycle="SKTreeMaker"
    #### JOB CONFIGURATION
    njobs=30
    data_lumi="C"
    loglevel="INFO"
    logstep=1000
    
    declare -a input_samples=("WZ" "ZZ" "WW" "WJets" "DY10to50" "DY50plus" "TTJets_MG5"  "QCD_mu20to30" "QCD_mu30to50" "QCD_mu50to80" "QCD_mu80to120" "QCD_mu120to170" "QCD_mu170to300" "QCD_mu300to470" "QCD_mu470to600" "QCD_mu600to800" "QCD_mu800to1000" "QCD_mu1000toINF" "QCD_em20to30" "QCD_em30to50" "QCD_em50to80" "QCD_em80to120" "QCD_em120to170" "QCD_em170to300" "QCD_em300toINF" "singletop_tbar" "singletop_t" "singletop_tbarW" "singletop_tW" "ttHtobb" "ttHnobb" "ttWJetsToLNu" "ttWJetsToQQ" "ttZToLLNuNu" "ttZToQQ")

    source submit.sh
fi    


################ DOUBLEMUON DATA
### submit this configured job (uses bin/submit.sh)
if [[ $runSingleMuon  == "true" ]];
then
    source functions.sh
    cycle="SKTreeMaker"

    njobs=30
    data_lumi="C" ### This should == input_samples
    loglevel="INFO"
    logstep=1000
    
    stream="singlemuon"
    declare -a input_samples=("C" )

    source submit.sh
fi


################ DOUBLEELECTRON DATA
### submit this configured job (uses bin/submit.sh)
if [[ $runDoubleElectron  == "true" ]];
then
    source functions.sh
    cycle="SKTreeMaker"
    njobs=30
    data_lumi="C"
    loglevel="INFO"
    logstep=1000
    
    stream="egamma"
    declare -a input_samples=("C")

    source submit.sh
fi    

################ DOUBLEMUON DATA
### submit this configured job (uses bin/submit.sh)
if [[ $runDoubleMuon  == "true" ]];
then
    source functions.sh
    cycle="SKTreeMaker"
    njobs=30
    data_lumi="C" ### This should == input_samples  
    loglevel="INFO"
    logstep=1000
    
    stream="muon"
    declare -a input_samples=("C" )

    source submit.sh
fi

if [[ $runElectronMuon  == "true" ]];
then
    source functions.sh
    cycle="SKTreeMaker"

    njobs=30
    data_lumi="C" ### This should == input_samples
    loglevel="INFO"
    logstep=1000


    stream="emu"
    declare -a input_samples=("C" )

    source submit.sh
fi


