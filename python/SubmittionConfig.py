############################################################    
### configure submisstion of CATANALYZER Jobs                                                                                                        #################################################################    
import os, getpass, sys,ROOT,time,curses,datetime
from functions import *
from datetime import timedelta
from optparse import OptionParser
from QuickHistCheck import *


path_jobpre="/data1/"
if "tamsa2.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    path_jobpre="/data2/"

an_jonpre="/data2/"
if "tamsa2.snu.ac.kr" in str(os.getenv("HOSTNAME")):
    an_jonpre="/data4/"




large_memory_check=800
large_file_size=100
file_increase_warning=2.
time_increase_warning=2.
######## set 10 process time for which email is sent
email_time_limit=600.

def   GetMonth(imonth):
    if imonth == 1:
        return "Jan"
    if imonth == 2:
        return "Feb"
    if imonth == 3:
        return "Mar"
    if imonth == 4:
        return "Apr"
    if imonth == 5:
        return "May"
    if imonth == 6:
        return "June"
    if imonth == 7:
        return "July"
    if imonth == 8:
        return "Aug"
    if imonth == 9:
        return "Sep"
    if imonth == 10:
        return "Oct"
    if imonth == 11:
        return "Nov"
    if imonth == 12:
        return "Dec"
    else:
        return "Dec"

def   MergeData(defrunnp,defruncf,defdata_lumi, defFinaloutputdir,  defcatversion, defuseskim, defcycle, defchannel):

    defoutput_file_skim_tag=defchannel
    if defuseskim == "FLATCAT":
        defoutput_file_skim_tag=defoutput_file_skim_tag+"_cat_"+defcatversion
    if defuseskim == "SKTree_LeptonSkim" :
        defoutput_file_skim_tag="SK"+defoutput_file_skim_tag+"_cat_"+defcatversion
    if defuseskim == "SKTree_DiLepSkim" :
        defoutput_file_skim_tag="SK"+defoutput_file_skim_tag+"_dilep_cat_"+defcatversion

    if defrunnp == "True":
        defoutput_file_skim_tag=defchannel
        foutname="nonprompt"
        if defuseskim == "FLATCAT":
            defoutput_file_skim_tag=defoutput_file_skim_tag+"_cat_"+defcatversion
            foutname=foutname+"_cat_"+defcatversion
        if defuseskim == "SKTree_LeptonSkim" :
            defoutput_file_skim_tag="SK"+defoutput_file_skim_tag+"_cat_"+defcatversion
            foutname="SK"+foutname+"_cat_"+defcatversion
        if defuseskim == "SKTree_DiLepSkim" :
            defoutput_file_skim_tag="SK"+defoutput_file_skim_tag+"_dilep_cat_"+defcatversion
            foutname="SK"+foutname+"_dilep_cat_"+defcatversion
        
            
        defFinaloutputdirMC=""

        if "/Fake/" in defFinaloutputdir:
            defFinaloutputdirMC=defFinaloutputdir
            defFinaloutputdirMC=defFinaloutputdirMC.replace("Fake/","")

        if defdata_lumi == "ALL" or defdata_lumi==os.getenv("catdatatag"):
            if not  "SKTreeMaker" in cycle:
                os.system("source hadd.sh " + defFinaloutputdir + " "+defcycle+"_"+defoutput_file_skim_tag+".root "+defFinaloutputdir+"/"+defcycle+"'*'"+defoutput_file_skim_tag+"'*'")
                os.system("mv "  + defFinaloutputdir+ "/"+ defcycle+"_"+defoutput_file_skim_tag+".root " + defFinaloutputdirMC+ "/"+defcycle+ "_"+defchannel+"_"+foutname+".root")


    elif defruncf == "True":
        print ""
    else:
        defoutput_file_skim_tag=defchannel
        if defuseskim == "FLATCAT":
            defoutput_file_skim_tag=defoutput_file_skim_tag+"_cat_"+defcatversion
        if defuseskim == "SKTree_LeptonSkim" :
            defoutput_file_skim_tag="SK"+defoutput_file_skim_tag+"_cat_"+defcatversion
        if defuseskim == "SKTree_DiLepSkim" :
            defoutput_file_skim_tag="SK"+defoutput_file_skim_tag+"_dilep_cat_"+defcatversion

        defFinaloutputdirMC=""

        if "/Data/" in defFinaloutputdir:
            defFinaloutputdirMC=defFinaloutputdir
            defFinaloutputdirMC=defFinaloutputdirMC.replace("Data/","")
        else:
            return 
        
        
        output_datafile=defFinaloutputdirMC+"/"+defcycle+"_data_cat_"+ defcatversion+".root"
        if defdata_lumi == "ALL" or defdata_lumi==os.getenv("catdatatag"):
            if not "SKTreeMaker" in cycle:

                os.system("source hadd.sh " + defFinaloutputdir + " "+defcycle+"_data_cat_"+defcatversion+".root "+defFinaloutputdir+"/"+defcycle+"'*'"+defoutput_file_skim_tag+"'*'")
                os.system("mv "  + defFinaloutputdir+ "/"+defcycle+"_data_cat_"+defcatversion+".root  " + defFinaloutputdirMC+ "/"+defcycle+"_data_" + defchannel+"_cat_"+defcatversion+".root")

def UpdateOutput(outputlist,outputlist_path):
    out_file = open(outputlist_path,"w")
    for x in outputlist:
        out_file.write(x+"\n")
    out_file.close()
    
def GetListOld():
    
    ### List of tags with stst logging but no memeory values in the Master File
    clist=[]
    clist.append("v8-0-1.7")
    clist.append("v8-0-1.6")
    return clist

def GetListNew():

    ### change in format (add memory of jobs since v8-0-1.8)
    ### This new format allows job to output WARNING if jobs are using more memory than previously.
    ### This should help fix any possible memory leaks when the happen
    clist=[]
    clist.append("v8-0-1.8")
    clist.append("v8-0-2")### this will take care of all 802 tags
    return clist

def GetList():
    
    clist=GetListNew() + GetListOld()
    #for ct in GetListNew():
    return clist

def NewForat(ct):
    
    clist = GetListNew()
    for ctag in clist:
        if ct == ctag:
            return True
        if ctag in ct:
            return True
    return False


def GetNFiles( deftagger,defsample,defcycle,defskim):

    if not os.path.exists(path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/MasterFile_"+ os.getenv("CATVERSION")+".txt"):
        return -999

    nit=2
    avg_time=-999
    checkdate = datetime.datetime.now()
    diff = datetime.timedelta(days=30)
    checkdate=checkdate+diff
    get_nfiles=0
    while avg_time < 0:
        if nit < 0:
            return 0.
        nit =nit-1
        #### get previous month                                                                                                                                                 
        checkdate=checkdate-diff
        month=checkdate.strftime("%m")
        day=checkdate.strftime("%d")
        year=checkdate.strftime("%y")

        file_jobsummary="/data1//LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/JobSummary"+str(GetMonth(int(month)))+"_20"+year+".txt"
        read_file_jobsummary = open(file_jobsummary,"r")
    

        for line in read_file_jobsummary:
            if os.getenv("USER") in line:
                if defsample in line and defskim in line:
                    splitline = line.split()
                    nthsplit=0
                    for s in splitline:
                        if nthsplit== 12:
                            get_nfiles=s
                        nthsplit=nthsplit+1
        read_file_jobsummary.close()

        if get_nfiles < 1:
            continue
        else:
            return get_nfiles

    return get_nfiles


def GetAverageTime( deftagger,defsample,defcycle,defskim):
    if not os.path.exists(path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/MasterFile_"+ os.getenv("CATVERSION")+".txt"):
        return 1.
    

    nit=2
    avg_time=-999
    checkdate = datetime.datetime.now()
    diff = datetime.timedelta(days=30)
    checkdate=checkdate+diff
    gettime_nfiles=0
    gettime_njobs=0
    gettime_jobtime=0
    
    while avg_time < 0:
        if nit < 0:
            return 0.
        nit =nit-1
        #### get previous month
        checkdate=checkdate-diff
        month=checkdate.strftime("%m")
        day=checkdate.strftime("%d")
        year=checkdate.strftime("%y")

        file_jobsummary="/data1//LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/JobSummary"+str(GetMonth(int(month)))+"_20"+year+".txt"
        read_file_jobsummary = open(file_jobsummary,"r")
        for line in read_file_jobsummary:
            if os.getenv("USER") in line:
                if defsample in line and defskim in line:
                    splitline = line.split()
                    nthsplit=0
                    for s in splitline:
                        if nthsplit==10:
                            gettime_njobs=s
                        if nthsplit==12:
                            gettime_nfiles=s
                        if nthsplit==14:
                            gettime_jobtime=s
                        nthsplit=nthsplit+1
        read_file_jobsummary.close()
        if gettime_jobtime < 1.:
            continue
        if gettime_nfiles < 1:
            continue
        if gettime_njobs < 1:
            continue
        gettime_jobtime = float(gettime_jobtime) / float(gettime_nfiles) 
        gettime_jobtime = float(gettime_jobtime) * float(gettime_njobs)
        avg_time=float(gettime_jobtime)
        return avg_time
    return avg_time
    

def DetermineNjobs(longestjobtime, ncores_job, deftagger,defsample,defcycle,defskim):

    if not os.path.exists(path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/MasterFile_"+ os.getenv("CATVERSION")+".txt"):
        return 1000

    if ncores_job == 1:
        return 1
    
    expectedjobtime=GetAverageTime(deftagger, defsample, defcycle,defskim)
    expectedjobnfiles=GetNFiles(deftagger, defsample, defcycle,defskim)
    if not longestjobtime == expectedjobtime:
        for i in range(2, 20):
            if float(expectedjobtime*float(expectedjobnfiles)/float(i)) < longestjobtime:
                return i
    

    path_clust_check_njobs=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(deftagger)+ "/clustercheck.txt"
    os.system("qstat -u '*'  > " +  path_clust_check_njobs)
    file_clust_check_njobs=open(path_clust_check_njobs ,"r")
    nrunning=0.
    for cline in file_clust_check_njobs:
        nrunning = nrunning + 1.
    file_clust_check_njobs.close()

    pre_job_time=GetTime(defsample,defcycle,deftagger, defskim)
    if pre_job_time == "None":
        pre_job_time=100.
    njobs_max=50
    if float(pre_job_time) < 100.:
        njobs_max=50
    elif float(pre_job_time) < 200.:
        njobs_max=25
    elif float(pre_job_time) < 300.:
        njobs_max=20
    elif float(pre_job_time) < 500.:
        njobs_max=15
    elif float(pre_job_time) < 1000.:
        njobs_max=10
    elif float(pre_job_time) < 1500.:
        njobs_max=5
    else:
        njobs_max=5

    if nrunning < 50:
        return njobs_max
    else:
        return njobs_max
        if njobs_max> 40:
            return njobs_max - 40
        else:
            return njobs_max

def CheckJobHistory(info_type, defsample, defcycle, tagger,defskim):

    if "FLATCAT" in defskim:
        defsample +="_lepton"        
    if "SKTree_LeptonSkim" in defskim:
        defsample +="_lepton"
    if "SKTree_DiLepSkim"in defskim:
        defsample +="_dilepton"


    
    jobline=""
    list_tags=GetList()
    
    itag=-1
    while not jobline:
        itag=itag+1
        cattag=list_tags[itag]
        newformat=NewForat(cattag)

        #info_file_master = path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/MasterFile_"+ os.getenv("CATVERSION")+".txt"
        info_file= path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + os.getenv("USER") + "/PreMasterFile_"+os.getenv("CATVERSION")+ str(tagger)+".txt"

        #os.system("cp " + info_file_master + " " + info_file)
        read_info_file = open(info_file,"r")
        isuser=False
        iscycle=False
        jobline=""
        for infoline in read_info_file:
            if "USER" in infoline:
                if  os.getenv("USER") in infoline:
                    isuser=True
                else:
                    isuser=False
            if isuser: 
                if "Code:" in infoline and defcycle in infoline:
                    iscycle=True
                if "Code:" in infoline and not defcycle in infoline:
                    iscycle=False
                if iscycle:
                    if defsample in infoline:
                        if cattag in infoline:
                            jobline=infoline
        read_info_file.close()            
        if itag == (len(list_tags) - 1):
            if not jobline:
                jobline="None"
    if jobline == "None":
        return "None"

    sinfoline = jobline.split()
    if newformat:
        if len(sinfoline) == 7:
            if info_type == "MemoryV":
                return sinfoline[5]
            if info_type == "MemoryP":
                return sinfoline[6]
            if info_type == "Time":
                return sinfoline[3]
            if info_type == "FileSize":
                return sinfoline[4]
        
    else:
        if len(sinfoline) == 5:
            if info_type == "MemoryV":
                return "None"
            if info_type == "MemoryP":
                return "None"
            if info_type == "Time":
                return sinfoline[3]
            if info_type == "FileSize":
                return sinfoline[4]

    return "None"        

def GetVirtualMemoryUsage(defsample, defcycle,tagger,defskim):
    return CheckJobHistory("MemoryV", defsample, defcycle, tagger,defskim)

def GetPhysicalMemoryUsage(defsample, defcycle, tagger,defskim):
    return CheckJobHistory("MemoryP", defsample, defcycle, tagger,defskim)

def GetTime(defsample, defcycle, tagger,defskim):
    return CheckJobHistory("Time", defsample, defcycle, tagger,defskim)

def GetFileSize(defsample, defcycle,tagger,defskim):
    return CheckJobHistory("FileSize", defsample, defcycle,tagger,defskim)
                                
                
def SendEmail(jobsummary, deftagger, e_subject, email_user, sendplots, plotlist):

    if not  os.getenv("USER") == "jalmond":
        if "jalmond" in email_user:
            print "Email could not be set since email address is not set correctly in bin/catconfig."
            return 

    plotstring =""
    for x in plotlist:
        plotstring = plotstring+" -a " + x + " "

    path_file_email=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(deftagger)+"/email.sh"
    file_email=open(path_file_email,"w")
    if not sendplots:
        file_email.write('cat "+an_jonpre+"CAT_SKTreeOutput/' + os.getenv("USER")  + '/CLUSTERLOG' + str(deftagger) + '/email.txt | mail -s "Job summary for job ' + str(deftagger) + " " + e_subject + '" '+str(email_user)+'')
    else:
        file_email.write('cat "+an_jonpre+"CAT_SKTreeOutput/' + os.getenv("USER")  + '/CLUSTERLOG' + str(deftagger) + '/email.txt | mail  '+ plotstring+' -s "Job summary for job ' + str(deftagger) + " " + e_subject + '" '+str(email_user)+'')
    file_email.close()

    filejobsummary = open(an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(deftagger)+"/email.txt","w")
    for eline in jobsummary:
        filejobsummary.write(eline)
    filejobsummary.close()    
    os.system("source " + path_file_email)
    
    
def CrashHelper(defcrashlog):
    
    crashmessage=[]
    for cline in defcrashlog:
        if "vector::_M_range_check" in cline:
            crashmessage=CrashHelpMessage(0)
        if "Did not find"  in cline:
            crashmessage=CrashHelpMessage(1)        
        if "Muon collection not found" in cline:
            crashmessage=CrashHelpMessage(2)
        if "Electron collection not found" in cline:
            crashmessage=CrashHelpMessage(3)
        if "Jet collection not found" in cline:
            crashmessage=CrashHelpMessage(4)
    for mline in crashmessage:
        print mline        
    return crashmessage

def CrashHelpMessage(crashtype):
    crashmessage=[]
    if crashtype == 0:
        crashmessage.append("Crash likedly caused by attempt to access an element of a vector that does not exist")
        crashmessage.append("Check in your code places where you access any vector using vec.at(x) or vec[x] that x < vec.size()")
        crashmessage.append("i.e.,")
        crashmessage.append("std::vector muons = GetMuons()")
        crashmessage.append("if(muons.at(0).Pt() > 20.) ...  --> this would cause crash for events with no muons. You must require muons vector size > 0"  )
        crashmessage.append("__________________________")
        crashmessage.append("if(muons.size() > 0){ ")
        crashmessage.append("   if(muons.at(0).Pt() > 20.) ... This will not crash")
        return crashmessage
    if crashtype == 1: 
        crashmessage.append("Crash caused due to missing file. Check that file is located in $LQANALYZER_FILE_DIR directory (use ls $LQANALYZER_FILE_DIR).")
        crashmessage.append("If file is located in $LQANALYZER_FILE_DIR do:")
        crashmessage.append("cp $LQANALYZER_FILE_DIR/file  $FILEDIR")
        crashmessage.append("If file is not located in $LQANALYZER_FILE_DIR email jalmond@cern.ch with crash message")
    if crashtype == 2:
        crashmessage.append("Crash caused due to Error in GetMuons input string")
    if crashtype == 3:
        crashmessage.append("Crash caused due to Error in GetElectrons input string")
    if crashtype == 4:
        crashmessage.append("Crash caused due to Error in GetJets input string")

    return crashmessage

def GetPartualName(defskim, ismc , defsample, defrunnp, defruncf, defchannel ,defcycle):

    if defskim == "SKTree_NoSkim":
        defskim ="NoCut"
    elif defskim == "SKTree_LeptonSkim":
        defskim ="Lepton"
    elif defskim == "SKTree_DiLepSkim":
        defskim="DiLep"
    elif defskim == "SKTree_TriLepSkim":
        defskim="TriLep"
    
    output_catversion=str(os.getenv("CATVERSION"))
                                              
    if not ismc:
        if defskim == "Lepton":
            defchannel="SK" + defchannel
        else:
            if defskim == "NoCut":
                defchannel="SK" + defchannel + "_nocut"
            else:
                if defskim == "DiLep":
                    defchannel="SK" + defchannel + "_dilep"
                else:
                    if defskim == "TriLep":
                        defchannel="SK" + defchannel + "_trilep"
    else:
        if defskim == "Lepton":
            defsample="SK" + defsample
        else:
            if defskim == "NoCut":
                defsample="SK" + defsample + "_nocut"
            else:
                if defskim == "DiLep":
                    defsample="SK" + defsample + "_dilep"
                else:
                    if defskim == "TriLep":
                        defsample="SK" + defsample + "_trilep"

    outsamplename=""                        
    if ismc:
        outsamplename=defsample
    else:
        outsamplename="period"+defsample

    if defrunnp == "True":
        outsamplename = "nonprompt_" + outsamplename
    if defruncf == "True":
        outsamplename = "chargeflip_" + outsamplename
    if not ismc:
        outsamplename = outsamplename +  "_" + defchannel
        outsamplename = outsamplename + "_cat_" + str(output_catversion)
    else:
        outsamplename = outsamplename + "_cat_"+ str(output_catversion)

    if  "SKTreeMaker" in defcycle:
        outsamplename = outsamplename +  str(os.getenv("CATTAG"))

    return outsamplename
                        

def GetLogFilePath(defskim, ismc , defsample, defrunnp, defruncf, defchannel ,defcycle):
    tmpname= GetPartualName(defskim, ismc , defsample, defrunnp, defruncf, defchannel ,defcycle)
    outlogname = str(os.getenv("LQANALYZER_LOG_PATH")) + "/" + tmpname + "/" + tmpname+"_Job_*"
    return outlogname

def GetOutFileName(defskim, ismc , defsample, defrunnp, defruncf, defchannel ,defcycle, mergedname):
    tmpname= GetPartualName(defskim, ismc , defsample, defrunnp, defruncf, defchannel ,defcycle)
    outsamplename=""
    if ismc:
        outsamplename=  defcycle +"_"+tmpname + ".root"
    elif defrunnp == "True":
        output_catversion=str(os.getenv("CATVERSION"))
        skimtag=""
        if defskim == "DiLep":
            skimtag= "_dilep"
        elif defskim == "TriLep":
            skimtag= "_trilep"
        outsamplename=  defcycle + "_" + defchannel + "_SKnonprompt"+skimtag+ "_cat_" +  str(output_catversion)+ ".root"
        if not mergedname:
            outsamplename=  defcycle + tmpname + ".root"

    elif defruncf == "True":
        output_catversion=str(os.getenv("CATVERSION"))
        skimtag=""
        if defskim == "DiLep":
            skimtag= "_dilep"
        elif defskim == "TriLep":
            skimtag= "_trilep"
        outsamplename=  defcycle + "_" + defchannel + "SKchargeflip_"+skimtag+ "_cat_" +  str(output_catversion)+ ".root"
        if not mergedname:
            outsamplename=  defcycle + tmpname + ".root"
                        

    else:
        output_catversion=str(os.getenv("CATVERSION"))
        outsamplename=  defcycle +"_data_"+ defchannel + "_cat_" +  str(output_catversion)+ ".root"
        if not  mergedname:
            outsamplename=defcycle +"_"+tmpname + ".root"
     
    if "SKTreeMaker" in defcycle:
        if ismc:
            if defskim == "DiLep":
                return an_jonpre+"/CatNtuples/"+str(os.getenv("CATVERSION"))+"/SKTrees/MCDiLep/"+defsample+"/"
            elif defskim == "TriLep":
                return an_jonpre+"/CatNtuples/"+str(os.getenv("CATVERSION"))+"/SKTrees/MCTriLep/"+defsample+"/"
            else:
                return an_jonpre+"/CatNtuples/"+str(os.getenv("CATVERSION"))+"/SKTrees/MC/"+defsample+"/"
                                                    
    return outsamplename    

def LargeFileSize(fsize):
    sfilesize=str(fsize)
    string_length= len(sfilesize)
    float_only_fsize= sfilesize[:-1]
    unit_only_fsize=  sfilesize[(int(string_length) - 1):]
    if unit_only_fsize == "G":
        return True
    if  unit_only_fsize == "M":
        if float(float_only_fsize)> large_file_size:
            return True
    return False    

def FileSizeInB(fsize):
    unit_remove=1
    if "B" in str(fsize):
        unit_remove=2
    sfilesize=str(fsize)
    string_length= len(sfilesize)
    float_only_fsize= sfilesize[:-unit_remove]
    unit_only_fsize=  sfilesize[(int(string_length) - unit_remove):]
    
    if "K" in  unit_only_fsize:
        return float(float_only_fsize)*1000.
    elif "M" in  unit_only_fsize:
        return float(float_only_fsize)*1000000.
    elif "G" in  unit_only_fsize:
        return float(float_only_fsize)*1000000000.
    else:
        return float(float_only_fsize)

def LargeFileSizeIncrease(fsize,defsample, defcycle, tagger,defskim):

    if GetFileSize(defsample, defcycle,tagger,defskim) == "None":
        return False
    file_size_prev=   FileSizeInB(GetFileSize(defsample, defcycle, tagger,defskim))
    file_size = FileSizeInB(fsize)
    if file_size > file_increase_warning*file_size_prev:
        return True
    else:
        return False
    

def TimeIncrease(stime,defsample, defcycle, tagger,defskim):
    time_prev=GetTime(defsample, defcycle, tagger,defskim)
    if time_prev == "None":
        return False

    if float(stime) > float(time_prev)*time_increase_warning:
        if float(stime) < 30.:
            ##### if time is small this can be an issue with machine not code
            return False
        else:
            return True
    else:
        return False
    
    
def LargeMemory(mem):
    if "GB" in str(mem):
        return True
    
    smem=str(mem)
    string_length= len(smem)
    float_only_mem= smem[:-2]
    if "MB" in smem:
        if float(float_only_mem) > large_memory_check:
            return True
    else:
        return False


def LargePhysicalMemoryIncrease(mem,defsample, defcycle, tagger,defskim):

    if GetPhysicalMemoryUsage(defsample, defcycle, tagger,defskim) == "None":
        return False
    
    if not "MB" in str(mem):
        return False
    smem=str(mem)
    string_length= len(smem)
    float_only_mem= smem[:-2]
    unit_only_mem=  smem[(int(string_length)-2):]

    rounded_float= str(round(float(float_only_mem),2))
    if unit_only_mem == "MB":
        if float(float_only_mem)  < 300:
            if (float(float_only_mem)*1000.) > FileSizeInB(GetPhysicalMemoryUsage(defsample, defcycle, tagger,defskim))*1.4:
                return True
            else:
                return False
        else:
             if (float(float_only_mem)*1000.)  > FileSizeInB(GetPhysicalMemoryUsage(defsample, defcycle, tagger,defskim))*1.3:
                 return True
             else:
                 return False
    return False 
    

def LargeVirtualMemoryIncrease(mem,defsample, defcycle, tagger,defskim):

    if GetVirtualMemoryUsage(defsample, defcycle,tagger,defskim) == "None":
        return False
    if not "MB" in str(mem):
        return False
    smem=str(mem)
    string_length= len(smem)
    float_only_mem= smem[:-2]
    unit_only_mem=  smem[(int(string_length)-2):]

    rounded_float= str(round(float(float_only_mem),2))
    if unit_only_mem == "MB":
        if float(float_only_mem) < 300:
            if (float(float_only_mem)*1000.) > (FileSizeInB(GetVirtualMemoryUsage(defsample, defcycle, tagger,defskim))*1.4):
                return True
            else:
                return False
        else:
             if (float(float_only_mem)*1000.) > ((FileSizeInB(GetVirtualMemoryUsage(defsample, defcycle,tagger,defskim))*1.3)):
                 return True
             else:
                 return False
    return False



def CheckRunningStatus(check_crashfile):
    if not os.path.exists(check_crashfile):
        return
    file_check_crashfile=open(check_crashfile,"r")
    for line in file_check_crashfile:
        if " No such file or directory:" in line:
            print line 
            sys.exit()
        if "is not defined" in line:
            print line 
            sys.exit()
        if "SyntaxError: invalid syntax" in line:
            print line 
            sys.exit()
    file_check_crashfile.close()            
    return

def RoundMemory(mem):
    if not "MB" in str(mem):
        return str(mem)
    smem=str(mem)
    string_length= len(smem)
    float_only_mem= smem[:-2]
    unit_only_mem=  smem[(int(string_length)-2):]
    
    rounded_float= str(round(float(float_only_mem),2)) + unit_only_mem
    return rounded_float
    

#Import parser to get options                                                                                                                                                 
parser = OptionParser()
parser.add_option("-p", "--period", dest="period", default="A",help="which data period or mc sample")
parser.add_option("-s", "--stream", dest="stream", default="NULL", help="Which data channel- ee,or mumu?")
parser.add_option("-j", "--jobs", dest="jobs", default=1, help="Name of Job")
parser.add_option("-c", "--cycle", dest="cycle", default="Analyzer", help="which cycle")
parser.add_option("-t", "--tree", dest="tree", default="ntuple/event", help="What is input tree name?")
parser.add_option("-o", "--logstep", dest="logstep", default=-1, help="How many events betwene log messages")
parser.add_option("-d", "--data_lumi", dest="data_lumi", default="A", help="How much data are you running on/ needed to weight mc?")
parser.add_option("-l", "--loglevel", dest="loglevel", default="INFO", help="Set Log output level")
parser.add_option("-n", "--nevents", dest="nevents", default=-1, help="Set number of events to process")
parser.add_option("-k", "--skipevent", dest="skipevent", default=-1, help="Set number of events to skip")
parser.add_option("-a", "--datatype", dest="datatype", default="", help="Is data or mc?")
parser.add_option("-e", "--totalev", dest="totalev", default=-1, help="How many events in sample?")
parser.add_option("-x", "--xsec", dest="xsec", default=-1., help="How many events in sample?")
parser.add_option("-X", "--tagger", dest="tagger", default="123", help="random number string?")
parser.add_option("-T", "--targetlumi", dest="targetlumi", default=-1., help="How many events in sample?")
parser.add_option("-E", "--efflumi", dest="efflumi", default=-1., help="How many events in sample?")
parser.add_option("-O", "--outputdir", dest="outputdir", default="${LQANALYZER_DIR}/data/output/", help="Where do you like output to go?")
parser.add_option("-w", "--remove", dest="remove", default=True, help="Remove the work space?")
parser.add_option("-S", "--skinput", dest="skinput", default=True, help="Use SKTree as input?")
parser.add_option("-R", "--runevent", dest="runevent", default=True, help="Run Specific Event?")
parser.add_option("-N", "--useCATv742ntuples", dest="useCATv742ntuples", default=True, help="' to run on these samples")
parser.add_option("-L", "--LibList", dest="LibList", default="NULL", help="Add extra lib files to load")
parser.add_option("-D", "--debug", dest="debug", default=False, help="Run submit script in debug mode?")
parser.add_option("-m", "--useskim", dest="useskim", default="Lepton", help="Run submit script in debug mode?")
parser.add_option("-P", "--runnp", dest="runnp", default="runnp", help="Run fake mode for np bkg?")
parser.add_option("-Q", "--runcf", dest="runcf", default="runcf", help="Run fake mode for np bkg?")
parser.add_option("-q", "--queue", dest="queue", default="", help="Which queue to use?")
parser.add_option("-v", "--catversion", dest="catversion", default="NULL", help="What cat version?")
parser.add_option("-f", "--skflag", dest="skflag", default="NULL", help="add input flag?")
parser.add_option("-b", "--usebatch", dest="usebatch", default="usebatch", help="Run in batch queue?")
parser.add_option("-u", "--useremail", dest="useremail", default="", help="Set user email")
parser.add_option("-B", "--bkg", dest="bkg", default="False", help="run in bkg")
parser.add_option("-A","--drawhists",dest="drawhists",default="False", help="draw nothing")

#curses.resizeterm

job_time=0.

now = datetime.datetime.now()
diff = datetime.timedelta(days=3)
diffweek = datetime.timedelta(days=7)
future = now + diff
future_week = now + diffweek
future=future.strftime("%m/%d/%Y")
future_week=future_week.strftime("%m/%d/%Y")


###################################################                                                                                                                            
#set the local variables using options                                                                                                                                         
###################################################                                                                                                                            
(options, args) = parser.parse_args()
number_of_cores = int(options.jobs)
sample = options.period
channel = options.stream
cycle = options.cycle
logstep = int(options.logstep)
loglevel = options.loglevel
runnp = options.runnp
runcf = options.runcf
queue = options.queue
tagger= options.tagger
useremail=options.useremail
### THESE ARE OPTIONS THAT CAN BE INCLUDED but not in example                                                                                                                  
tree = options.tree
number_of_events_per_job= int(options.nevents)
skipev = int(options.skipevent)
dataType = options.datatype
totalev = int(options.totalev)
xsec = float(options.xsec)
tar_lumi = float(options.targetlumi)
eff_lumi = float(options.efflumi)
data_lumi = options.data_lumi
catversion = options.catversion
Finaloutputdir = options.outputdir
remove_workspace=options.remove
useskinput=options.skinput
runevent= options.runevent
useCATv742ntuples = options.useCATv742ntuples
tmplist_of_extra_lib=options.LibList
DEBUG = options.debug
useskim = options.useskim
skflag = options.skflag
usebatch =options.usebatch
runinbkg=options.bkg
quick_draw=options.drawhists


queuepath=path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CattupleConfig/QUEUE/ForceQueue.txt"
file_queuepath = open(queuepath,"r")
for line in file_queuepath:
    if "#" in line:
        continue
    sline = line.split()
    for s in sline:
        if "@" in s:
            queue= s
file_queuepath.close()


printedqueue=queue
printedqueue=printedqueue.replace("all.q@cms-0-","node")
printedqueue=printedqueue.replace(".local","")
printedqueue=printedqueue.replace("node","")
if printedqueue == "None":
    printedqueue = "fastq"

DoSendEmail=False
run_in_bkg=False
if runinbkg == "True":
    run_in_bkg=True
    DoSendEmail=True
    print "Running in background"
    if not  os.getenv("USER") == "jalmond":
        if "jalmond" in useremail:
            print "Email could not be sent since email address is not set correctly in bin/catconfig."
            print "Please set email address and resubmit job"
            sys.exit()

skim_print=useskim
if useskim == "SKTree_LeptonSkim":
    skim_print ="Lepton   "
elif useskim == "SKTree_DiLepSkim":
    skim_print="DiLep. "
elif useskim == "SKTree_TriLepSkim":
    skim_print="TriLep. "


sample = sample.replace(":", " ")
sample = sample.replace("!!", " ")
sample = sample.split()

for s in sample:
    if not s:
        sys.exit()

istatus_message=2
winx = 5*len(sample) + 22 +  istatus_message
winy = 180

crash_output=[]
crash_outputjob=[]

output_bkg=[]
for out_x in range(1,winx):
    output_bkg.append(" "*170)

start_time = time.time()  
screen = curses.initscr()


screen.refresh()
stdscr = curses.newwin(winx, winy)
if not run_in_bkg:
    curses.noecho()
    curses.cbreak()

box_shift=0
summary_block0=4 + box_shift
summary_block1=summary_block0+ len(options.cycle) + 5
summary_block2=summary_block1+20
summary_block3=summary_block2+10
summary_block4=summary_block3+15
summary_block5=summary_block4+27
summary_block6=summary_block5+8
summary_block7=summary_block6+60

stdscr.addstr(0, box_shift,  "Job[" + str(tagger)+"] Summary:",curses.A_UNDERLINE)
stdscr.addstr(1, box_shift, "Job ",curses.A_STANDOUT)
stdscr.addstr(1, summary_block0, "| Code Name                   ",curses.A_STANDOUT)
stdscr.addstr(1, summary_block1, "| Sample Name                   ",curses.A_STANDOUT)
stdscr.addstr(1, summary_block2, "| Skim    ",curses.A_STANDOUT)
stdscr.addstr(1, summary_block3, "| Status                  ",curses.A_STANDOUT)
stdscr.addstr(1, summary_block4, "| Cluster ID Range | njobs  ",curses.A_STANDOUT)
stdscr.addstr(1, summary_block5, "| Queue  ",curses.A_STANDOUT)
stdscr.addstr(1, summary_block6, "| Cluster Progess                                            ",curses.A_STANDOUT)
stdscr.addstr(1, summary_block7, "|" ,curses.A_STANDOUT)
stdscr.refresh()


#### list 2 is more technical properties                                                                                                                                                                                                                                        
summary2_block0=4+box_shift
summary2_block1=summary2_block0+20
summary2_block2=summary2_block1+20
summary2_block3=summary2_block2+20
summary2_block4=summary2_block3+20
summary2_block5=summary2_block4+20


list2= istatus_message + 2 + len(sample)
stdscr.addstr(list2, box_shift,  "PostJob Summary:",curses.A_UNDERLINE)
stdscr.addstr(1+list2, box_shift, "Job ",curses.A_STANDOUT)
stdscr.addstr(1+list2, summary2_block0, "| Cum.Process. Time  ",curses.A_STANDOUT)
stdscr.addstr(1+list2, summary2_block1, "| Job Process. Time  ",curses.A_STANDOUT)
stdscr.addstr(1+list2, summary2_block2, "| Physical Memory             ",curses.A_STANDOUT)
stdscr.addstr(1+list2, summary2_block3, "| Virtual Memory         ",curses.A_STANDOUT)
stdscr.addstr(1+list2, summary2_block4, "| OutputFile Size    ",curses.A_STANDOUT)
stdscr.addstr(1+list2, summary2_block5, "|",curses.A_STANDOUT)

list3=2*len(sample)+11 + istatus_message 
list3b=3*len(sample)+16 + istatus_message 
list3c=4*len(sample)+18 + istatus_message 
list4=5*len(sample)+20 + istatus_message 


remove_from_end= 3 

stdscr.addstr(list4,box_shift,  "Job Status (Summary of latest background process):" ,curses.A_UNDERLINE)
#### make a random number for the job ID. This will help make a directory to work in

stdscr.addstr(list3, box_shift,  "Job Terminal Output: (will be deleted "+ future +")" ,curses.A_UNDERLINE)
stdscr.addstr(list3b, box_shift,  "Log Files: (will be deleted "+ future_week +")" ,curses.A_UNDERLINE)
stdscr.addstr(list3c, box_shift,  "Output Files:" ,curses.A_UNDERLINE)

#### make working directorr
if not os.path.exists(an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)):
    os.system("mkdir " + an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  +"/CLUSTERLOG" +str(tagger))
else:
    os.system("rm -r " + an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  +"/CLUSTERLOG" +str(tagger))
    os.system("mkdir " + an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  +"/CLUSTERLOG" +str(tagger))

### Make directory for job stats
statdir=path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/"+ os.getenv("USER")
if not os.path.exists(statdir):
    os.system("mkdir " +statdir)

if channel == "":
    channel = NULL

if not os.path.exists(path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/MasterFile_"+ os.getenv("CATVERSION")+".txt"):
    os.system("cp  " + path_jobpre +"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/MasterFileSkeleton.txt " + path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/MasterFile_"+ os.getenv("CATVERSION")+".txt")


jobinfo_file_master = path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/MasterFile_"+ os.getenv("CATVERSION")+".txt"
jobinfo_file= path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + os.getenv("USER") + "/PreMasterFile_"+os.getenv("CATVERSION")+ str(tagger)+".txt"
os.system("cp " + jobinfo_file_master + " " + jobinfo_file)



#http://stackoverflow.com/questions/14300770/how-to-noutrefresh-the-multi-line-output-dynamically
### split sample and check cluster queue

#an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/" + str(tagger)+ "/jobid.txt"
samples_inbackground=[]
samples_complete=[]
jobidcrash=[]
isample=-1

useskimtmp=""
ismctmp=  True
sampletmp=""
runnptmp=""
runcftmp=""
channeltmp=""

output_warning=[]

#### Loop over all samples in job
quickdraw=False
if quick_draw == "True":
    quickdraw=True

fileoutputlist=[]

for s in sample:
    if not s:
        continue
    ncopies=0
    for s2 in sample:
        if s == s2:
            ncopies=ncopies+1
    if ncopies > 1:
        if not run_in_bkg:
            curses.echo()
            curses.nocbreak()
            curses.endwin()
        print "%"*45
        print "List of samples contains duplicates. Please fix list"
        print "%"*45
        sys.exit()

runningData=False
longestjob=0.
for s in sample:
    stime=GetAverageTime(tagger, s, cycle,useskim)
    if stime > longestjob:
        longestjob=stime

for s in sample:
    #### Get number of subjobs from DetermineNjobs function. Unless number_of_cores is set to 1 this will check the processing time of the cycle and batch queue to determin the number of jobs to run
    njobs_for_submittion=DetermineNjobs(longestjob,number_of_cores, tagger, s, cycle,useskim)
    isample=isample+1

    ## set MC bool from the sample length. This is the letter of the data period for data
    isMC = len(s) > 1
    if s == "H_v2" or s == "H_v3":
        isMC= False

    runningData= not isMC
    
    if run_in_bkg:
        del output_bkg[:]
        for out_x in range(1,winx):
            output_bkg.append(stdscr.instr(out_x, 0))
        UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")    
    #### if sample is submitted this for loop will check job process on cluster and update the screen/terminal
    for x in range(0, len(sample)):
        
        if run_in_bkg:
            del output_bkg[:]
            for out_x in range(1,winx):
                output_bkg.append(stdscr.instr(out_x, 0))
            UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")

        ### check if sample is already complet. Only update terminal for samples not complete
        sample_complete=False
        for finsamp in samples_complete:
            if sample[x] == finsamp:
                sample_complete=True
        if sample_complete:
            continue
        
        #### Check that the sample is running. Do not update terminal for samples that are not submitted to batch queue
        backgroundsamples=False
        for bs in samples_inbackground:
            if bs == sample[x]:
                backgroundsamples=True
        ### the following variables are used to print the output sample name that is determined in the localsubmit.py

        useskimtmp=useskim
        ismctmp=  len(sample[x])>1
        if sample[x] == "H_v2" or sample[x] == "H_v3":
            ismctmp= False

        sampletmp=sample[x]
        runnptmp=runnp
        runcftmp=runcf
        channeltmp=channel

        ### if sample has been submitted to cluster backgroundsamples will be True        
        if backgroundsamples:
            #### in case the job id file is not yet filled add while condition
            jobid_exists=True
            while jobid_exists:        
                path_job_check=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] + "jobid.txt"  
                if  os.path.exists(path_job_check):
                    jobid_exists=False
                    


            #CheckRunningStatus(an_jonpre+"/CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) +"/" + tagger + "/" + sample[x] +".txt")

            file_job_check=open(path_job_check ,"r")
            jobid1=0
            jobid2=0
            #### check status of job on cluster
            njobs_finished=0.
            nrunning=0.
            nqueue=0.
            njobs_in_total=0.  ### should equal three above
            path_clust_check=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] + "clust.txt"
            os.system("qstat -u " + os.getenv("USER") + " > " +  path_clust_check)
            ijob=0

            ### loop over job IDS
            for sline in file_job_check:
                ssline = sline.split()
                if len(ssline) < 1:
                    continue
                line=ssline[0]
                if jobid1 == 0:
                    jobid1=line
                jobid2=line    
                ### line == job ID
                ijob=ijob+1.
                njobs_in_total=njobs_in_total+1.
                file_clust_check=open(path_clust_check,"r")
                job_inqueue=False
                for cline in file_clust_check:
                    splitline  = cline.split()
                    if len(splitline) < 6:
                        continue
                    if not os.getenv("USER")  in cline:
                        continue
                    if line == splitline[0]:
                        job_inqueue=True
                        if splitline[4] == "r":
                            nrunning = nrunning + 1.
                        else:
                            nqueue = nqueue + 1.
                if not job_inqueue:
                    njobs_finished=njobs_finished+1.
                file_clust_check.close()
            path_clust_check=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] + "clust.txt"
            os.system("rm " + path_clust_check)
    
            if njobs_in_total == 0:
                continue
            file_job_check.close()
            if njobs_finished == njobs_in_total:
                samples_complete.append(sample[x])
                
            ### add finished jobs to running jobs status since they have ran    
            nrunning = nrunning + njobs_finished    
            nrun_per = round( (nrunning/float(njobs_in_total)), 2)
            nfin_per = round( (njobs_finished/float(njobs_in_total)), 2)

            nscreen_run= "["
            nscreen_fin= "["
            for xsr in range(0, 9):
                if nrun_per > float(xsr)* 0.1:
                    nscreen_run+="#"
                else:
                    nscreen_run+="-"
            if nrun_per == 1.:
                nscreen_run+="#"
            else:
                nscreen_run+="-"
            nscreen_run+="]"
            for xsr in range(0, 9):
                if nfin_per > float(xsr)* 0.1:
                    nscreen_fin+="#"
                else:
                    nscreen_fin+="-"
            if nfin_per == 1.:
                nscreen_fin+="#"
            else:
                nscreen_fin+="-"
            nscreen_fin+="]"


            if nrun_per == 0.:
                stdscr.addstr(int(x)+istatus_message, box_shift,  str(int(x+1)), curses.A_NORMAL )
                stdscr.addstr(int(x)+istatus_message, summary_block0, "| " +  cycle, curses.A_NORMAL )
                if isMC:
                    stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + sample[x], curses.A_NORMAL )
                else:
                    stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + channel+"_"+sample[x], curses.A_NORMAL )
                stdscr.addstr(int(x)+istatus_message, summary_block2, "| " + skim_print, curses.A_NORMAL )
                stdscr.addstr(int(x)+istatus_message, summary_block3,"| Submitted " , curses.A_NORMAL) 
                stdscr.refresh()
            else:
                stdscr.addstr(int(x)+istatus_message, box_shift,  str(int(x+1)), curses.A_BOLD )
                stdscr.addstr(int(x)+istatus_message, summary_block0,  "| " + cycle, curses.A_BOLD )
                if isMC:
                    stdscr.addstr(int(x)+istatus_message, summary_block1,  "| " + sample[x], curses.A_BOLD )
                else:
                    stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + channel+"_"+sample[x], curses.A_NORMAL )

                stdscr.addstr(int(x)+istatus_message, summary_block2,  "| " + skim_print, curses.A_BOLD )
                stdscr.addstr(int(x)+istatus_message, summary_block3,"| Running    " , curses.A_BLINK) 
                if not jobid1 == 0:
                    nblanks=2*(8-len(str(jobid2)))
                    stdscr.addstr(int(x)+istatus_message, summary_block4,"| " + str(jobid1) + "-" + str(jobid2) + " "*nblanks +   "| " + str(1+int(jobid2)-int(jobid1)) , curses.A_BOLD)
                    stdscr.addstr(int(x)+istatus_message, summary_block5,"| " + printedqueue  , curses.A_BOLD)
                else:
                    stdscr.addstr(int(x)+istatus_message, summary_block4,"|" + " "*18 + "|"  , curses.A_BOLD)
                    stdscr.addstr(int(x)+istatus_message, summary_block5,"| " + printedqueue  , curses.A_BOLD)
                stdscr.addstr(int(x)+istatus_message, summary_block6,  "| Running   " + nscreen_run + " " + str(100*nrun_per)+ "%  Complete" + nscreen_fin + " " + str(100*nfin_per) + "%",curses.A_BOLD)
                stdscr.addstr(int(x)+istatus_message, summary_block7 ,"|    ",curses.A_BOLD)    
                stdscr.refresh()
                if run_in_bkg:
                    del output_bkg[:]
                    for out_x in range(1,winx):
                        output_bkg.append(stdscr.instr(out_x, 0))
                        UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")
            if nfin_per == 1.:
                stdscr.addstr(int(x)+istatus_message, box_shift,  str(int(x+1)), curses.A_DIM )
                stdscr.addstr(int(x)+istatus_message, summary_block0,  "| " + cycle, curses.A_DIM )
                if isMC:
                    stdscr.addstr(int(x)+istatus_message, summary_block1,  "| " + sample[x], curses.A_DIM )
                else:
                    stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + channel+"_"+sample[x], curses.A_NORMAL )

                stdscr.addstr(int(x)+istatus_message, summary_block2,  "| " + skim_print, curses.A_DIM )
                stdscr.addstr(int(x)+istatus_message, summary_block3,"| Complete " , curses.A_DIM) 
                job_time = time.time() - start_time
                if not jobid1 == 0:    
                    nblanks=2*(8-len(str(jobid2)))
                    stdscr.addstr(int(x)+istatus_message, summary_block4,"| " + str(jobid1) + "-" + str(jobid2)  + " "*nblanks +   "| " + str(1+int(jobid2)-int(jobid1)), curses.A_DIM)
                    stdscr.addstr(int(x)+istatus_message, summary_block5,"| " + printedqueue , curses.A_BOLD)
                else:
                    stdscr.addstr(int(x)+istatus_message, summary_block6,"| "  , curses.A_DIM)
 
                stdscr.addstr(int(x)+istatus_message, summary_block6,  "| Running   " + nscreen_run + " " + str(100*nrun_per)+ "%  Complete" + nscreen_fin + " " + str(100*nfin_per) + "%",curses.A_DIM) 
                stdscr.refresh()

                path_job=path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser() + "/" + str(tagger)+ "/statlog_time_"+sample[x] + tagger + ".txt"
                ismerging=True
                while not os.path.exists(path_job):
                    if ismerging:
                        stdscr.addstr(2+list2+int(x), box_shift , str(int(x+1)),curses.A_DIM)
                        stdscr.addstr(2+list2+int(x), summary2_block1 ,"| MERGING OUTPUT ", curses.A_BLINK)
                        stdscr.refresh()
                        ismerging=False
                    check_crash_stat=False
                    if check_crash_stat:
                        check_crashfile=an_jonpre+"/CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) +"/" + tagger + "/" + sample[x] +".txt"
                        if not os.path.exists(check_crashfile):
                            continue
                        file_check_crashfile=open(check_crashfile,"r")
                        crash_in_job=False
                        crashlog=[]
                        for line in file_check_crashfile:
                            if " No such file or directory:" in line:
                                crash_in_job=True
                                crashlog.append(line)
                        file_check_crashfile.close()
                        if crash_in_job:
                            crash_log= an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] +"_crash/crashlog.txt"
                            if os.path.exists(crash_log):
                                file_crash=open(crash_log,"a")
                                for linex in  crashlog:
                                    file_crash.write(linex)
                                file_crash.close()
                            else:
                                file_crash=open(crash_log,"w")
                                for linex in  crashlog:
                                    file_crash.write(linex)
                                file_crash.close()
                            break

                    continue
                

                crash_log= an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] +"_crash/crashlog.txt"

                if os.path.exists(crash_log):
                    jobidcrash.append(int(x))
                    stdscr.addstr(2+list2+int(x), box_shift , str(int(x+1)),curses.A_DIM)
                    stdscr.addstr(2+list2+int(x), summary2_block1 ,"| JOB CRASHED    ", curses.A_BLINK)
                    logpath=GetLogFilePath(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle)
                    outfilepath=str(Finaloutputdir)  +GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,False)
                    fileoutputlist.append(outfilepath)
                    if "SKTreeMaker" in cycle:
                        outfilepath=GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,False)
                    stdscr.addstr(list3b + 1+int(x), box_shift,  str(1+int(x)) +": Log files for " + sample[x] + " found at " + str(logpath))
                    stdscr.addstr(list3c + 1+int(x), box_shift,   str(1+int(x)) +":###CRASH###: No OutputFile for " + sample[x])
                    crash_output.append( logpath)
                    crash_outputjob.append(str(jobid1))
                    stdscr.refresh()
                else:
                    stdscr.addstr(2+list2+int(x), box_shift , str(int(x+1)),curses.A_DIM) 
                    stdscr.addstr(2+list2+int(x), summary2_block0 ,"| " + str(round(job_time,2)) + "[s]",curses.A_DIM) 
                    stdscr.refresh()
                    logpath=GetLogFilePath(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle)
                    outfilepath=str(Finaloutputdir) + GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,False)
                    fileoutputlist.append(outfilepath)
                    if "SKTreeMaker" in cycle:
                        outfilepath=GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,False)

                    stdscr.addstr(list3b + 1+int(x), box_shift,  str(1+int(x)) +": Log files for " + sample[x] + " found at " + str(logpath))
                    stdscr.addstr(list3c + 1+int(x), box_shift,   str(1+int(x)) +": OutputFile for " + sample[x] + " = " + str(outfilepath))
                    stdscr.refresh()
                    
                    filesize=0.
                    memoryusage_p=0.
                    memoryusage_v=0.
                    outputfile_size=0.
                    first_job_time=0.
                    last_job_time=0.
                    file_job=open(path_job,"r")
                    for line in file_job:
                        splitline  = line.split()
                        if len(splitline) < 2:
                            continue
                        if "job_time"  == splitline[0]:
                            first_job_time=float(splitline[1])
                        if "last_job_time"  == splitline[0]:
                            last_job_time=float(splitline[1])
                        if "outputfile_size" == splitline[0]:
                            outputfile_size=splitline[1]
                        if "memoryusage_p" ==  splitline[0]:
                            memoryusage_p=splitline[1]    
                        if "memoryusage_v" ==  splitline[0]:
                            memoryusage_v=splitline[1]
                    file_job.close()
                    stdscr.addstr(2+list2+int(x), summary2_block1 ,"| " + str(round(first_job_time,2)) + "-" +  str(round(last_job_time,2)) + "[s]   ",curses.A_DIM) 

                    if LargeMemory(memoryusage_p):
                        output_warning.append("WARNING:: Physical memory of job " + sample[x]+ "  is excess of " + str(large_memory_check) + " MB.")
                    if LargeMemory(memoryusage_v):
                        output_warning.append("WARNING:: Virtual memory of job " + sample[x] + " is excess of " + str(large_memory_check) + " MB.")

                    if isMC:
                        if LargePhysicalMemoryIncrease(memoryusage_p,sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Physical memory of job " + sample[x] + " is " + RoundMemory(memoryusage_p) + " compared to previous jobs " + GetPhysicalMemoryUsage(sample[x], cycle,tagger,useskim) + ". This is likely a memory leak introduced recently")
                        if LargeVirtualMemoryIncrease(memoryusage_v,sample[x], cycle,tagger,useskim):
                                output_warning.append("WARNING:: Virtual memory of job " + sample[x] + " is " + RoundMemory(memoryusage_v) + " compared to previous jobs " + GetVirtualMemoryUsage(sample[x], cycle,tagger,useskim) + ". This is likely a memory leak introduced recently")    
                        if LargeFileSize(outputfile_size):
                            output_warning.append("WARNING:: Size of output rootfile for " + sample[x] + " is " + str(outputfile_size) + " greater than " + str(large_file_size) + " MB")
                        if LargeFileSizeIncrease(outputfile_size,sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Size of output rootfile for " + sample[x] + " is " + str(outputfile_size) + " compared to previous jobs " + GetFileSize(sample[x], cycle,tagger,useskim))
                        if TimeIncrease(first_job_time,sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Job time per input file increased for " + sample[x] + "  from " + str(first_job_time) + " compared to previous job time " + GetTime(sample[x], cycle,tagger,useskim))
                    else:
                        if LargePhysicalMemoryIncrease(memoryusage_p,channel + "_" + sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Physical memory of job " + sample[x] + " is " + RoundMemory(memoryusage_p) + " compared to previous jobs " + GetPhysicalMemoryUsage(channel + "_" + sample[x], cycle,tagger,useskim) + ". This is likely a memory leak introduced recently")
                        if LargeVirtualMemoryIncrease(memoryusage_v,channel + "_" + sample[x], cycle,tagger,useskim):
                                output_warning.append("WARNING:: Virtual memory of job " +sample[x] +" is " + RoundMemory(memoryusage_v) + " compared to previous jobs " + GetVirtualMemoryUsage(channel + "_" + sample[x], cycle,tagger,useskim) + ". This is likely a memory leak introduced recently")
                        if LargeFileSize(outputfile_size):
                            output_warning.append("WARNING:: Size of output rootfile for " + sample[x] + " is " + str(outputfile_size) + " greater than " + str(large_file_size) + " MB")
                        if LargeFileSizeIncrease(outputfile_size,channel + "_" + sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Size of output rootfile for " + sample[x] + " is " + str(outputfile_size) + " compared to previous jobs " + GetFileSize(channel + "_" + sample[x], cycle,tagger,useskim))
                        if TimeIncrease(first_job_time,channel + "_" + sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Job time per input file for  " + sample[x] + " increased from " + str(first_job_time) + " compared to previous job time " + GetTime(vsample[x], cycle,tagger,useskim))
     
                        
                    stdscr.addstr(2+list2+int(x), summary2_block2 ,"| " + RoundMemory(memoryusage_p),curses.A_DIM) 
                    stdscr.addstr(2+list2+int(x), summary2_block3 ,"| " + RoundMemory(memoryusage_v),curses.A_DIM) 
                    stdscr.addstr(2+list2+int(x), summary2_block4 ,"| " + str(outputfile_size),curses.A_DIM) 
                    stdscr.addstr(2+list2+int(x), summary2_block5 ,"| ",curses.A_DIM)
                    stdscr.refresh()

                    #logpath=GetLogFilePath(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle)
                    #outfilepath=str(Finaloutputdir)  +GetOutFilName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle)
                    #stdscr.addstr(list3b + 1+int(x), box_shift,  str(1+int(x)) +": Log files for " + sample[x] + " found at " + str(logpath))
                    #stdscr.addstr(list3c + 1+int(x), box_shift,   str(1+int(x)) +": OutputFile for " + sample[x] + " = " + str(outfilepath))

        else:
            stdscr.addstr(int(x)+istatus_message, box_shift,  str(int(x+1)) )
            stdscr.addstr(int(x)+istatus_message, summary_block0,  "| " + cycle )
            if isMC:
                stdscr.addstr(int(x)+istatus_message, summary_block1,  "| " + sample[x] )
            else:
                stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + channel+"_"+sample[x], curses.A_NORMAL )
            stdscr.addstr(int(x)+istatus_message, summary_block2,  "| " + skim_print )
            stdscr.addstr(int(x)+istatus_message, summary_block3,"| In Queue " )
            stdscr.addstr(int(x)+istatus_message, summary_block4,"|" + " "*18 + "|" )
            stdscr.addstr(int(x)+istatus_message, summary_block5,"|    ")
            stdscr.addstr(int(x)+istatus_message, summary_block6,"|    ")
            stdscr.addstr(int(x)+istatus_message, summary_block7,"|    ")
            stdscr.refresh()

    ##### New Code looks at submittnig jobs to batch queue        
    logfile=path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + os.getenv("USER") + "/" + str(tagger)+ "/statlog_"+ s + tagger + ".txt"
    if not os.path.exists(path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + os.getenv("USER") + "/" + str(tagger)):
        os.system("mkdir " + path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + os.getenv("USER") + "/" + str(tagger))

    filestatlog=open(logfile,"w")
    filestatlog.write("user " +  os.getenv("USER")+ " \n")
    filestatlog.write(" -c " + cycle + " \n")
    filestatlog.write(" -v " + catversion + " \n")
    filestatlog.write(" -s " + channel + " \n")
    filestatlog.write(" -d " + data_lumi + " \n")
    filestatlog.write("sample "+ s + " \n")
    filestatlog.write(" -m " + useskim + " \n")
    filestatlog.write("cattag " + os.getenv("CATTAG") + " \n")
    filestatlog.write(time.strftime("%c")  + " \n")
    filestatlog.write("############################" + " \n")
    filestatlog.close()
    blankbuffer = "         "
    if not queue:
        queue="None"
    command1= "python  " +  os.getenv("LQANALYZER_DIR")+  "/python/CATConfig.py -p " + s + "  -s " + str(channel) + "  -j " + str(njobs_for_submittion) + " -c  " + str(cycle)+ " -o " + str(logstep)+ "  -d " + str(data_lumi) + " -O " + str(Finaloutputdir) + "  -w " + str(remove_workspace)+ " -l  " + str(loglevel) + "  -k " + str(skipev) + "  -n " + str(number_of_events_per_job) + "  -e " + str(totalev) + "  -x " + str(xsec) + "  -T " + str(tar_lumi) + " -E " + str(eff_lumi) + "  -S " + str(useskinput) + " -R " + str(runevent)+ "  -N " + str(useCATv742ntuples) + " -L " + str(tmplist_of_extra_lib) + " -D " + str(DEBUG) + " -m " + str(useskim) + " -P  " + str(runnp) + " -Q " + str(runcf) + " -v " + str(catversion) + " -f " + str(skflag) + " -b " + str(usebatch) + "  -X " + str(tagger) +" -q " + str(queue)
    command2=command1
    command2 = command2.replace("CATConfig.py", "localsubmit.py")
    command2_background=command2 + "&>  "+an_jonpre+"CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) +"/" + tagger + "/" + s+".txt&"
    checkqueue=True
    stdscr.addstr(list4+1, box_shift,  "Initialise:: sample " + s +  blankbuffer)
    stdscr.refresh()
    while checkqueue:
        os.system(command1)
        if not os.path.exists(an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER") +"/"+tagger):
            stdscr.addstr(list4+1, box_shift,  "Queue busy.. please wait")
            stdscr.refresh()
            
            for x in range(0, len(sample)):
                sample_complete=False
                for finsamp in samples_complete:
                    if sample[x] == finsamp:
                        sample_complete=True
            if sample_complete:
                continue
            backgroundsamples=False
            for bs in samples_inbackground:
                if bs == sample[x]:
                    backgroundsamples=True
            useskimtmp=useskim
            ismctmp=  len(sample[x])>1
            if sample[x] == "H_v2" or sample[x] == "H_v3":
                ismctmp= False

            sampletmp=sample[x]
            runnptmp=runnp
            runcftmp=runcf
            channeltmp=channel
            if backgroundsamples:
                jobid_exists=True
                while jobid_exists:
                    path_job_check=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] + "jobid.txt"
                    if  os.path.exists(path_job_check):
                        jobid_exists=False
                file_job_check=open(path_job_check ,"r")
                jobid1=0
                jobid2=0
                njobs_finished=0.
                nrunning=0.
                nqueue=0.
                njobs_in_total=0.  ### should equal three above                                                                                                                                                                                                       
                path_clust_check=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] + "clust.txt"
                os.system("qstat -u " + os.getenv("USER") + " > " +  path_clust_check)
                ijob=0
                for sline in file_job_check:
                    ssline = sline.split()
                    if len(ssline) < 1:
                        continue
                    line=ssline[0]
                    if jobid1 == 0:
                        jobid1=line
                    jobid2=line
                    ijob=ijob+1.
                    njobs_in_total=njobs_in_total+1.
                    file_clust_check=open(path_clust_check,"r")
                    job_inqueue=False
                    for cline in file_clust_check:
                        splitline  = cline.split()
                        if len(splitline) < 6:
                            continue
                        if not os.getenv("USER")  in cline:
                            continue
                        if line == splitline[0]:
                            job_inqueue=True
                            if splitline[4] == "r":
                                nrunning = nrunning + 1.
                            else:
                                nqueue = nqueue + 1.
                    if not job_inqueue:
                        njobs_finished=njobs_finished+1.
                    file_clust_check.close()
                path_clust_check=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] + "clust.txt"
                os.system("rm " + path_clust_check)
                if njobs_in_total == 0:
                    continue
                file_job_check.close()
                if njobs_finished == njobs_in_total:
                    samples_complete.append(sample[x])
                nrunning = nrunning + njobs_finished
                nrun_per = round( (nrunning/float(njobs_in_total)), 2)
                nfin_per = round( (njobs_finished/float(njobs_in_total)), 2)
                nscreen_run= "["
                nscreen_fin= "["
                for xsr in range(0, 9):
                    if nrun_per > float(xsr)* 0.1:
                        nscreen_run+="#"
                    else:
                        nscreen_run+="-"
                if nrun_per == 1.:
                    nscreen_run+="#"
                else:
                    nscreen_run+="-"
                nscreen_run+="]"
                for xsr in range(0, 9):
                    if nfin_per > float(xsr)* 0.1:
                        nscreen_fin+="#"
                    else:
                        nscreen_fin+="-"
                if nfin_per == 1.:
                    nscreen_fin+="#"
                else:
                    nscreen_fin+="-"
                nscreen_fin+="]"
                
                if nrun_per == 0.:
                    stdscr.addstr(int(x)+istatus_message, box_shift,  str(int(x+1)), curses.A_NORMAL )
                    stdscr.addstr(int(x)+istatus_message, summary_block0, "| " +  cycle, curses.A_NORMAL )
                    if isMC:
                        stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + sample[x], curses.A_NORMAL )
                    else:
                        stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + channel+"_"+sample[x], curses.A_NORMAL )
                    stdscr.addstr(int(x)+istatus_message, summary_block2, "| " + skim_print, curses.A_NORMAL )
                    stdscr.addstr(int(x)+istatus_message, summary_block3,"| Submitted " , curses.A_NORMAL)
                    stdscr.refresh()
                else:
                    stdscr.addstr(int(x)+istatus_message, box_shift,  str(int(x+1)), curses.A_BOLD )
                    stdscr.addstr(int(x)+istatus_message, summary_block0,  "| " + cycle, curses.A_BOLD )
                    if isMC:
                        stdscr.addstr(int(x)+istatus_message, summary_block1,  "| " + sample[x], curses.A_BOLD )
                    else:
                        stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + channel+"_"+sample[x], curses.A_NORMAL )

                    stdscr.addstr(int(x)+istatus_message, summary_block2,  "| " + skim_print, curses.A_BOLD )
                    stdscr.addstr(int(x)+istatus_message, summary_block3,"| Running    " , curses.A_BLINK)
                    if not jobid1 == 0:
                        nblanks=2*(8-len(str(jobid2)))
                        stdscr.addstr(int(x)+istatus_message, summary_block4,"| " + str(jobid1) + "-" + str(jobid2) +  " "*nblanks +   "| " + str(1+int(jobid2)-int(jobid1)) , curses.A_BOLD)
                        stdscr.addstr(int(x)+istatus_message, summary_block5,"| " + printedqueue , curses.A_BOLD)
                                            
                    else:
                        stdscr.addstr(int(x)+istatus_message, summary_block4,"|" + " "*18 + "|"  , curses.A_BOLD)
                        stdscr.addstr(int(x)+istatus_message, summary_block5,"| " + printedqueue  , curses.A_BOLD)

                    stdscr.addstr(int(x)+istatus_message, summary_block6,  "| Running   " + nscreen_run + " " + str(100*nrun_per)+ "%  Complete" + nscreen_fin + " " + str(100*nfin_per) + "%",curses.A_BOLD)
                    stdscr.addstr(int(x)+istatus_message, summary_block7 ,"|    ",curses.A_BOLD)
                    stdscr.refresh()
                if nfin_per == 1.:
                    stdscr.addstr(int(x)+istatus_message, box_shift,  str(int(x+1)), curses.A_DIM )
                    stdscr.addstr(int(x)+istatus_message, summary_block0,  "| " + cycle, curses.A_DIM )
                    if isMC:
                        stdscr.addstr(int(x)+istatus_message, summary_block1,  "| " + sample[x], curses.A_DIM )
                    else:
                        stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + channel+"_"+sample[x], curses.A_NORMAL )
                    stdscr.addstr(int(x)+istatus_message, summary_block2,  "| " + skim_print, curses.A_DIM )
                    stdscr.addstr(int(x)+istatus_message, summary_block3,"| Complete " , curses.A_DIM)
                    if run_in_bkg:
                        del output_bkg[:]
                        for out_x in range(1,winx):
                            output_bkg.append(stdscr.instr(out_x, 0))
                            UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")
                    job_time = time.time() - start_time
                    if not jobid1 == 0:
                        nblanks=2*(8-len(str(jobid2)))
                        stdscr.addstr(int(x)+istatus_message, summary_block4,"| " + str(jobid1) + "-" + str(jobid2)  + " "*nblanks +   "| " + str(1+int(jobid2)-int(jobid1)), curses.A_DIM)
                        stdscr.addstr(int(x)+istatus_message, summary_block5,"| " + printedqueue  , curses.A_DIM)

                    else:
                        stdscr.addstr(int(x)+istatus_message, summary_block6,"| "  , curses.A_DIM)
                        
                    stdscr.refresh()
                    
        else:
            os.system("rm -r "+an_jonpre+"CAT_SKTreeOutput/" + os.getenv("USER") +"/"+tagger)
            stdscr.addstr(list4+1, box_shift,  "Submitting sample to queue  " + s)
            stdscr.refresh()
            
            checkqueue=False
            if not os.path.exists(an_jonpre+"/CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) +"/" + tagger):
                os.system("mkdir  "+an_jonpre+"CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) +"/" + tagger)

            stdscr.addstr(list3 + 1+int(isample), box_shift,  "Running " + s + " in background: terminal output sent to "+an_jonpre+"CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" +  str(tagger) + "/" + tagger + "/" + s + ".txt")
            stdscr.refresh()
            stdscr.addstr(int(isample)+istatus_message, box_shift,  str(int(isample+1)) )
            stdscr.addstr(int(isample)+istatus_message, summary_block0,  "| " + cycle  )
            stdscr.addstr(int(isample)+istatus_message, summary_block1,  "| " + s  )
            stdscr.addstr(int(isample)+istatus_message, summary_block2,  "| " + skim_print  )
            stdscr.addstr(int(isample)+istatus_message, summary_block3,"| In Queue " , curses.A_NORMAL)
            stdscr.addstr(int(isample)+istatus_message, summary_block4,"|" + " "*18 + "|"  , curses.A_NORMAL)
            stdscr.addstr(int(isample)+istatus_message, summary_block5,"|    ")
                        
            stdscr.addstr(int(isample)+istatus_message, summary_block6,"| Running[---------]   Complete[----------]   ")
            stdscr.addstr(int(isample)+istatus_message, summary_block7,"|    ")
            stdscr.refresh()
            os.system(command2_background)
            samples_inbackground.append(s)
            if run_in_bkg:
                del output_bkg[:]
                for out_x in range(1,winx):
                    output_bkg.append(stdscr.instr(out_x, 0))
                    UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")

stdscr.refresh()
StillRunning=True
while StillRunning:
    StillRunning=False

    if run_in_bkg:
        del output_bkg[:]
        for out_x in range(1,winx):
            output_bkg.append(stdscr.instr(out_x, 0))
        UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")

    for x in range(0, len(sample)):
        sample_complete=False
        for finsamp in samples_complete:
            if sample[x] == finsamp:
                sample_complete=True
        if sample_complete:
            continue
        ## Need to sleep here so that cursor is not moving on screen too fast
        time.sleep(2.)
        
        useskimtmp=useskim
        ismctmp=  len(sample[x])>1
        if sample[x] == "H_v2" or sample[x] == "H_v3":
            ismctmp= False

        sampletmp=sample[x]
        runnptmp=runnp
        runcftmp=runcf
        channeltmp=channel

        if run_in_bkg:
            del output_bkg[:]
            for out_x in range(1,winx):
                output_bkg.append(stdscr.instr(out_x, 0))
            UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")

        backgroundsamples=False
        for bs in samples_inbackground:
            if bs == sample[x]:
                backgroundsamples=True
        if backgroundsamples:
            jobid_exists=True
            while jobid_exists:
                path_job_check=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] + "jobid.txt"
                if  os.path.exists(path_job_check):
                    jobid_exists=False

            #CheckRunningStatus(an_jonpre+"/CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) +"/" + tagger + "/" + sample[x] +".txt")


            file_job_check=open(path_job_check ,"r")
            njobs_finished=0.
            nrunning=0.
            nqueue=0.
            njobs_in_total=0.
            path_clust_check=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] + "clust.txt"
            os.system("qstat -u " + os.getenv("USER") + " > " +  path_clust_check)
            ijob=0
            jobid1=0
            jobid2=0
            for sline in file_job_check:
                ssline = sline.split()
                if len(ssline) < 1:
                    continue
                line=ssline[0]
                if jobid1==0:
                    jobid1=line
                jobid2=line
                ijob=ijob+1.
                njobs_in_total=njobs_in_total+1.
                file_clust_check=open(path_clust_check,"r")
                job_in_queue=False
                for cline in file_clust_check:
                    splitline  = cline.split()
                    if len(splitline) < 6:
                        continue
                    if not os.getenv("USER")  in cline:
                        continue
                    if line == splitline[0]:
                        job_in_queue=True
                        if splitline[4] == "r":
                            nrunning = nrunning + 1.
                        else:
                            nqueue = nqueue + 1.
                if not job_in_queue:
                    njobs_finished=njobs_finished+1.
                    
                file_clust_check.close()
            path_clust_check=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] + "clust.txt"
            os.system("rm " + path_clust_check)
        

            if njobs_in_total         == 0.:
                continue
            file_job_check.close()
            if njobs_finished == njobs_in_total:
                samples_complete.append(sample[x])
            else:
                StillRunning=True

            nrunning = nrunning + njobs_finished    
            nrun_per = round( (nrunning/float(njobs_in_total)), 2)
            nfin_per = round( (njobs_finished/float(njobs_in_total)), 2)
            
            nscreen_run= "["
            nscreen_fin= "["
            for xsr in range(0, 9):
                if nrun_per > float(xsr)* 0.1:
                    nscreen_run+="#"
                else:
                    nscreen_run+="-"
            if nrun_per == 1.:
                nscreen_run+="#"
            else:
                nscreen_run+="-"
            nscreen_run+="]"
            for xsr in range(0, 9):
                if nfin_per > float(xsr)* 0.1:
                    nscreen_fin+="#"
                else:
                    nscreen_fin+="-"
            if nfin_per == 1.:
                nscreen_fin+="#"
            else:
                nscreen_fin+="-"
            nscreen_fin+="]"
            
            if nrun_per == 0.:
                stdscr.addstr(int(x)+istatus_message, box_shift,  str(int(x+1)), curses.A_NORMAL )
                stdscr.addstr(int(x)+istatus_message, summary_block0,  "| " + cycle, curses.A_NORMAL )
                if isMC:
                    stdscr.addstr(int(x)+istatus_message, summary_block1, "| " +  sample[x], curses.A_NORMAL )
                else:
                    stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + channel+"_"+sample[x], curses.A_NORMAL )
                stdscr.addstr(int(x)+istatus_message, summary_block2, "| " +  skim_print, curses.A_NORMAL )
                stdscr.addstr(int(x)+istatus_message, summary_block3,"| Submitted " , curses.A_NORMAL)
                stdscr.refresh()
            else:
                stdscr.addstr(int(x)+istatus_message, box_shift,  str(int(x+1)), curses.A_BOLD )
                stdscr.addstr(int(x)+istatus_message, summary_block0,  "| " + cycle, curses.A_BOLD )
                if isMC:
                    stdscr.addstr(int(x)+istatus_message, summary_block1,  "| " + sample[x], curses.A_BOLD )
                else:
                    stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + channel+"_"+sample[x], curses.A_NORMAL )
                stdscr.addstr(int(x)+istatus_message, summary_block2,  "| " + skim_print, curses.A_BOLD )
                stdscr.addstr(int(x)+istatus_message, summary_block3,"| Running    " , curses.A_BLINK)
                if not jobid1 == 0:
                    nblanks=2*(8-len(str(jobid2)))
                    stdscr.addstr(int(x)+istatus_message, summary_block4,"| " + str(jobid1) + "-" + str(jobid2)  + " "*nblanks +   "| " + str(1+int(jobid2)-int(jobid1)), curses.A_BOLD)
                    stdscr.addstr(int(x)+istatus_message, summary_block5,"| " + printedqueue  , curses.A_BOLD)

                else:
                    stdscr.addstr(int(x)+istatus_message, summary_block4,"|" + " "*18 + "|" , curses.A_BOLD)
                    stdscr.addstr(int(x)+istatus_message, summary_block5,"| " + printedqueue  , curses.A_BOLD)

                stdscr.addstr(int(x)+istatus_message, summary_block6,  "| Running   " + nscreen_run + " " + str(100*nrun_per)+ "%  Complete" + nscreen_fin + " " + str(100*nfin_per) + "%",curses.A_BOLD)
                
                stdscr.addstr(int(x)+istatus_message, summary_block7 ,"|    ",curses.A_BOLD)
                stdscr.refresh()
            if nfin_per == 1.:
                if run_in_bkg:
                    del output_bkg[:]
                    for out_x in range(1,winx):
                        output_bkg.append(stdscr.instr(out_x, 0))
                        UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")

                stdscr.addstr(int(x)+istatus_message, box_shift,  str(int(x+1)), curses.A_DIM )
                stdscr.addstr(int(x)+istatus_message, summary_block0,  "| " + cycle, curses.A_DIM )
                if isMC:
                    stdscr.addstr(int(x)+istatus_message, summary_block1,  "| " + sample[x], curses.A_DIM )
                else:
                    stdscr.addstr(int(x)+istatus_message, summary_block1, "| " + channel+"_"+sample[x], curses.A_NORMAL )
                stdscr.addstr(int(x)+istatus_message, summary_block2,  "| " + skim_print, curses.A_DIM )
                stdscr.addstr(int(x)+istatus_message, summary_block3,"| Complete " , curses.A_DIM)
                job_time = time.time() - start_time
                if not jobid1 == 0:
                    nblanks=2*(8-len(str(jobid2)))
                    stdscr.addstr(int(x)+istatus_message, summary_block4,"| " + str(jobid1) + "-" + str(jobid2)  + " "*nblanks +   "| " + str(1+int(jobid2)-int(jobid1)), curses.A_DIM)
                    stdscr.addstr(int(x)+istatus_message, summary_block5,"| " + printedqueue, curses.A_BOLD)

                else:
                    stdscr.addstr(int(x)+istatus_message, summary_block4,"|" + " "*18 + "|"  , curses.A_DIM)
                    stdscr.addstr(int(x)+istatus_message, summary_block5,"| " + printedqueue , curses.A_BOLD)

                stdscr.addstr(int(x)+istatus_message, summary_block6,  "| Running   " + nscreen_run + " " + str(100*nrun_per)+ "%  Complete" + nscreen_fin + " " + str(100*nfin_per) + "%",curses.A_DIM)
                stdscr.addstr(int(x)+istatus_message, summary_block7 ,"|    ",curses.A_DIM)
                stdscr.refresh()
                path_job=path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser()  + "/" + str(tagger)+ "/statlog_time_"+sample[x] + tagger + ".txt"
                ismerging=True
                while not os.path.exists(path_job):
                    if ismerging:
                        stdscr.addstr(2+list2+int(x), box_shift , str(int(x+1)),curses.A_DIM)
                        stdscr.addstr(2+list2+int(x), summary2_block1 ,"| MERGING OUTPUT ", curses.A_BLINK)
                        stdscr.refresh()
                        ismerging=False
                        
                    check_crash_stat=False
                    if check_crash_stat:
                        check_crashfile=an_jonpre+"/CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) +"/" + tagger + "/" + sample[x] +".txt"    
                        if not os.path.exists(check_crashfile):
                            continue
                        file_check_crashfile=open(check_crashfile,"r")
                        crash_in_job=False
                        crashlog=[]
                        for line in file_check_crashfile:
                            if " No such file or directory:" in line:
                                crash_in_job=True
                                crashlog.append(line)
                        file_check_crashfile.close()
                        if crash_in_job:
                            crash_log= an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] +"_crash/crashlog.txt"
                            if os.path.exists(crash_log):
                                file_crash=open(crash_log,"a")
                                for linex in  crashlog:
                                    file_crash.write(linex)
                                file_crash.close()
                            else:
                                file_crash=open(crash_log,"w")
                                for linex in  crashlog:
                                    file_crash.write(linex)
                                file_crash.close()
                            break


                    continue

                crash_log= an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/" + sample[x] +"_crash/crashlog.txt"
                if os.path.exists(crash_log):
                    jobidcrash.append(int(x))
                    stdscr.addstr(2+list2+int(x), box_shift , str(int(x+1)),curses.A_DIM)
                    stdscr.addstr(2+list2+int(x), summary2_block1 ,"| JOB CRASHED      ", curses.A_BLINK)
                    logpath=GetLogFilePath(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle)
                    outfilepath=str(Finaloutputdir)  +GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,False)
                    fileoutputlist.append(outfilepath)

                    if "SKTreeMaker" in cycle:
                        outfilepath=GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,False)

                    stdscr.addstr(list3b + 1+int(x), box_shift,  str(1+int(x)) +": Log files for " + sample[x] + " found at " + str(logpath))
                    stdscr.addstr(list3c + 1+int(x), box_shift,   str(1+int(x)) +":###CRASH###: No OutputFile for " + sample[x])
                    crash_output.append(logpath)
                    crash_outputjob.append(jobid1)
                    stdscr.refresh()
                    if run_in_bkg:
                        del output_bkg[:]
                        for out_x in range(1,winx):
                            output_bkg.append(stdscr.instr(out_x, 0))
                            UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")
                else:
                    stdscr.addstr(2+list2+int(x), box_shift , str(int(x+1)),curses.A_DIM)
                    stdscr.addstr(2+list2+int(x), summary2_block0 ,"| " + str(round(job_time,2)) + "[s]",curses.A_DIM)
                    stdscr.refresh()
                
                    logpath=GetLogFilePath(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle)
                    outfilepath=str(Finaloutputdir) + GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,False)
                    fileoutputlist.append(outfilepath)
                    if "SKTreeMaker" in cycle:
                        outfilepath=GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,False)
                        

                    stdscr.addstr(list3b + 1+int(x), box_shift,  str(1+int(x)) +": Log files for " + sample[x] + " found at " + str(logpath))
                    stdscr.addstr(list3c + 1+int(x), box_shift,   str(1+int(x)) +": OutputFile for " + sample[x] + " = " + str(outfilepath))

                    filesize=0.
                    memoryusage_p=0.
                    memoryusage_v=0.
                    outputfile_size=0.
                    first_job_time=0.
                    last_job_time=0.
                    file_job=open(path_job,"r")
                    for line in file_job:
                        splitline  = line.split()
                        if len(splitline) < 2:
                            continue
                        if "job_time"  == splitline[0]:
                            first_job_time=float(splitline[1])
                        if "last_job_time"  == splitline[0]:
                            last_job_time=float(splitline[1])
                        if "outputfile_size" == splitline[0]:
                            outputfile_size=splitline[1]
                        if "memoryusage_p" ==  splitline[0]:
                            memoryusage_p=splitline[1]
                        if "memoryusage_v" ==  splitline[0]:
                            memoryusage_v=splitline[1]
                    file_job.close()
                    stdscr.addstr(2+list2+int(x), summary2_block1 ,"| " + str(round(first_job_time,2)) + "-" +  str(round(last_job_time,2)) + "[s]     ",curses.A_DIM) 


                    if LargeMemory(memoryusage_p):
                        output_warning.append("WARNING:: Physical memory of job " + sample[x] + "  is excess of " + str(large_memory_check) + " MB.")
                    if LargeMemory(memoryusage_v):
                        output_warning.append("WARNING:: Virtual memory of job " + sample[x] + "  is excess of " + str(large_memory_check) + " MB.")

                    if isMC:
                        if LargePhysicalMemoryIncrease(memoryusage_p,sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Physical memory of job " + sample[x] + "  is " + RoundMemory(memoryusage_p) + " compared to previous jobs " + GetPhysicalMemoryUsage(sample[x], cycle,tagger,useskim))
                        if LargeVirtualMemoryIncrease(memoryusage_v,sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Virtual memory of job " + sample[x] + "  is " + RoundMemory(memoryusage_v) + " compared to previous jobs " + GetVirtualMemoryUsage(sample[x], cycle,tagger,useskim))
                        if LargeFileSize(outputfile_size):
                            output_warning.append("WARNING:: Size of output rootfile for " + sample[x] + " is " + str(outputfile_size) + " greater than " + str(large_file_size) + " MB")
                        if LargeFileSizeIncrease(outputfile_size,sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Size of output rootfile for " + sample[x] + " is " + str(outputfile_size) + " compared to previous jobs " + GetFileSize(sample[x], cycle,tagger,useskim))
                        if TimeIncrease(first_job_time,sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Job time per input file for " +  sample[x] +"  increased from " + str(first_job_time) + " compared to previous job time " + GetTime(sample[x], cycle,tagger,useskim))
                    else:
                        if LargePhysicalMemoryIncrease(memoryusage_p,channel + "_" +sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Physical memory of job " + sample[x] + "  is " + RoundMemory(memoryusage_p) + " compared to previous jobs " + GetPhysicalMemoryUsage(channel + "_" +sample[x], cycle,tagger,useskim))
                        if LargeVirtualMemoryIncrease(memoryusage_v,channel + "_" +sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Virtual memory of job " + sample[x] + "  is " + RoundMemory(memoryusage_v) + " compared to previous jobs " + GetVirtualMemoryUsage(channel + "_" +sample[x], cycle,tagger,useskim))
                        if LargeFileSize(outputfile_size):
                            output_warning.append("WARNING:: Size of output rootfile for " + sample[x] + " is " + str(outputfile_size) + " greater than " + str(large_file_size) + " MB")
                        if LargeFileSizeIncrease(outputfile_size,channel + "_" +sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Size of output rootfile for " + sample[x] + " is " + str(outputfile_size) + " compared to previous jobs " + GetFileSize(channel + "_" +sample[x], cycle,tagger,useskim))
                        if TimeIncrease(first_job_time,channel + "_" +sample[x], cycle,tagger,useskim):
                            output_warning.append("WARNING:: Job time per input file " + sample[x] +" increased from " + str(first_job_time) + " compared to previous job time " + GetTime(channel + "_" +sample[x], cycle,tagger,useskim))
               

                    stdscr.addstr(2+list2+int(x), summary2_block2 ,"| " + RoundMemory(memoryusage_p),curses.A_DIM)
                    stdscr.addstr(2+list2+int(x), summary2_block3 ,"| " + RoundMemory(memoryusage_v),curses.A_DIM)
                    stdscr.addstr(2+list2+int(x), summary2_block4 ,"| " + str(outputfile_size),curses.A_DIM)
                    stdscr.addstr(2+list2+int(x), summary2_block5 ,"| ",curses.A_DIM)
                    stdscr.refresh()
                    if run_in_bkg:
                        del output_bkg[:]
                        for out_x in range(1,winx):
                            output_bkg.append(stdscr.instr(out_x, 0))
                            UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")
            stdscr.refresh()
            if run_in_bkg:
                del output_bkg[:]
                for out_x in range(1,winx):
                    output_bkg.append(stdscr.instr(out_x, 0))
                    UpdateOutput(output_bkg,an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt")

            
    stdscr.refresh()

#del output_bkg[:]
#for out_x in range(1,winx):
#    output_bkg.append(stdscr.instr(out_x, 0))


############## END OF JOB ################3

mypad_contents = []
for i in range(0, winx - remove_from_end):
    mypad_contents.append(stdscr.instr(i, 0))

redopad=False
if redopad:
    time.sleep(5.)
    for i in range(0, winx-remove_from_end):
        stdscr.addstr(i,0, mypad_contents[i])
        stdscr.refresh()
    time.sleep(5.)
if not run_in_bkg:    
    curses.echo()
    curses.nocbreak()
    curses.endwin()


iplotsample=-1
listofplots=[]

collist=[]
#fcol    870       
#zcol    64
#wcol    92
#wwcol   68
#zzcol   80
#wzcol   74
#ttvcol  88
#tcol    88 
#ttcol   97
collist.append(64) 
collist.append(92) 
collist.append(68) 

if quickdraw:
    plotallhist(fileoutputlist,an_jonpre+"/CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) +"/"+ str(tagger)  + "_hist.pdf", collist, tagger)
    if not DoSendEmail:
        os.system("display " + an_jonpre+"/CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) +"/"+ str(tagger)  + "_hist.pdf&")
    listofplots.append(an_jonpre+"/CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) +"/"+ str(tagger) + "_hist.pdf")

mergedoutfilepath=""
if not ismctmp:
    if runnptmp == "True":
        datadir=str(Finaloutputdir)
        datadir=datadir.replace("Fake/","")
        mergedoutfilepath=str(datadir) + GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,True)
    elif runcftmp == "True":
        datadir=str(Finaloutputdir)
        datadir=datadir.replace("CF/","")
        mergedoutfilepath=str(datadir) + GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,True)
    else:
        datadir=str(Finaloutputdir)
        datadir=datadir.replace("Data/","")
        mergedoutfilepath=str(datadir) + GetOutFileName(useskimtmp, ismctmp , sampletmp, runnptmp, runcftmp, channeltmp , cycle,True)

        
job_summary=[]
print "\n"
for i in range(0, winx-remove_from_end):
    if "Job Terminal Output:" in mypad_contents[i]:
        job_summary.append("*"*summary_block6+"\n")
        print "*"*summary_block6
        print  "Job Terminal Output:(will be deleted "+ future +")" + " " *20
        job_summary.append("Job Terminal Output:(will be deleted "+ future +")\n")
    elif  "Log Files:"  in mypad_contents[i]:
        print  "Log Files:(will be deleted "+ future_week +")" + " " *20
        job_summary.append("Log Files:(will be deleted "+ future_week +")\n")
    elif "Output Files:" in mypad_contents[i]:
        print mypad_contents[i]
        job_summary.append(mypad_contents[i]+"\n")
        if not ismctmp:
            if data_lumi== "ALL":
                print "Merged output = " +mergedoutfilepath
                job_summary.append("Merged output = " + mergedoutfilepath+"\n")
    else:
        print mypad_contents[i]
        job_summary.append(mypad_contents[i]+"\n")
    if "PostJob " in mypad_contents[i]:
        print "_"*summary2_block5    
        job_summary.append("_"*summary2_block5+"\n")
    elif "Cum.Process" in mypad_contents[i]:
        print "_"*summary2_block5
        job_summary.append("_"*summary2_block5+"\n")
    elif "Terminal" in mypad_contents[i]:
        print "_"*50
        job_summary.append("_"*50+"\n")
    elif "Job " in mypad_contents[i]:
        print "_"*summary_block6
        job_summary.append("_"*summary_block6+"\n")
    if "Log Files:"  in mypad_contents[i]:
        print "_"*40
        job_summary.append("_"*40+"\n")
    if "Output Files" in mypad_contents[i]:
        print "_"*40
        job_summary.append("_"*40+"\n")


path_stat_dir=path_jobpre+"/LQAnalyzer_rootfiles_for_analysis/CATAnalyzerStatistics/" + getpass.getuser()  + "/" + str(tagger)+ "/"
os.system("rm -r " + path_stat_dir)


end_job_time = time.time() - start_time

email_subject=""
if end_job_time > email_time_limit:
    DoSendEmail=True
    email_subject=email_subject+"Job time>600s "

if run_in_bkg:
    os.system("mv "+an_jonpre+"CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_bkg.txt "+an_jonpre+"CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/output_finished.txt")
if runningData:
    MergeData(runnp,runcf,data_lumi, Finaloutputdir, catversion, useskim, cycle, channel)

if len(output_warning) > 0:
    print "\n"
    print "The job produced the following WARNING MESSAGES:"
    job_summary.append( "\n")
    job_summary.append("The job produced the following WARNING MESSAGES:")
    DoSendEmail=True
    email_subject=email_subject+"Job WARNING "
    for s in output_warning:
        print s
        job_summary.append(s + "\n")


crashlog_printout=[]


if len(crash_output) > 0:
    print "\n"
    print "%"*summary_block6
    print "%"*summary_block6
    print "Crash Job Error Message:"
    print "_"*summary_block6
    job_summary.append("\n")
    job_summary.append("%"*summary_block6+"\n")
    job_summary.append("%"*summary_block6+"\n")
    job_summary.append("Crash Job Error Message:\n")
    job_summary.append("_"*summary_block6+"\n")
    jobids_failed=""
    for x in jobidcrash:
        jobids_failed=jobids_failed+str(x+1)+","
    jobids_failed= jobids_failed[:-1]    
    print "Job IDs " + jobids_failed + " failed with following error:"
    print "#"*summary_block6
    job_summary.append("Job IDs " + jobids_failed + " failed with following error:\n")
    job_summary.append("#"*summary_block6+"\n")
    Crash_Printed=False
    for c in crash_output:
        if Crash_Printed:
            continue
        index_crash=0
        found_crash=False
        njobs_tried=-1
        while not found_crash:
            njobs_tried=njobs_tried+1
            errlogpath=c
            errlogpath=errlogpath.replace("Job_*","Job_"+str(njobs_tried+1) +".e"+str(int(crash_outputjob[index_crash])+njobs_tried))
            while not os.path.exists(errlogpath):
                index_crash=index_crash+1
                if index_crash == len(crash_outputjob):
                    break
                errlogpath=errlogpath.replace("Job_*","Job_" + str(njobs_tried+1) +".e"+str(int(crash_outputjob[index_crash])+njobs_tried))

            if not os.path.exists(errlogpath):
                found_crash=True
                print "Error in locating log file " + errlogpath + ". Could not print out error message from batch jobs"
                job_summary.append("Error in locating log file " + errlogpath + ". Could not print out error message from batch jobs\n")
            else:
                file_read_err = open(errlogpath,"r")
                for rline in file_read_err:
                    print rline    
                    job_summary.append(rline+"\n")
                    crashlog_printout.append(rline)
                    found_crash=True
                    Crash_Printed=True
        print "#"*summary_block6
        job_summary.append("#"*summary_block6+"\n")
    list_crash=CrashHelper(crashlog_printout)
    for linecrash in list_crash:
        job_summary.append(linecrash+"\n")
    
    print " "*summary_block6
    print "Run the following command to help debug job error: command runs crashed job in terminal instead of on batch machine."  
    print "sktree -a " + cycle + " -i " + sample[jobidcrash[0]] + " -s " + useskim + " -d DEBUG -n 1"
    job_summary.append(" "*summary_block6+"\n")
    job_summary.append("Run the following command to help debug job error: command runs crashed job in terminal instead of on batch machine.\n")
    job_summary.append("sktree -a " + cycle + " -i " + sample[jobidcrash[0]] + " -s " + useskim + " -d DEBUG -n 1 \n")
    email_subject=email_subject+"Job CRASH "
    DoSendEmail=True

if DoSendEmail:
    SendEmail(job_summary,tagger,email_subject,useremail,quickdraw, listofplots)


##### CODE WRITEEN TO REMOVE DIRECTORY... THIS WILL BE KEPT UNLESS USER WISHES TO SET RDIR+TRUE

remdir=False

if remdir:
    njobs_in_total=0.
    path_clust_check2=an_jonpre+"/CAT_SKTreeOutput/" + os.getenv("USER")  + "/CLUSTERLOG" + str(tagger)+ "/clusterjobs.txt"
    os.system("qstat -u " + os.getenv("USER") + " > " +  path_clust_check2)
    file_clust_check2=open(path_clust_check2,"r")
    for sline in file_clust_check2:
        ssline = sline.split()
        if len(ssline) < 1:
            continue
        line=ssline[0]
        njobs_in_total=njobs_in_total+1.

    file_clust_check2.close()
    
    if njobs_in_total == 0 :    
        os.system("rm -r "+an_jonpre+"/CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG*")
    else:
        os.system("rm -r "+an_jonpre+"CAT_SKTreeOutput/"+os.getenv("USER")+"/CLUSTERLOG" + str(tagger) )
