import pickle
import numpy as np
from datetime import datetime
from multiprocessing import Pool
import pandas as pd
import statsmodels.api as sm
import sys
import re
import warnings
from scipy.stats import ks_2samp
import glob
from pathlib import Path
import os
import json
import pickle


#This script checks if U, v0,v1,v2, eta_H values reach equilibrium and writes summary file of dist
#This file checks pkl files

argErrCode=2
sameErrCode=3
missingErrCode=4
if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit(argErrCode)

# print("entering")
jsonFromSummaryLast=json.loads(sys.argv[1])
jsonDataFromConf=json.loads(sys.argv[2])
# print(jsonFromSummaryLast)
TDirRoot=jsonFromSummaryLast["TDirRoot"]
U_dist_dataDir=jsonFromSummaryLast["U_dist_dataDir"]
effective_data_num_required=int(jsonDataFromConf["effective_data_num_required"])
N=int(jsonDataFromConf["N"])

summary_U_distFile=TDirRoot+"/summary_U_dist.txt"


lastFileNum=10
def sort_data_files_by_flushEnd(oneDir):
    dataFilesAll=[]
    flushEndAll=[]
    for oneDataFile in glob.glob(oneDir+"/flushEnd*.pkl"):
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"flushEnd(\d+)",oneDataFile)
        if matchEnd:
            indTmp=int(matchEnd.group(1))
            flushEndAll.append(indTmp)
    endInds=np.argsort(flushEndAll)
    sortedDataFiles=[dataFilesAll[i] for i in endInds]
    return sortedDataFiles



def parseSummaryU_Dist():
    startingFileInd=-1
    startingVecPosition=-1

    summaryFileExists=os.path.isfile(summary_U_distFile)
    if summaryFileExists==False:
        return startingFileInd,startingVecPosition
    with open(summary_U_distFile,"r") as fptr:
        lines=fptr.readlines()

    for oneLine in lines:
        #match startingFileInd
        matchStartingFileInd=re.search(r"startingFileInd=(\d+)",oneLine)
        if matchStartingFileInd:
            startingFileInd=int(matchStartingFileInd.group(1))

        #match startingVecPosition
        matchStartingVecPosition=re.search(r"startingVecPosition=(\d+)",oneLine)
        if matchStartingVecPosition:
            startingVecPosition=int(matchStartingVecPosition.group(1))

    return startingFileInd, startingVecPosition



def auto_corrForOneVec(vec):
    """

    :param colVec: a vector of data
    :return:
    """
    same=False
    eps=1e-2
    NLags=int(len(vec)*1/4)

    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVec=sm.tsa.acf(vec,nlags=NLags)
    except Warning as w:
        same=True

    acfOfVecAbs=np.abs(acfOfVec)
    minAutc=np.min(acfOfVecAbs)

    lagVal=-1
    if minAutc<=eps:
        lagVal=np.where(acfOfVecAbs<=eps)[0][0]
    # np.savetxt("autc.txt",acfOfVecAbs[lagVal:],delimiter=',')
    return same,lagVal

def ksTestOneColumn(vec,lag):
    """

    :param vec: a vector of data
    :param lag: auto-correlation length
    :return:
    """
    vecSelected=vec[::lag]

    lengthTmp=len(vecSelected)
    if lengthTmp%2==1:
        lengthTmp-=1
    lenPart=int(lengthTmp/2)

    vecToCompute=vecSelected[-lengthTmp:]

    #ks test
    selectedVecPart0=vecToCompute[:lenPart]
    selectedVecPart1=vecToCompute[lenPart:]
    result=ks_2samp(selectedVecPart0,selectedVecPart1)

    return result.pvalue,result.statistic, lenPart*2



def checkUDataFilesForOneT(UData_dir):
    U_sortedDataFilesToRead=sort_data_files_by_flushEnd(UData_dir)
    if len(U_sortedDataFilesToRead)==0:
        print("no data for U.")
        exit(0)

    startingFileInd,startingVecPosition=parseSummaryU_Dist()

    if startingFileInd<0:
        #we guess that the equilibrium starts at this file
        startingFileInd=len(U_sortedDataFilesToRead)-lastFileNum
    startingFileName=U_sortedDataFilesToRead[startingFileInd]
    # print("U, startingFileName="+str(startingFileName))
    with open(startingFileName,"rb") as fptr:
        inArrStart=pickle.load(fptr)

    in_nRowStart=len(inArrStart)
    if startingVecPosition<0:
        #we guess equilibrium starts at this position
        startingVecPosition=0
    arr=inArrStart[startingVecPosition:]

    #read the rest of the pkl files
    # print(U_sortedDataFilesToRead[(startingFileInd+1):])
    for pkl_file in U_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            inArr=pickle.load(fptr)
            # print("len(inArr)="+str(len(inArr)))
        arr=np.append(arr,inArr)

    avg_Uarr=arr

    sameUTmp,lagUTmp=auto_corrForOneVec(avg_Uarr)

    #if one lag==-1, then the auto-correlation is too large

    if sameUTmp==True or lagUTmp==-1:
        return [sameUTmp,lagUTmp,-1,-1,-1,-1,-1]

    pUTmp,statUTmp,lengthUTmp=ksTestOneColumn(avg_Uarr,lagUTmp)
    numDataPoints=lengthUTmp

    return [sameUTmp,lagUTmp,pUTmp,statUTmp,numDataPoints,startingFileInd,startingVecPosition]



def check_oneDistDataFilesForOneT(v0Dir,v1Dir,v2Dir,N0,N1,N2,eta_HDir):

    v0_sortedDataFilesToRead=sort_data_files_by_flushEnd(v0Dir)
    v1_sortedDataFilesToRead=sort_data_files_by_flushEnd(v1Dir)
    v2_sortedDataFilesToRead=sort_data_files_by_flushEnd(v2Dir)

    if len(v0_sortedDataFilesToRead)!= len(v1_sortedDataFilesToRead):
        print("data missing.")
        exit(missingErrCode)

    if len(v1_sortedDataFilesToRead)!= len(v2_sortedDataFilesToRead):
        print("data missing.")
        exit(missingErrCode)
    startingFileInd,startingVecPosition=parseSummaryU_Dist()

    if startingFileInd<0:
        #we guess that the equilibrium starts at this file
        startingFileInd=len(v0_sortedDataFilesToRead)-lastFileNum

    v0StartingFileName=v0_sortedDataFilesToRead[startingFileInd]
    v1StartingFileName=v1_sortedDataFilesToRead[startingFileInd]
    v2StartingFileName=v2_sortedDataFilesToRead[startingFileInd]

    # print("v0StartingFileName="+str(v0StartingFileName))
    # print("v1StartingFileName="+str(v1StartingFileName))
    # print("v2StartingFileName="+str(v2StartingFileName))

    #read v0
    with open(v0StartingFileName,"rb") as fptr:
        v0_inArrStart=pickle.load(fptr)

    v0Arr=v0_inArrStart
    #read the rest of v0 pkl files
    for pkl_file in v0_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            v0_inArr=pickle.load(fptr)
            v0Arr=np.append(v0Arr,v0_inArr)

    #read v1
    with open(v1StartingFileName,"rb") as fptr:
        v1_inArrStart=pickle.load(fptr)
    v1Arr=v1_inArrStart
    #read the rest of v1 pkl files
    for pkl_file in v1_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            v1_inArr=pickle.load(fptr)
            v1Arr=np.append(v1Arr,v1_inArr)

    #read v2
    with open(v2StartingFileName,"rb") as fptr:
        v2_inArrStart=pickle.load(fptr)

    v2Arr=v2_inArrStart
    #read the rest of v2 pkl files
    for pkl_file in v2_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            v2_inArr=pickle.load(fptr)
            v2Arr=np.append(v2Arr,v2_inArr)

    v0Arr=v0Arr.reshape((-1,5))
    v1Arr=v1Arr.reshape((-1,5))
    v2Arr=v2Arr.reshape((-1,5))

    print(v0Arr.shape)





UDataDir=U_dist_dataDir+"/U/"

sameVec=[]
lagVec=[]
pVec=[]
statVec=[]
numDataVec=[]
print("checking U")
sameUTmp,lagUTmp,pUTmp,statUTmp,numDataPointsU,startingFileInd,startingVecPosition=checkUDataFilesForOneT(UDataDir)

print("lagU="+str(lagUTmp))
sameVec.append(sameUTmp)
lagVec.append(lagUTmp)
pVec.append(pUTmp)
statVec.append(statUTmp)
numDataVec.append(numDataPointsU)


discrete_vals=list(range(0,N))
N0=np.random.choice(discrete_vals, size=1, replace=True)[0]
N1=np.random.choice(discrete_vals, size=1, replace=True)[0]
N2=np.random.choice(discrete_vals, size=1, replace=True)[0]

v0_dataDir=U_dist_dataDir+"/v0/"
v1_dataDir=U_dist_dataDir+"/v1/"
v2_dataDir=U_dist_dataDir+"/v2/"
eta_H_dataDir=U_dist_dataDir+"/eta_H/"

check_oneDistDataFilesForOneT(v0_dataDir,v1_dataDir,v2_dataDir,N0,N1,N2,eta_H_dataDir)