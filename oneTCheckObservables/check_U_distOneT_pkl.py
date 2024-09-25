import pickle
import numpy as np
from datetime import datetime
from multiprocessing import Pool
import pandas as pd
import statsmodels.api as sm
import sys
import re
import warnings

from oauthlib.uri_validate import pct_encoded
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
sweep_to_write=int(jsonDataFromConf["sweep_to_write"])
summary_U_distFile=TDirRoot+"/summary_U_dist.txt"
xi_Ba=0.2
xi_Ti=0.76
xi_O_parallel=-0.53
xi_O_perpendicular=-0.21

xiVec=np.array([xi_Ba,xi_Ti,xi_O_parallel,xi_O_perpendicular,xi_O_perpendicular])
lastFileNum=20
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
    eps=5e-2
    NLags=int(len(vec)*3/4)

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



def check_oneDistDataFilesForOneT(v0Dir,v1Dir,v2Dir,i_val,j_val,k_val):

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
        v0_inArrStart=np.array(pickle.load(fptr))

    v0Arr=v0_inArrStart.reshape((sweep_to_write,-1,5))
    #read the rest of v0 pkl files
    for pkl_file in v0_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            v0_inArr=np.array(pickle.load(fptr))
            v0_inArr=v0_inArr.reshape((sweep_to_write,-1,5))

            v0Arr=np.concatenate((v0Arr,v0_inArr),axis=0)



    #read v1
    with open(v1StartingFileName,"rb") as fptr:
        v1_inArrStart=np.array(pickle.load(fptr))

    v1Arr=v1_inArrStart.reshape((sweep_to_write,-1,5))
    #read the rest of v1 pkl files
    for pkl_file in v1_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            v1_inArr=np.array(pickle.load(fptr))
            v1_inArr=v1_inArr.reshape((sweep_to_write,-1,5))
            v1Arr=np.concatenate((v1Arr,v1_inArr),axis=0)



    #read v2
    with open(v2StartingFileName,"rb") as fptr:
        v2_inArrStart=np.array(pickle.load(fptr))

    v2Arr=v2_inArrStart.reshape((sweep_to_write,-1,5))
    #read the rest of v2 pkl files
    for pkl_file in v2_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            v2_inArr=np.array(pickle.load(fptr))
            v2_inArr=v2_inArr.reshape((sweep_to_write,-1,5))
            v2Arr=np.concatenate((v2Arr,v2_inArr),axis=0)

    ind=i_val*N*N+j_val*N+k_val

    v0_oneUnitCell=v0Arr[:,ind,:]
    v1_oneUnitCell=v1Arr[:,ind,:]
    v2_oneUnitCell=v2Arr[:,ind,:]

    u0_oneUnitCell=v0_oneUnitCell@xiVec
    u1_oneUnitCell=v1_oneUnitCell@xiVec
    u2_oneUnitCell=v2_oneUnitCell@xiVec
    print("u0_oneUnitCell[-20:]:\n"+str(u0_oneUnitCell[-20:]))
    print("u1_oneUnitCell[-20:]:\n"+str(u1_oneUnitCell[-20:]))
    print("u2_oneUnitCell[-20:]:\n"+str(u2_oneUnitCell[-20:]))

    print("mean(u0_oneUnitCell)="+str(np.mean(u0_oneUnitCell)))
    print("mean(u1_oneUnitCell)="+str(np.mean(u1_oneUnitCell)))
    print("mean(u2_oneUnitCell)="+str(np.mean(u2_oneUnitCell)))


    # sameVec_v0v1v2=[]
    # lagVec_v0v1v2=[]

    sameVec_u0u1u2=[]
    lagVec_u0u1u2=[]


    _,nCol=v0_oneUnitCell.shape

    sameTmp_u0,lagTmp_u0=auto_corrForOneVec(u0_oneUnitCell)
    sameTmp_u1,lagTmp_u1=auto_corrForOneVec(u1_oneUnitCell)
    sameTmp_u2,lagTmp_u2=auto_corrForOneVec(u2_oneUnitCell)

    sameVec_u0u1u2.append(sameTmp_u0)
    sameVec_u0u1u2.append(sameTmp_u1)
    sameVec_u0u1u2.append(sameTmp_u2)

    lagVec_u0u1u2.append(lagTmp_u0)
    lagVec_u0u1u2.append(lagTmp_u1)
    lagVec_u0u1u2.append(lagTmp_u2)


    # print("lagTmp_u0="+str(lagTmp_u0))
    # print("lagTmp_u1="+str(lagTmp_u1))
    # print("lagTmp_u2="+str(lagTmp_u2))

    # for j in range(0,nCol):
    #
    #     #v0
    #     sameTmp,lagTmp=auto_corrForOneVec(v0_oneUnitCell[:,j])
    #     sameVec_v0v1v2.append(sameTmp)
    #     lagVec_v0v1v2.append(lagTmp)
    #
    #     #v1
    #     sameTmp,lagTmp=auto_corrForOneVec(v1_oneUnitCell[:,j])
    #     sameVec_v0v1v2.append(sameTmp)
    #     lagVec_v0v1v2.append(lagTmp)
    #
    #
    #     #v2
    #     sameTmp,lagTmp=auto_corrForOneVec(v2_oneUnitCell[:,j])
    #     sameVec_v0v1v2.append(sameTmp)
    #     lagVec_v0v1v2.append(lagTmp)
    # print(lagVec_v0v1v2)
    if any(sameVec_u0u1u2) or -1 in lagVec_u0u1u2:
        return [-2],[-1],-1,[]

    lagMax=np.max(lagVec_u0u1u2)
    # print("v0 selected: \n"+str(v0_oneUnitCell[-20:,:]))
    # print("v1 selected: \n"+str(v1_oneUnitCell[-20:,:]))
    # print("v2 selected: \n"+str(v2_oneUnitCell[-20:,:]))
    ##############################
    # print("u0 selected: \n"+str(v0_oneUnitCell[-20:,:]@xiVec))
    # print("u1 selected: \n"+str(v1_oneUnitCell[-20:,:]@xiVec))
    # print("u2 selected: \n"+str(v2_oneUnitCell[-20:,:]@xiVec))

    pVec=[]
    statVec=[]

    #u0
    pTmp,statTmp,lengthTmp=ksTestOneColumn(u0_oneUnitCell,lagMax)
    pVec.append(pTmp)
    statVec.append(statTmp)

    # u1
    pTmp,statTmp,lengthTmp=ksTestOneColumn(u1_oneUnitCell,lagMax)
    pVec.append(pTmp)
    statVec.append(statTmp)

    #u2
    pTmp,statTmp,lengthTmp=ksTestOneColumn(u2_oneUnitCell,lagMax)
    pVec.append(pTmp)
    statVec.append(statTmp)

    # for j in range(0,3):
    #
    #     # v0
    #     pTmp,statTmp,lengthTmp=ksTestOneColumn(v0_oneUnitCell[:,j],lagMax)
    #     pVec.append(pTmp)
    #     statVec.append(statTmp)
    #
    #
    #     # v1
    #     pTmp,statTmp,lengthTmp=ksTestOneColumn(v1_oneUnitCell[:,j],lagMax)
    #     pVec.append(pTmp)
    #     statVec.append(statTmp)
    #
    #
    #
    #     # v2
    #     pTmp,statTmp,lengthTmp=ksTestOneColumn(v2_oneUnitCell[:,j],lagMax)
    #     pVec.append(pTmp)
    #     statVec.append(statTmp)
    numDataPoints=lengthTmp
    return pVec,statVec,numDataPoints,lagVec_u0u1u2


def check_eta_H(eta_H_dir):

    eta_H_sortedDataFilesToRead=sort_data_files_by_flushEnd(eta_H_dataDir)
    # print(eta_H_sortedDataFilesToRead)
    startingFileInd,startingVecPosition=parseSummaryU_Dist()
    if startingFileInd<0:
        #we guess that the equilibrium starts at this file
        startingFileInd=len(eta_H_sortedDataFilesToRead)-lastFileNum

    eta_HStartingFileName=eta_H_sortedDataFilesToRead[startingFileInd]

    with open(eta_HStartingFileName,"rb") as fptr:
        eta_H_inArrStart=np.array(pickle.load(fptr))

    eta_H_Arr=eta_H_inArrStart.reshape((-1,6))

    #read the rest of eta_H pkl files
    for pkl_file in eta_H_sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            eta_H_inArr=np.array(pickle.load(fptr))
            eta_H_inArr=eta_H_inArr.reshape((-1,6))
            eta_H_Arr=np.concatenate((eta_H_Arr,eta_H_inArr),axis=0)

    sameVec_eta_H_1_6=[]
    lagVec_eta_H_1_6=[]

    _,nCol=eta_H_Arr.shape

    for j in range(0,nCol):
        sameTmp,lagTmp=auto_corrForOneVec(eta_H_Arr[:,j])
        sameVec_eta_H_1_6.append(sameTmp)
        lagVec_eta_H_1_6.append(lagTmp)
    # print(sameVec)
    # print(lagVec_eta_H_1_6)
    if any(sameVec_eta_H_1_6) or -1 in lagVec_eta_H_1_6:
        return [-2],[-1],-1,[]

    lagMax=np.max(lagVec_eta_H_1_6)
    pVec=[]
    statVec=[]

    for j in range(0,nCol):
        pTmp,statTmp,lengthTmp=ksTestOneColumn(eta_H_Arr[:,j],lagMax)
        pVec.append(pTmp)
        statVec.append(statTmp)

    numDataPoints=lengthTmp
    return pVec,statVec,numDataPoints,lagVec_eta_H_1_6














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
i_val=np.random.choice(discrete_vals, size=1, replace=True)[0]
j_val=np.random.choice(discrete_vals, size=1, replace=True)[0]
k_val=np.random.choice(discrete_vals, size=1, replace=True)[0]

v0_dataDir=U_dist_dataDir+"/v0/"
v1_dataDir=U_dist_dataDir+"/v1/"
v2_dataDir=U_dist_dataDir+"/v2/"
eta_H_dataDir=U_dist_dataDir+"/eta_H/"

pVec_u0u1u2,statVec_u0u1u2,numDataPoints_u0u1u2,lagVec_u0u1u2=check_oneDistDataFilesForOneT(v0_dataDir,v1_dataDir,v2_dataDir,i_val,j_val,k_val)
# print(pVec_v0v1v2)
# print(statVec)
# print(numDataPoints_v0v1v2)
print("lagVec_u0u1u2="+str(lagVec_u0u1u2))
pVec+=pVec_u0u1u2
statVec+=statVec_u0u1u2
pVec_eta_H_1_6,statVec_eta_H_1_6,numDataPoints_eta_H_1_6,lagVec_eta_H_1_6=check_eta_H(eta_H_dataDir)
# print(pVec_eta_H_1_6)
# print(numDataPoints_eta_H_1_6)
pVec+=pVec_eta_H_1_6
statVec+=statVec_eta_H_1_6
lagVecAll=[lagUTmp]+lagVec_u0u1u2+lagVec_eta_H_1_6
lagMax=np.max(lagVecAll)
numDataPoints=np.min([numDataPointsU,numDataPoints_u0u1u2,numDataPoints_eta_H_1_6])
print("lagMax="+str(lagMax))
print("numDataPoints="+str(numDataPoints))

############################################
summary_U_distFile=TDirRoot+"/summary_U_dist.txt"
statThreshhold=0.1
if pVec[0]==-2:
    with open(summary_U_distFile,"w+") as fptr:
        msg="error: same\n"
        fptr.writelines(msg)
        exit(sameErrCode)

if numDataPoints<0:
    msg="high correlation"
    with open(summary_U_distFile,"w+") as fptr:
        fptr.writelines(msg)
    exit(0)


if (np.min(pVec)>=0.01 or np.max(statVec)<=statThreshhold) and numDataPoints>=200:
    if numDataPoints>=effective_data_num_required:
        newDataPointNum=0
    else:
        newDataPointNum=effective_data_num_required-numDataPoints

    msg="equilibrium\n" \
        +"lag="+str(lagMax)+"\n" \
        +"numDataPoints="+str(numDataPoints)+"\n" \
        +"startingFileInd="+str(startingFileInd)+"\n" \
        +"startingVecPosition="+str(startingVecPosition)+"\n" \
        +"newDataPointNum="+str(newDataPointNum)+"\n"
    print(msg)
    with open(summary_U_distFile,"w+") as fptr:
        fptr.writelines(msg)
    exit(0)

#continue
continueMsg="continue\n"

if not (np.min(pVec)>=0.01 or np.max(statVec)<=statThreshhold):
    continueMsg+="stat value: "+str(np.max(statVec))+"\n"
    continueMsg+="p value: "+str(np.min(pVec))+"\n"

if numDataPoints<200:
    #not enough data number

    continueMsg+="numDataPoints="+str(numDataPoints)+" too low\n"
    continueMsg+="lag="+str(lagMax)+"\n"
print(continueMsg)
with open(summary_U_distFile,"w+") as fptr:
    fptr.writelines(continueMsg)
exit(0)
