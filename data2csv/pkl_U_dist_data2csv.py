import numpy as np
from datetime import datetime
import sys
import re
import glob
import os
import json
from pathlib import Path
import pandas as pd
import pickle
#this script extracts effective data from pkl files

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()


# rowNum=int(sys.argv[1])

N=int(sys.argv[1])
dataRoot="../dataAll/dataAllUnitCell"+str(N)+"/"
obs_U_dist="U_dist"

#search directory
TVals=[]
TFileNames=[]
TStrings=[]
for TFile in glob.glob(dataRoot+"/T*"):
    # print(TFile)
    matchT=re.search(r"T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",TFile)
    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))
        TStrings.append("T"+matchT.group(1))


#sort T values
sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]
sortedTStrings=[TStrings[ind] for ind in sortedInds]


def parseSummary(oneTFolder,obs_name):

    startingFileInd=-1
    startingVecPosition=-1
    lag=-1
    sweep_to_write=-1

    smrFile=oneTFolder+"/summary_"+obs_name+".txt"
    summaryFileExists=os.path.isfile(smrFile)
    if summaryFileExists==False:
        return startingFileInd,startingVecPosition,-1

    with open(smrFile,"r") as fptr:
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

        #match lag
        matchLag=re.search(r"lag=(\d+)",oneLine)
        if matchLag:
            lag=int(matchLag.group(1))
        #match sweep_to_write
        match_sweep_to_write=re.search(r"sweep_to_write=(\d+)",oneLine)
        if match_sweep_to_write:
            sweep_to_write=int(match_sweep_to_write.group(1))
    return startingFileInd, startingVecPosition,lag,sweep_to_write



def sort_data_files_by_flushEnd(oneTFolder,obs_name,varName):
    """

    :param oneTFolder: Txxx
    :param obs_name: data files sorted by flushEnd
    :return:
    """

    dataFolderName=oneTFolder+"/"+obs_name+"_dataFiles/"+varName+"/"
    dataFilesAll=[]
    flushEndAll=[]

    for oneDataFile in glob.glob(dataFolderName+"/flushEnd*.pkl"):
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"flushEnd(\d+)",oneDataFile)
        if matchEnd:
            flushEndAll.append(int(matchEnd.group(1)))


    endInds=np.argsort(flushEndAll)
    # sweepStartSorted=[sweepStartAll[i] for i in startInds]
    sortedDataFiles=[dataFilesAll[i] for i in endInds]

    return sortedDataFiles

def U_extract_ForOneT(oneTFolder,oneTStr,startingFileInd,startingVecPosition,lag):
    TRoot=oneTFolder
    sortedUDataFilesToRead=sort_data_files_by_flushEnd(TRoot,obs_U_dist,"U")
    # print(sortedUDataFilesToRead)
    startingUFileName=sortedUDataFilesToRead[startingFileInd]

    with open(startingUFileName,"rb") as fptr:
        inUStart=pickle.load(fptr)
    # print(len(sortedUDataFilesToRead))
    # print(startingUFileName)
    UVec=inUStart[startingVecPosition:]
    for pkl_file in sortedUDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            # print(pkl_file)
            in_UArr=pickle.load(fptr)
            UVec=np.append(UVec,in_UArr)

    UVecSelected=UVec[::lag]
    # print("len(UVecSelected)="+str(len(UVecSelected)))

    return UVecSelected




def one_v_extract_ForOneT(oneTFolder,oneTStr,startingFileInd,startingVecPosition,lag,component_name,sweep_to_write):
    TRoot=oneTFolder
    sorted_one_v_DataFilesToRead=sort_data_files_by_flushEnd(TRoot,obs_U_dist,component_name)

    one_v_StaringFileName=sorted_one_v_DataFilesToRead[startingFileInd]

    with open(one_v_StaringFileName,"rb") as fptr:
        one_v_inArrStart=np.array(pickle.load(fptr))

    one_v_Arr=one_v_inArrStart.reshape((sweep_to_write,-1,5))
    #read the rest of one_v pkl files
    for pkl_file in sorted_one_v_DataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            one_v_inArr=np.array(pickle.load(fptr))
            one_v_inArr=one_v_inArr.reshape((sweep_to_write,-1,5))
            one_v_Arr=np.concatenate((one_v_Arr,one_v_inArr),axis=0)

    one_v_ArrSelected=one_v_Arr[::lag,:,:]
    return one_v_ArrSelected


def unflatten_index(ind, N):
    i_val = ind // (N * N)  # Get the i index
    j_val = (ind % (N * N)) // N  # Get the j index
    k_val = ind % N  # Get the k index
    return i_val, j_val, k_val

def save_one_unitCell_data_one_v(one_v_ArrSelected,ind,component_name,oneTStr):
    arr_one_slice=one_v_ArrSelected[:,ind,:]
    i_val, j_val, k_val=unflatten_index(ind,N)
    outCsvDataRoot=dataRoot+"/csvOutAll/"
    outCsvFolder=outCsvDataRoot+"/"+oneTStr+"/"
    Path(outCsvFolder).mkdir(exist_ok=True,parents=True)

    outFileName=component_name+"_"+str(i_val)+"_"+str(j_val)+"_"+str(k_val)+".csv"

    outCsvFile=outCsvFolder+outFileName

    df=pd.DataFrame(arr_one_slice)
    # Save to CSV
    df.to_csv(outCsvFile, index=False, header=False)

def save_all_unitCell_data_one_v(one_v_ArrSelected,component_name,oneTStr):
    indLength=one_v_ArrSelected.shape[1]
    for ind in range(0,indLength):
        save_one_unitCell_data_one_v(one_v_ArrSelected,ind,component_name,oneTStr)

def eta_H_extract_ForOneT(oneTFolder,oneTStr,startingFileInd,startingVecPosition,lag):
    TRoot=oneTFolder
    sorted_eta_H_DataFilesToRead=sort_data_files_by_flushEnd(TRoot,obs_U_dist,"eta_H")

    eta_H_StaringFileName=sorted_eta_H_DataFilesToRead[startingFileInd]

    with open(eta_H_StaringFileName,"rb") as fptr:
        eta_H_inArrStart=np.array(pickle.load(fptr))

    eta_H_Arr=eta_H_inArrStart.reshape((-1,6))
    #read the rest of eta_H pkl files

    for pkl_file in sorted_eta_H_DataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            eta_H_inArr=np.array(pickle.load(fptr))
            eta_H_inArr=eta_H_inArr.reshape((-1,6))
            eta_H_Arr=np.concatenate((eta_H_Arr,eta_H_inArr),axis=0)

    eta_H_ArrSelected=eta_H_Arr[::lag,:]

    return eta_H_ArrSelected

def save_eta_H_data(eta_H_ArrSelected,oneTStr):
    outCsvDataRoot=dataRoot+"/csvOutAll/"
    outCsvFolder=outCsvDataRoot+"/"+oneTStr+"/"
    Path(outCsvFolder).mkdir(exist_ok=True,parents=True)
    outFileName="eta+H.csv"
    outCsvFile=outCsvFolder+outFileName
    df=pd.DataFrame(eta_H_ArrSelected)

    # Save to CSV
    df.to_csv(outCsvFile, index=False, header=False)
def save_U_data(UVecSelected,oneTStr):
    outCsvDataRoot=dataRoot+"/csvOutAll/"
    outCsvFolder=outCsvDataRoot+"/"+oneTStr+"/"
    Path(outCsvFolder).mkdir(exist_ok=True,parents=True)

    outFileName="U.csv"
    outCsvFile=outCsvFolder+outFileName
    df=pd.DataFrame(UVecSelected)

    # Save to CSV
    df.to_csv(outCsvFile, index=False, header=False)






for k in range(0,len(sortedTFiles)):
    tStart=datetime.now()
    oneTFolder=sortedTFiles[k]
    oneTStr=sortedTStrings[k]
    startingfileIndTmp,startingVecIndTmp,lagTmp,sweep_to_writeTmp=parseSummary(oneTFolder,obs_U_dist)
    if startingfileIndTmp<0:
        print("summary file does not exist for "+oneTStr+" "+obs_U_dist)
        continue

    UVecSelected=U_extract_ForOneT(oneTFolder,oneTStr,startingfileIndTmp,startingVecIndTmp,lagTmp)
    save_U_data(UVecSelected,oneTStr)

    v0_ArrSelected= one_v_extract_ForOneT(oneTFolder,oneTStr,startingfileIndTmp,startingVecIndTmp,lagTmp,"v0",sweep_to_writeTmp)
    v1_ArrSelected= one_v_extract_ForOneT(oneTFolder,oneTStr,startingfileIndTmp,startingVecIndTmp,lagTmp,"v1",sweep_to_writeTmp)
    v2_ArrSelected= one_v_extract_ForOneT(oneTFolder,oneTStr,startingfileIndTmp,startingVecIndTmp,lagTmp,"v2",sweep_to_writeTmp)

    save_all_unitCell_data_one_v(v0_ArrSelected,"v0",oneTStr)

    save_all_unitCell_data_one_v(v1_ArrSelected,"v1",oneTStr)
    save_all_unitCell_data_one_v(v2_ArrSelected,"v2",oneTStr)

    eta_H_ArrSelected=eta_H_extract_ForOneT(oneTFolder,oneTStr,startingfileIndTmp,startingVecIndTmp,lagTmp)

    save_eta_H_data(eta_H_ArrSelected,oneTStr)
    tEnd=datetime.now()
    print("processed T="+str(sortedTVals[k])+": ",tEnd-tStart)


