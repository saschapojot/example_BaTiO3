import pickle
import numpy as np
from datetime import datetime
import pandas as pd
import statsmodels.api as sm
import sys
import re
import warnings
import subprocess


from scipy.stats import ks_2samp
import glob
from pathlib import Path
import os
import json
import pickle
import matplotlib.pyplot as plt



#this file plots U
argErrCode=2
confErrCode=4
# if (len(sys.argv)!=2):
#     print("wrong number of arguments")
#     exit(argErrCode)

# confFileName=str(sys.argv[1])
#
# #parse conf, get jsonDataFromConf
# confResult=subprocess.run(["python3", "./init_run_scripts/parseConf.py", confFileName], capture_output=True, text=True)
# confJsonStr2stdout=confResult.stdout
# print(confJsonStr2stdout)
# if confResult.returncode !=0:
#     print("Error running parseConf.py with code "+str(confResult.returncode))
#     # print(confResult.stderr)
#     exit(confErrCode)
# match_confJson=re.match(r"jsonDataFromConf=(.+)$",confJsonStr2stdout)
# if match_confJson:
#     jsonDataFromConf=json.loads(match_confJson.group(1))
# else:
#     print("jsonDataFromConf missing.")
#     exit(confErrCode)
# # print(jsonDataFromConf)

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
TStr=str(290)

N=2

dataRoot="./dataAll/dataAllUnitCell"+str(N)+"/T"+TStr+"/U_dist_dataFiles/"

UData_dir=dataRoot+"/U/"

U_sortedDataFilesToRead=sort_data_files_by_flushEnd(UData_dir)

startingFileInd=3
startingFileName=U_sortedDataFilesToRead[startingFileInd]
# print("U, startingFileName="+str(startingFileName))
with open(startingFileName,"rb") as fptr:
    arr=pickle.load(fptr)

for pkl_file in U_sortedDataFilesToRead[(startingFileInd+1):]:
    with open(pkl_file,"rb") as fptr:
        inArr=pickle.load(fptr)
        # print("len(inArr)="+str(len(inArr)))
    arr=np.append(arr,inArr)

plt.figure()
plt.plot(arr)
plt.savefig("tmp.png")