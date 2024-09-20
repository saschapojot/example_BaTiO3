import sys
import glob
import re
import json
from decimal import Decimal, getcontext
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path
import pickle
#this script loads previous data
numArgErr=4
valErr=5
if (len(sys.argv)!=3):
    print("wrong number of arguments.")
    exit(numArgErr)


jsonDataFromConf =json.loads(sys.argv[1])
jsonFromSummary=json.loads(sys.argv[2])

potential_function_name=jsonDataFromConf["potential_function_name"]
U_dist_dataDir=jsonFromSummary["U_dist_dataDir"]
startingFileInd=jsonFromSummary["startingFileInd"]
startingVecPosition=jsonFromSummary["startingVecPosition"]
N=int(jsonDataFromConf["N"])

if N<=0:
    print("N="+str(N)+"<=0")
    exit(valErr)
#search and read U_dist files

#give arbitrary values to v0,v1,v2 if they don't exist

#search flushEnd
pklFileList=[]
flushEndAll=[]
#assume that the v0, v1, v2, eta_H files are intact, we only check v0 directory
for file in glob.glob(U_dist_dataDir+"/v0/flushEnd*.pkl"):
    pklFileList.append(file)
    matchEnd=re.search(r"flushEnd(\d+)",file)
    if matchEnd:
        flushEndAll.append(int(matchEnd.group(1)))
# print(U_dist_dataDir)
flushLastFile=-1

def format_using_decimal(value, precision=10):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)
def create_v_init_alpha_direction(U_dist_dataDir,alpha):
    """
    create initial values for coordinates along direction alpha
    :return:
    """
    #
    random_perturbations=np.random.uniform(-0.1,0.1,5*N**3)

    outPath=U_dist_dataDir+"/v"+format_using_decimal(alpha)+"/"

    Path(outPath).mkdir(exist_ok=True,parents=True)

    outFileName=outPath+"/v"+format_using_decimal(alpha)+"_init.pkl"

    with open(outFileName,"wb") as fptr:
        pickle.dump(random_perturbations,fptr)

def create_eta_H(U_dist_dataDir):
    """
    create initial values for eta_H
    :param U_dist_dataDir:
    :return:
    """
    random_perturbations=np.random.uniform(-0.1,0.1,6)
    outPath=U_dist_dataDir+"/eta_H/"
    Path(outPath).mkdir(exist_ok=True,parents=True)
    outFileName=outPath+"/eta_H_init.pkl"
    with open(outFileName,"wb") as fptr:
        pickle.dump(random_perturbations,fptr)




def create_loadedJsonData(flushLastFileVal):

    initDataDict={

        "flushLastFile":str(flushLastFileVal)
    }
    # print(initDataDict)
    return json.dumps(initDataDict)

#if no data found, return flush=-1
if len(pklFileList)==0:
    create_v_init_alpha_direction(U_dist_dataDir,0)
    create_v_init_alpha_direction(U_dist_dataDir,1)
    create_v_init_alpha_direction(U_dist_dataDir,2)
    create_eta_H(U_dist_dataDir)
    loadedJsonDataStr=create_loadedJsonData(flushLastFile)
    loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
    print(loadedJsonData_stdout)
    exit(0)


#if found pkl data with flushEndxxxx
sortedEndInds=np.argsort(flushEndAll)
sortedflushEnd=[flushEndAll[ind] for ind in sortedEndInds]
loadedJsonDataStr=create_loadedJsonData(sortedflushEnd[-1])
loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
print(loadedJsonData_stdout)
exit(0)
