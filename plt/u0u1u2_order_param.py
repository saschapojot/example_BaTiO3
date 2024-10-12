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
from decimal import Decimal, getcontext

# this script concatenates configurations of u0, u1 ,u2

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()

xi_Ba=0.2
xi_Ti=0.76
xi_O_parallel=-0.53
xi_O_perpendicular=-0.21


xiVec=np.array([xi_Ba,xi_Ti,xi_O_parallel,xi_O_perpendicular,xi_O_perpendicular])

N=int(sys.argv[1])
csvDataFolderRoot="../dataAll/dataAllUnitCell"+str(N)+"/csvOutAll/"
TVals=[]
TFileNames=[]


for TFile in glob.glob(csvDataFolderRoot+"/T*"):

    matchT=re.search(r"T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",TFile)
    # if float(matchT.group(1))<1:
    #     continue

    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))


sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]

def format_using_decimal(value, precision=10):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

def load_one_unitCell(oneTFile,component_name,i,j,k):
    """

    :param component_name: v0,v1,v2
    :param i:
    :param j:
    :param k:
    :return:
    """
    one_v_file=oneTFile+"/"+component_name+"_"+str(i)+"_"+str(j)+"_"+str(k)+".csv"

    csv_arr=np.array(pd.read_csv(one_v_file,header=None))

    u_oneUnitCell=csv_arr@xiVec

    # print(np.mean(u_oneUnitCell))
    return np.mean(u_oneUnitCell)

def v0_to_u0_and_combine(oneTFile):
    component_name="v0"

    abs_vals=[]
    for i in range(0,2):
        for j in range(0,2):
            for k in range(0,2):
                val=load_one_unitCell(oneTFile,component_name,i,j,k)
                abs_vals.append(np.abs(val))

    return np.mean(abs_vals)



def v1_to_u1_and_combine(oneTFile):
    component_name="v1"

    abs_vals=[]
    for i in range(0,2):
        for j in range(0,2):
            for k in range(0,2):
                val=load_one_unitCell(oneTFile,component_name,i,j,k)
                abs_vals.append(np.abs(val))

    return np.mean(abs_vals)

def v2_to_u2_and_combine(oneTFile):
    component_name="v2"

    abs_vals=[]
    for i in range(0,2):
        for j in range(0,2):
            for k in range(0,2):
                val=load_one_unitCell(oneTFile,component_name,i,j,k)
                abs_vals.append(np.abs(val))

    return np.mean(abs_vals)

u0u1u2_vals=[]
for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    oneTVal=sortedTVals[k]
    oneTStr=format_using_decimal(oneTVal)
    u0Tmp=v0_to_u0_and_combine(oneTFile)
    u1Tmp=v1_to_u1_and_combine(oneTFile)
    u2Tmp=v2_to_u2_and_combine(oneTFile)

    oneRow=[u0Tmp,u1Tmp,u2Tmp]
    u0u1u2_vals.append(oneRow)


u0u1u2_vals=np.array(u0u1u2_vals)
colNames = ['u0', 'u1', 'u2']

rowNames=[format_using_decimal(elem) for elem in sortedTVals]



df = pd.DataFrame(u0u1u2_vals, columns=colNames, index=rowNames)


outCsvName=csvDataFolderRoot+"/order_params.csv"

df.to_csv(outCsvName, index=True)  # index=True is the default





# u0u1u2ValsArr=[]
# for k in range(0,len(sortedTFiles)):
#     oneTFile=sortedTFiles[k]
#     oneTVal=sortedTVals[k]
#     oneTStr=format_using_decimal(oneTVal)
#     mean_abs_u0_unitAvg=v0_to_u0_and_combine(oneTFile)
#     mean_abs_u1_unitAvg=v1_to_u1_and_combine(oneTFile)
#     mean_abs_u2_unitAvg=v2_to_u2_and_combine(oneTFile)
#
#     oneRow=[mean_abs_u0_unitAvg,mean_abs_u1_unitAvg,mean_abs_u2_unitAvg]
#     u0u1u2ValsArr.append(oneRow)
#
# print(u0u1u2ValsArr)



# k=0
# oneTFile=sortedTFiles[k]
# oneTVal=sortedTVals[k]
# oneTStr=format_using_decimal(oneTVal)
# print("T="+str(oneTStr))
# print(v0_to_u0_and_combine(oneTFile))
# load_one_unitCell(oneTFile,"v0",0,1,1)
# sumTmp=0
# for i in range(0,2):
#     for j in range(0,2):
#         for k in range(0,2):
#             sumTmp+=load_one_unitCell(oneTFile,"v0",i,j,k)
#
#
# print("sumTmp="+str(sumTmp))