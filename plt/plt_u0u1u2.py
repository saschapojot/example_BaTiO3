import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd
import scipy.stats as stats
from scipy.stats import gaussian_kde
from pathlib import Path
from decimal import Decimal, getcontext
#This script loads csv data and plot u0, u1, u2, with confidence interval


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
def v0_to_u0(oneTFile):
    matchT=re.search(r'T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',oneTFile)
    TVal=float(matchT.group(1))

    v0_folder=oneTFile+"/v0/"

    v0_filesAll=[]
    for one_v0File in glob.glob(v0_folder+"/v0*.csv"):
        v0_filesAll.append(one_v0File)

    csv_arr=np.array(pd.read_csv(v0_filesAll[0],header=None))

    for csv_file in v0_filesAll[1:]:
        in_csv_arr=pd.read_csv(csv_file,header=None)
        csv_arr=np.concatenate((csv_arr,in_csv_arr),axis=0)

    u0Vec=csv_arr@xiVec

    return u0Vec


def v1_to_u1(oneTFile):
    matchT=re.search(r'T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',oneTFile)
    TVal=float(matchT.group(1))

    v1_folder=oneTFile+"/v1/"

    v1_filesAll=[]
    for one_v1File in glob.glob(v1_folder+"/v1*.csv"):
        v1_filesAll.append(one_v1File)

    csv_arr=np.array(pd.read_csv(v1_filesAll[0],header=None))

    for csv_file in v1_filesAll[1:]:
        in_csv_arr=pd.read_csv(csv_file,header=None)
        csv_arr=np.concatenate((csv_arr,in_csv_arr),axis=0)

    u1Vec=csv_arr@xiVec

    return u1Vec

def v2_to_u2(oneTFile):
    matchT=re.search(r'T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',oneTFile)
    TVal=float(matchT.group(1))

    v2_folder=oneTFile+"/v2/"

    v2_filesAll=[]
    for one_v2File in glob.glob(v2_folder+"/v2*.csv"):
        v2_filesAll.append(one_v2File)

    csv_arr=np.array(pd.read_csv(v2_filesAll[0],header=None))

    for csv_file in v2_filesAll[1:]:
        in_csv_arr=pd.read_csv(csv_file,header=None)
        csv_arr=np.concatenate((csv_arr,in_csv_arr),axis=0)

    u2Vec=csv_arr@xiVec

    return u2Vec

u0Folder=csvDataFolderRoot+"/u0/"
u1Folder=csvDataFolderRoot+"/u1/"

u2Folder=csvDataFolderRoot+"/u2/"

Path(u0Folder).mkdir(parents=True,exist_ok=True)
Path(u1Folder).mkdir(parents=True,exist_ok=True)
Path(u2Folder).mkdir(parents=True,exist_ok=True)
# plot u0
for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    oneTVal=sortedTVals[k]
    oneTStr=format_using_decimal(oneTVal)
    fig, ax = plt.subplots()
    u0Vec=v0_to_u0(oneTFile)
    mean_u0=np.mean(u0Vec)

    ax.hist(u0Vec, bins=30, density=True, alpha=0.6, color='skyblue', edgecolor='black')
    density = gaussian_kde(u0Vec)
    x_vals = np.linspace(min(u0Vec), max(u0Vec), 1000)
    den=density(x_vals)
    ax.plot(x_vals, den, color='red')
    ax.set_xlabel('$u_{0}$')
    ax.set_ylabel('Density')
    ax.set_title(rf'$u_{0}$ histogram with Density Curve, T={oneTStr}K')
    # half_x=(min(u0Vec)+ max(u0Vec))/2
    # half_y=(min(den)+max(den))/2
    # ax.text(half_x, half_y, 'mean u0='+str(np.round(mean_u0,2)), fontsize=12, color='red')
    plt.savefig(u0Folder+"/"+oneTStr+"_u0.png")
    plt.close()


# plot u1
for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    oneTVal=sortedTVals[k]
    oneTStr=format_using_decimal(oneTVal)
    fig, ax = plt.subplots()
    u1Vec=v1_to_u1(oneTFile)
    mean_u1=np.mean(u1Vec)

    ax.hist(u1Vec, bins=60, density=True, alpha=0.6, color='skyblue', edgecolor='black')
    density = gaussian_kde(u1Vec)
    x_vals = np.linspace(min(u1Vec), max(u1Vec), 1000)
    den=density(x_vals)
    ax.plot(x_vals, den, color='red')
    ax.set_xlabel('$u_{1}$')
    ax.set_ylabel('Density')
    ax.set_title(rf'$u_{1}$ histogram with Density Curve, T={oneTStr}K')
    # half_x=(min(u1Vec)+ max(u1Vec))/2
    # half_y=(min(den)+max(den))/2
    # ax.text(half_x, half_y, 'mean u1='+str(np.round(mean_u1,2)), fontsize=12, color='red')
    plt.savefig(u1Folder+"/"+oneTStr+"_u1.png")
    plt.close()



# plot u2
for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    oneTVal=sortedTVals[k]
    oneTStr=format_using_decimal(oneTVal)
    fig, ax = plt.subplots()
    u2Vec=v2_to_u2(oneTFile)
    mean_u2=np.mean(u2Vec)

    ax.hist(u2Vec, bins=60, density=True, alpha=0.6, color='skyblue', edgecolor='black')
    density = gaussian_kde(u2Vec)
    x_vals = np.linspace(min(u2Vec), max(u2Vec), 1000)
    den=density(x_vals)
    ax.plot(x_vals, den, color='red')
    ax.set_xlabel('$u_{2}$')
    ax.set_ylabel('Density')
    ax.set_title(rf'$u_{2}$ histogram with Density Curve, T={oneTStr}K')
    # half_x=(min(u2Vec)+ max(u2Vec))/2
    # half_y=(min(den)+max(den))/2
    # ax.text(half_x, half_y, 'mean u2='+str(np.round(mean_u2,2)), fontsize=12, color='red')
    plt.savefig(u2Folder+"/"+oneTStr+"_u2.png")
    plt.close()


