import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd
import scipy.stats as stats
from matplotlib.ticker import ScalarFormatter

#This script loads csv data and plot eta_H, with confidence interval
if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()

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


def generate_one_eta_H(oneTFile):
    matchT=re.search(r'T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',oneTFile)
    TVal=float(matchT.group(1))

    eta_path=oneTFile+"/eta_H.csv"
    df=pd.read_csv(eta_path,header=None)

    num_rows,_ = df.shape

    eta_H1Vec=df.iloc[:,0]
    eta_H2Vec=df.iloc[:,1]
    eta_H3Vec=df.iloc[:,2]
    eta_H4Vec=df.iloc[:,3]
    eta_H5Vec=df.iloc[:,4]
    eta_H6Vec=df.iloc[:,5]

    mean_eta_H1=np.mean(eta_H1Vec)
    mean_eta_H2=np.mean(eta_H2Vec)
    mean_eta_H3=np.mean(eta_H3Vec)
    mean_eta_H4=np.mean(eta_H4Vec)
    mean_eta_H5=np.mean(eta_H5Vec)
    mean_eta_H6=np.mean(eta_H6Vec)

    return [mean_eta_H1,mean_eta_H2,mean_eta_H3,mean_eta_H4,mean_eta_H5,mean_eta_H6]



mean_eta_H_all=[]#each row corresponds to one T

for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    one_row=generate_one_eta_H(oneTFile)
    mean_eta_H_all.append(one_row)

sortedTVals=np.array(sortedTVals)
mean_eta_H_all=np.array(mean_eta_H_all)
TInds=np.where(sortedTVals<300)
TToPlt=sortedTVals[TInds]

for j in range(0,6):
    eta_H_one_vec=mean_eta_H_all[:,j]
    eta_H_one_vec_toPlot=eta_H_one_vec[TInds]
    fig,ax=plt.subplots()
    ax.errorbar(TToPlt,eta_H_one_vec_toPlot,fmt='o',color="blue", ecolor='r', capsize=5)
    ax.set_xlabel('$T$')
    ax.set_ylabel(rf"$\eta_{j+1}$")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.set_title(rf"$\eta_{j+1}$")
    plt.savefig(csvDataFolderRoot+"/eta_H"+str(j+1)+".png")
    plt.close()