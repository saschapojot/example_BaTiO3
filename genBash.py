from pathlib import Path
from decimal import Decimal, getcontext

import numpy as np
import pandas as pd


#this script creates slurm bash files

def format_using_decimal(value, precision=10):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

outPath="./bashFiles/"
Path(outPath).mkdir(exist_ok=True,parents=True)

TVals=[110,120,130,140,150,160,170,180]
TStrAll=[]
# print(TDirsAll)
for k in range(0,len(TVals)):
    T=TVals[k]
    # print(T)

    TStr=str(T)#format_using_decimal(T)
    TStrAll.append(TStr)

def contents_to_bash(k):
    TStr=TStrAll[k]

    contents=[
    "#!/bin/bash\n",
        "#SBATCH -n 24\n",
        "#SBATCH -N 1\n",
        "#SBATCH -t 0-60:00\n",
        "#SBATCH -p hebhcnormal01\n"
        "#SBATCH --mem=40GB\n",
        f"#SBATCH -o outmcT{TStr}.out\n",
        f"#SBATCH -e outmcT{TStr}.out\n",
        "cd /public/home/hkust_jwliu_1/liuxi/Document/cppCode/example_BaTiO3\n",
        f"python3 -u exec_checking.py {TStr} 2"
        ]

    outBashName=outPath+f"/run_mcT{TStr}.sh"
    with open(outBashName,"w+") as fptr:
        fptr.writelines(contents)

# for k in range(0,len(TStrAll)):
#     contents_to_bash(k)

N=2
def check_after1run_bash(k):
    TStr=TStrAll[k]
    contents=[
        "#!/bin/bash\n",
        "#SBATCH -n 24\n",
        "#SBATCH -N 1\n",
        "#SBATCH -t 0-60:00\n",
        "#SBATCH -p hebhcnormal01\n"
        "#SBATCH --mem=40GB\n",
        f"#SBATCH -o outmcT{TStr}.out\n",
        f"#SBATCH -e outmcT{TStr}.out\n",
        "cd /public/home/hkust_jwliu_1/liuxi/Document/cppCode/example_BaTiO3\n",
        f"python3 -u exec_noChecking.py  {TStr} {N}"
    ]

    outBashName=outPath+f"/run_T{TStr}.sh"
    with open(outBashName,"w+") as fptr:
        fptr.writelines(contents)

for k in range(0,len(TStrAll)):
    check_after1run_bash(k)