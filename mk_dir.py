from pathlib import Path
from decimal import Decimal, getcontext

import numpy as np
import pandas as pd


#This script creates directories and conf files for mc

######################################################
#Table II values

#on-site
kappa2_val= 0.0568
alpha_val= 0.320
gamma_val= -0.473

# intersite
j1_val= -0.02734
j3_val= 0.00927
j6_val= 0.00370

j2_val= 0.04020
j4_val= -0.00815
j7_val= 0.00185
j5_val= 0.00580

# elastic
B11_val= 4.64
B12_val= 1.65
B44_val= 1.85

#coupling
B1xx_val= -2.18
B1yy_val= -0.20
B4yz_val= -0.08

B100=4.64
B111=-0.2
B412=-0.08

# dipole

ZStar_val= 9.956
epsilon_infty= 5.24

#xi values

xi_Ba=0.2
xi_Ti=0.76
xi_O_parallel=-0.53
xi_O_perpendicular=-0.21

vals2Potential=[kappa2_val,alpha_val,gamma_val,
                j1_val,j2_val,
                j3_val,j4_val,j5_val,
                j6_val,j7_val,
                B11_val,B12_val,B44_val,
                B1xx_val,B1yy_val,B4yz_val,
                ZStar_val,epsilon_infty,
                xi_Ba,xi_Ti,xi_O_parallel,xi_O_perpendicular]


######################################################



def format_using_decimal(value, precision=10):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

formated_vals2Potential=[format_using_decimal(val) for val in vals2Potential]
str_formated_vals2Potential=','.join(map(str, formated_vals2Potential))


TVals=[290,300]
N=2
dataRoot="./dataAll/"
dataOutDir=dataRoot+"/dataAllUnitCell"+str(N)+"/"

TDirsAll=[]
TStrAll=[]
# print(TDirsAll)
for k in range(0,len(TVals)):
    T=TVals[k]
    # print(T)

    TStr=str(T)#format_using_decimal(T)
    TStrAll.append(TStr)
    TDir=dataOutDir+"/T"+TStr+"/"
    TDirsAll.append(TDir)
    Path(TDir).mkdir(exist_ok=True,parents=True)



def contents_to_conf(k):
    """

    :param k: index of T
    :return:
    """

    contents=[
        "#This is the configuration file for mc computations\n",
        "\n"
        "potential_function_name=V_BaTiO3\n",
        "\n" ,
        "#parameters of coefficients\n",

        "#the following row is the values in Table II\n"
        "coefs=["+str_formated_vals2Potential+"]\n",
        "\n",
        "#Temperature\n",
        "T="+TStrAll[k]+"\n",
        "\n",
        "#unit cell number along each direction\n",
        "N="+str(N)+"\n",
        "\n",
        "erase_data_if_exist=False\n",
        "\n",
        "search_and_read_summary_file=True\n"
        "\n",
        "#For the observable name, only digits 0-9, letters a-zA-Z, underscore _ are allowed\n",
        "\n",
        "observable_name=U_dist\n",
        "\n",
        "effective_data_num_required=1000\n",
        "\n",
        "sweep_to_write=100\n",
        "\n",
        "#within each flush,  sweep_to_write mc computations are executed\n",
        "\n",
        "default_flush_num=15\n",
        "\n",
        "h=5e-2\n",
        "\n",
        "sweep_multiple=100\n",
        "\n",
        "lambda="+str(1/N*0.2)+"\n",
        "\n",
        "B100="+format_using_decimal(B100)+"\n",
        "\n",
        "B111="+format_using_decimal(B111)+"\n",
        "\n",
        "B412="+format_using_decimal(B412)+"\n"



    ]

    outConfName=TDirsAll[k]+"/run_T"+TStrAll[k]+".mc.conf"
    with open(outConfName,"w+") as fptr:
        fptr.writelines(contents)



for k in range(0,len(TDirsAll)):
    contents_to_conf(k)