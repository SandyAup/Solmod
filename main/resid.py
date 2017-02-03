#!/usr/bin/python2.7
# -*-coding:Utf-8 -*

import time
import sys
sys.path.insert(0, '/ams/aupetit/Documents/Pheno/Solmod/Lib/')
from Physics import *
from ExtractData import *
from PlotLib import *


import numpy as np
import matplotlib.pyplot as plt

data_R    = np.loadtxt("results/FF/H_He/Rig_Spline/ISFlux_H_He.txt")
data_Rcs  = np.loadtxt("results/FF/H_He/Rig_CS/ISFlux_H_He.txt")
data_ekn  = np.loadtxt("results/FF/H_He/Ekn_Spline/ISFlux_H_He.txt")

ratio = data_R / data_ekn
ratioCS = data_Rcs / data_R

plt.xscale("log")
plt.yscale("log")
#plt.plot( data_R[:,0], data_R[:,1]*np.power(data_R[:,0], 2.8), "b-", label="R_{H}" )
#plt.plot( data_ekn[:,0], data_ekn[:,1]*np.power(data_ekn[:,0], 2.8), "r-", label="Ekn_{H}" )
#plt.plot( data_R[:,0], 20*data_R[:,2]*np.power(data_R[:,0], 2.8), "b--", label="R_{He}" )
#plt.plot( data_ekn[:,0], 20*data_ekn[:,2]*np.power(data_ekn[:,0], 2.8), "r--", label="Ekn_{He}" )
#plt.plot( data_ekn[:,0], ratio[:,1], label="R_ekn" )
#plt.plot( data_ekn[:,0], ratio[:,2], label="R_ekn" )

#plt.plot( data_Rcs[:,0], 20*data_Rcs[:,1]*np.power(data_Rcs[:,0], 2.8), "b--", label="Rcs_{H}" )
#plt.plot( data_Rcs[:,0], 20*data_Rcs[:,2]*np.power(data_Rcs[:,0], 2.8), "r--", label="Rcs_{He}" )

plt.plot( data_Rcs[:,0], ratioCS[:,1],"b-", label="R-Rcs_{H}" )
plt.plot( data_Rcs[:,0], ratioCS[:,2],"r-", label="R-Rcs_{He}" )



plt.legend()


plt.show()
