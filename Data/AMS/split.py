import os
import glob
import sys
import numpy as np
from math import *

f = np.loadtxt("AMS02_monthly.USINE", dtype="str")

j = 1
for i in range(0,1246) :
	data = open(("data_AMS%i.dat" % j), 'a')
	if os.path.getsize(("data_AMS%i.dat" % j)) == 0:
		data.write("# H ")
		data.write("%s\n" % f[i][1])
		data.write("# <E>  Elo  Eup   y   ystat_lo  ystat_up  ysyst_lo  ysyst_up  yerrtot_lo  yerrtot_up\n")
	if f[i][1] == f[i+1][1] :
		for k in range (3,11):
			data.write(f[i][k])
			data.write(" ")
		yerr_tot_lo = sqrt(pow(float(f[i][7]),2) + pow(float(f[i][9]),2))
		yerr_tot_up = sqrt(pow(float(f[i][8]),2) + pow(float(f[i][10]),2))
		data.write("%.6e %.6e" % (yerr_tot_lo, yerr_tot_up))
		data.write("\n")
	else :
		for k in range (3,11):
			data.write(f[i][k])
			data.write(" ")
		yerr_tot_lo = sqrt(pow(float(f[i][7]),2) + pow(float(f[i][9]),2))
		yerr_tot_up = sqrt(pow(float(f[i][8]),2) + pow(float(f[i][10]),2))
		data.write("%.6e %.6e" % (yerr_tot_lo, yerr_tot_up))
		data.write("\n")
		j += 1
