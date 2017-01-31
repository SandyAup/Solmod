#!/usr/bin/python2.7
# -*-coding:Utf-8 -*

from math import *
import numpy as np
import datetime

from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.optimize import fmin_slsqp
import scipy.interpolate as inter

from datetime import datetime

import time
import sys
sys.path.insert(0, '../Lib/')
from Physics import *
from ExtractData import *
from PlotLib import *

import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, YEARLY, MonthLocator, RRuleLocator, rrulewrapper, DateFormatter
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

np.set_printoptions(precision=4, linewidth=200)

#---------------------------------------------------------------
def main():

	list_CR 	= ['H','He']
	list_exp 	= RefDataset()
	Nexp = len(list_exp)
	directory	= "results/"
	powr = 0

	#------------------------------------------
	# Comparison between J_IS(FF) and J_IS(1D)
	#------------------------------------------

	# Loading ISFlux file

	directory_FF = directory + "FF" + "/" + "_".join(list_CR[i] for i in xrange(len(list_CR))) + "/"
	directory_1D = directory + "1D" + "/" + "_".join(list_CR[i] for i in xrange(len(list_CR))) + "/"
	if (powr is not 0) : directory_1D += "Kdiff/"
	
	filename_IS = "ISFlux_" + "_".join(list_CR[i] for i in xrange(len(list_CR))) 
	if powr is not 0 : stK = str(powr) ; filename_IS += "_K" + stK.replace(".", "_")
	else : filename_IS += ".txt"

	file_IS_FF  = np.loadtxt(directory_FF + filename_IS)
	file_IS_1D  = np.loadtxt(directory_1D + filename_IS)

	dict_IS = {}
	dict_IS["E_IS_FF"] = file_IS_FF[:,0]
	dict_IS["E_IS_1D"] = file_IS_1D[:,0]

	if np.array_equal(dict_IS["E_IS_FF"], dict_IS["E_IS_1D"])  : print "E_IS (1D/FF) consistent - Calculation of the ratio JIS..."
	for i in xrange(len(list_CR)) : dict_IS["ISFlux_FF_%s" % list_CR[i]] = file_IS_FF[:,i+1]
	for i in xrange(len(list_CR)) : dict_IS["ISFlux_1D_%s" % list_CR[i]] = file_IS_1D[:,i+1]

	#print dict_IS["ISFlux_1D_H"]
	#print dict_IS["ISFlux_FF_H"]

	ratio_JIS_H   = dict_IS["ISFlux_1D_H"] / dict_IS["ISFlux_FF_H"]
	ratio_JIS_He = dict_IS["ISFlux_1D_He"] / dict_IS["ISFlux_FF_He"]

	#print ratio_JIS_H

	max_H = max(ratio_JIS_H)
	max_He = max(ratio_JIS_He)
	ratio_JIS_H /= max_H 		# Renormalisation a cause du pb sur J_IS(1D) a haute energie
	ratio_JIS_He /= max_He
	
	# Bandes d'incertitudes
	uncert_IS = np.loadtxt("../Data/Uncertainties_Ghelfi/uncertainties_IS.dat")
	x1 		= uncert_IS[:,0]
	ymed 	= uncert_IS[:,1] 
	yup1 	= uncert_IS[:,2]
	ylow1 	= uncert_IS[:,3]
	
	ymed1 	= ymed /ymed
	yup1 	/= ymed
	ylow1 	/= ymed

	# Plot
	fig_IS_H = SetFig()
	plot_ratio_H = plt.plot(dict_IS["E_IS_1D"], ratio_JIS_H, 'r--', lw = 2.5, ms=3., label=r"H")
	plot_ratio_He = plt.plot(dict_IS["E_IS_1D"], ratio_JIS_He, 'b--', lw = 2.5, ms=3., label=r"He")
	plt.plot(x1, yup1, c='black', lw = 2., label=r"68\% confidence-level from [1]")
	plt.plot(x1, ylow1, c='black', lw = 2.)
	SetLegend(localisation = 4)
 	SetAxis(r"${\rm{}E}_{k/n} \rm{[GeV/n]}$", r"$\rm{J^{IS}_{1D} / J^{IS}_{FF}}$", 0.5, 1.5e3, "None", "None", "xlog")
 	#plt.annotate(r"[H]", fontsize=20, xy=(0.03, 0.94), xycoords='axes fraction')
 	#SaveFig("plots","test", "png")

 	'''
 	fig_IS_He = SetFig()
	plot_ratio = plt.plot(dict_IS["E_IS_1D"], ratio_JIS_He, 'r--', lw = 2.5, ms=3.)
	plt.plot(x1, yup1, c='black', lw = 2., label=r"68\% confidence-level from [1]")
	plt.plot(x1, ylow1, c='black', lw = 2.)
	SetLegend()
 	SetAxis(r"${\rm{}E}_{k/n} \rm{[GeV/n]}$", r"$\rm{J^{IS}_{1D} / J^{IS}_{FF}}$", 0.5, 1.5e3, "None", "None", "xlog")
 	plt.annotate(r"[He]", fontsize=20, xy=(0.03, 0.94), xycoords='axes fraction')
	'''
 	#plt.show()

 	
 	

 	#----------------------------------------------------
	# Reidual between J_TOA / Jdata for FF and 1D models
	#----------------------------------------------------

	filename_FF_TOA = []
	filename_1D_TOA = []
	for i in xrange(len(list_CR)) : 
		filename = "TOAFlux_%s" % list_CR[i]
		filename_FF = filename + ".txt"
		if powr is not 0 : stK = str(powr) ; filename_1D = filename + "_K" + stK.replace(".", "_") + ".txt"
		else : filename_1D = filename + ".txt"
		filename_FF_TOA.append(filename_FF)
		filename_1D_TOA.append(filename_1D)

	file_FF_TOA = []
	file_1D_TOA = []
	for i in range(0, len(list_CR)):
		file_FF_TOA.append(np.loadtxt(directory_FF + filename_FF_TOA[i]))
		file_1D_TOA.append(np.loadtxt(directory_1D + filename_1D_TOA[i]))

	list_exp_red = []

	for i in xrange(len(list_CR)) :
		#print directory + filename_TOA[i]
		search_exp = open(directory_FF + filename_FF_TOA[i])
		for l in search_exp :
			if "Experiment" in l : line = l		
		list_exp_tmp = []
		for name_exp in list_exp :
			if name_exp in line : list_exp_tmp.append(name_exp)
		list_exp_red.append(list_exp_tmp)

	#print list_exp_red 

	dict_FF_TOA = {}
	dict_1D_TOA = {}
	for i in xrange(len(list_CR)) : 
		dict_FF_TOA["E_TOA_%s" % list_CR[i]] = file_FF_TOA[i][:,0]
		dict_1D_TOA["E_TOA_%s" % list_CR[i]] = file_1D_TOA[i][:,0]
		for j in xrange(len(list_exp_red[i])) :
			dict_FF_TOA["TOAFlux_%s_%s" % (list_CR[i], list_exp_red[i][j])] = file_FF_TOA[i][:,j+1]
			dict_1D_TOA["TOAFlux_%s_%s" % (list_CR[i], list_exp_red[i][j])] = file_1D_TOA[i][:,j+1]

	#print dict_1D_TOA["TOAFlux_H_AMS02"]

	


	ratio_list_FF_TOA, ratio_list_1D_TOA = ([] for _ in xrange(2)) 
	
	Nexp, Ndata, dict_Data, list_exp_CR = ExtractDataDict(list_exp, list_CR)


	
	dict_ratio_TOA = {}

	for i in range(0, len(list_CR)):
		for name in list_exp_red[i] :
	
			flux_1D_TOA_interp = np.interp(dict_Data["Edata_%s_%s" % (list_CR[i],name)], dict_1D_TOA["E_TOA_%s" % list_CR[i]], dict_1D_TOA["TOAFlux_%s_%s" % (list_CR[i], name)])
			flux_FF_TOA_interp = np.interp(dict_Data["Edata_%s_%s" % (list_CR[i],name)], dict_FF_TOA["E_TOA_%s" % list_CR[i]], dict_FF_TOA["TOAFlux_%s_%s" % (list_CR[i], name)])

			dict_ratio_TOA["Ratio_1D_%s_%s" % (list_CR[i], name)] = flux_1D_TOA_interp / dict_Data["ydata_%s_%s" % (list_CR[i],name)]
			dict_ratio_TOA["Ratio_FF_%s_%s" % (list_CR[i], name)] = flux_FF_TOA_interp / dict_Data["ydata_%s_%s" % (list_CR[i],name)]

			'''
			ratio_1D_TOA = np.zeros(len(Edata[i-1]))
			ratio_FF_TOA = np.zeros(len(Edata[i-1]))
			for j in range(0,len(Edata[i-1])):
				ratio_1D_TOA[j] = flux_1D_TOA_interp[j] / ydata[i-1][j]
				ratio_FF_TOA[j] = flux_FF_TOA_interp[j] / ydata[i-1][j]
			ratio_list_1D_TOA.append(ratio_1D_TOA)
			ratio_list_FF_TOA.append(ratio_FF_TOA)
			'''

	'''
	uncert_pamela = np.loadtxt(directory + "final_results/uncertainties_PAMELA2006.dat")
	x_pam = uncert_pamela[:,0]
	yup_pam = uncert_pamela[:,1]
	ylow_pam = uncert_pamela[:,2]
	uncert_ams = np.loadtxt(directory + "final_results/uncertainties_AMS02.dat")
	x_ams = uncert_ams[:,0]
	yup_ams = uncert_ams[:,1]
	ylow_ams = uncert_ams[:,2]
	'''

	# Plot
	for j in range(0, len(list_CR)):
		f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
		f.set_facecolor('white')	
		for name in list_exp_red[j] :
			if name == "AMS02" or name == "PAMELA2006" :
				if name == "AMS02" :
					ax1.plot(dict_Data["Edata_%s_%s" % (list_CR[j],name)], dict_ratio_TOA["Ratio_1D_%s_%s" % (list_CR[j], name)], 'bo', ms=5., label = "1D model")
					ax1.plot(dict_Data["Edata_%s_%s" % (list_CR[j],name)], dict_ratio_TOA["Ratio_FF_%s_%s" % (list_CR[j], name)], 'rD', ms=5., label= "Force field")
				#ax1.plot(Edata[i], ratio_list_FF_TOA[i], 'rD', ms=5., label= "Force field")
				elif name == "PAMELA2006" :
					ax2.plot(dict_Data["Edata_%s_%s" % (list_CR[j],name)], dict_ratio_TOA["Ratio_1D_%s_%s" % (list_CR[j], name)], 'bo', ms=5., label = "1D model")
					ax2.plot(dict_Data["Edata_%s_%s" % (list_CR[j],name)], dict_ratio_TOA["Ratio_FF_%s_%s" % (list_CR[j], name)], 'rD', ms=5., label= "Force field")

		SetLegendSub(ax1, 2)
		SetAxisSub(ax1, ax2, r"$\rm{E_{k/n} [GeV/n]}$", r"$\rm{J^{TOA}_{best-fit} / J_{data}}$", 5e-2, 2e3, 0.85, 1.2, "sharex", "xlog")
		plt.locator_params(axis='y',nbins=4)
	
	plt.show()
		
	'''
	if list_exp[i] == "PAMELA2006" :
		ax2.plot(Edata[i], ratio_list_1D_TOA[i], 'bo', ms=5., label = "1D model")
		ax2.plot(Edata[i], ratio_list_FF_TOA[i], 'rD', ms=5., label= "Force field")
	# Fine-tune figure; make subplots close to each other and hide x ticks for
	# all but bottom plot.
	if list_exp[i] == "PAMELA2006":
		ax2.plot(x_pam, yup_pam, c='black', lw = 2.5, label="Data uncertainties")
		ax2.plot(x_pam, ylow_pam, c='black', lw = 2.5)
	if list_exp[i] == "AMS02":
		ax1.plot(x_ams, yup_ams, c='black', lw = 2.5, label="Data uncertainties")
		ax1.plot(x_ams, ylow_ams, c='black', lw = 2.5)
	'''

	
#---------------------------------------------------------------	

main()
