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
sys.path.insert(0, '/ams/aupetit/Documents/Pheno/Solmod/Lib/')
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

	directory	= "/ams/aupetit/Documents/Pheno/Solmod/main/results/"

	list_CR 	= ['H','He']
	#list_exp 	= ['AMS02', 'BESS00', 'PAMELA2008']
	list_exp 	= RefDataset()
	MODE  		= "FF"
	powr  		= 0

	Nexp, Ndata, Edata_exp, flux_data_exp, sigma_exp, Edata = ExtractData(list_exp, list_CR, "light")

	if (MODE == "FF") :	directory += MODE + "/" + "_".join(list_CR[i] for i in xrange(len(list_CR))) + "/"
	elif (MODE == "1D") :
		directory += MODE + "/" + "_".join(list_CR[i] for i in xrange(len(list_CR))) + "/"
		if (powr is not 0) : directory += "Kdiff/"
	
	filename_IS = "ISFlux_" + "_".join(list_CR[i] for i in xrange(len(list_CR))) + ".txt"

	filename_TOA = []
	for i in xrange(len(list_CR)) : 
		filename = "TOAFlux_%s" % list_CR[i]
		if (MODE == "FF") : filename += ".txt"
		elif (MODE == "1D") : 
			if powr is not 0 : stK = str(powr) ; filename += "_K" + stK.replace(".", "_")
			filename += ".txt"
		filename_TOA.append(filename)

	#print filename_TOA
	 

	#------------------------------------------
	# Loading ISFlux file
	#------------------------------------------

	file_IS  = np.loadtxt(directory + filename_IS)

	E_IS   = file_IS[:,0]
	ISFlux = []
	for i in xrange(len(list_CR)) : ISFlux.append(file_IS[:,i+1])
	

	### TEST dictionnaires en  python

	dict_IS = {}
	dict_IS["E_IS"] = file_IS[:,0]
	for i in xrange(len(list_CR)) : dict_IS["ISFlux_%s" % list_CR[i]] = file_IS[:,i+1]

	#print dict_IS
	
	#------------------------------------------
	# Loading TOAFlux files for each CR species
	#------------------------------------------

	file_TOA = []
	for i in xrange(len(list_CR)) : file_TOA.append(np.loadtxt(directory + filename_TOA[i]))

	list_exp_red = []

	for i in xrange(len(list_CR)) :
		print directory + filename_TOA[i]
		search_exp = open(directory + filename_TOA[i])
		for l in search_exp :
			if "Experiment" in l : line = l		
		list_exp_tmp = []
		for name_exp in list_exp :
			if name_exp in line : list_exp_tmp.append(name_exp)
		list_exp_red.append(list_exp_tmp)

	dict_TOA = {}
	for i in xrange(len(list_CR)) : 
		dict_TOA["E_TOA_%s" % list_CR[i]] = file_TOA[i][:,0]
		for j in xrange(len(list_exp_red[i])) :
			dict_TOA["TOAFlux_%s_%s" % (list_CR[i], list_exp_red[i][j])] = file_TOA[i][:,j+1]

	print dict_TOA["E_TOA_H"]

	#sys.exit()

	#------------------------------------------
	# One plot for each species
	#------------------------------------------

	for i in range(0, len(list_CR)):

		fig_all, ax =plt.subplots(1,1,figsize=(12,8))
		fig_all.set_facecolor('white')

		# IS FLUX
		#----------
		plot_IS = plt.plot(dict_IS["E_IS"], dict_IS["ISFlux_%s" % list_CR[i]]*np.power(dict_IS["E_IS"], 2.), c='black', lw = 2.5, label='Interstellar flux')

		#plot_IS = plt.plot(dict_IS["E_IS"], flux_IS[i]*np.power(E_IS, 2.), c='black', lw = 2.5, label='Interstellar flux')


		# DATA + TOA FLUXES
		#--------------------
		for j in xrange(len(list_exp_red[i])) :
			print list_CR[i],list_exp_red[i][j]
			plot_TOA = plt.plot(dict_TOA["E_TOA_%s" % list_CR[i]], dict_TOA["TOAFlux_%s_%s" % (list_CR[i], list_exp_red[i][j])]*np.power(dict_TOA["E_TOA_%s" % list_CR[i]], 2.), ls ='-',ms=3.) # color = col

		'''
		for j in range(0,Nexp):  
			if list_CR[i] in list_exp_CR[j] :
				index = list_exp_CR[j].index(list_CR[i])
				plot_exp = plt.errorbar(Edata_exp[j][index], flux_data_exp[j][index]*np.power(Edata_exp[j][index], 2.), xerr = None, yerr=sigma_exp[j][index]*np.power(Edata_exp[j][index], 2), fmt='o',ms=3., label=list_exp[j] + (', $\phi$ = %.2f' % Best_phi[j]) + ' GV')
				col = plot_exp[0].get_color()
				#plot_TOA = plt.plot(E_TOA, flux_TOA[j][index]*np.power(E_TOA, 2.), color = col, ls ='-',ms=3.)
		'''

		SetAxis(r"$\rm{E_{k/n}[GeV/n]}$", r"$\rm{J \times E^2_{k/n} [GeV/n.m^{-2}.s^{-1}.sr^{-1}]}$", min(Edata), 1e4, 0.1, 5e5, "xylog")
		SetLegend(ncolumn = 2)
		
		if (MODE == "FF"):
			plt.annotate(r"[FF - "+ list_CR[i] + r"]", fontsize=20, xy=(0.03, 0.94), xycoords='axes fraction')
		elif (MODE == "1D"):
			plt.annotate(r"[1D - "+ list_CR[i] + r"]", fontsize=20, xy=(0.03, 0.94), xycoords='axes fraction')
			plt.annotate(r"$K(r) = K_{0}\times r^{%.2f}$" % powr, fontsize=20, xy=(0.03, 0.88), xycoords='axes fraction')

	plt.show()
		

#---------------------------------------------------------------	

main()
