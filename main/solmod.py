#!/usr/bin/python2.7
# -*-coding:Utf-8 -*

from math import *
import numpy as np

from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.optimize import fmin_slsqp
import scipy.interpolate as inter

import time
import sys
sys.path.insert(0, '/ams/aupetit/Documents/Pheno/Solmod/Lib/')
from Physics import *
from ExtractData import *
from PlotLib import *

import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, YEARLY, RRuleLocator, rrulewrapper, DateFormatter
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

np.set_printoptions(precision=4, linewidth=200)

#-------------------------------------------------------------------------------#
#						  	    MACROS SECONDAIRES								#
#-------------------------------------------------------------------------------#

def InitParameters1D(MODE, Nexp, list_CR) :   

	init_pars = []
	for CR in list_CR :
		if CR == 'H':
			init_pars_H = [3.5, 1.9, -0.5, -1.3, -3., -3.8]
			init_pars += init_pars_H
		elif CR == 'He' :
			init_pars_He = [1.8,1.,-1.,-2.3,-3.8,-4.4]  
			init_pars += init_pars_He
		elif CR == 'Electron' :
			init_pars_Elec = [1.5,0.8,-0.2,-1.,-2.,-3.0]
			init_pars += init_pars_Elec

	N_IS = len(init_pars)

	if (MODE == "FF"):
		for i in range (0, Nexp) :
			phi = 0.8                     
			init_pars.append(phi) 
		#par_IS  = init_pars[0:N_IS]
		#par_phi = init_pars[N_IS:N_IS+Nexp]
		return init_pars, N_IS #par_IS, par_phi

	elif (MODE == "1D"):
		for i in range (0, Nexp) :
			K0 = 3.4      											# log10(K0) default value [AU2/yr]
			init_pars.append(K0) 
		#par_IS = init_pars[0:N_IS]
		#par_K0 = init_pars[N_IS:N_IS+Nexp]
		return init_pars, N_IS #par_IS, par_K0,

#--------------------------------

def IS_Spline_flux(ekn, parIS, name_exp_CR) :

	a_cr, z_cr, mGeV = CREntry(name_exp_CR)
	ekn_log = np.log10(ekn)
	xk = np.array([1., 7., 50., 100., 400., 800.])

	for i in range(0, xk.size):
		xk[i] = R_to_Ekn(xk[i], a_cr, mGeV, z_cr)
	xk_log = np.log10(xk)
	xk_log_min = min(xk_log)
	xk_log_max = max(xk_log)

	logA = 0. ; B = 0.
	ekn_tmp = 0. ; y_tmp = 0.
	result = 0.

	if ekn_log > xk_log_min and ekn_log < xk_log_max :
		result = inter.spline(xk_log, parIS, ekn_log, order=3, kind='smoothest')

	# Extrapolation loi de puissance : y = A x^B
	elif ekn_log < xk_log_min :
		ekn_tmp = np.log10(R_to_Ekn(1.005, a_cr, mGeV, z_cr))
		y_tmp = inter.spline(xk_log, parIS, ekn_tmp, order=3, kind='smoothest')
		B = (y_tmp - parIS[0]) / (ekn_tmp - xk_log[0])
		logA = parIS[0] - (B * xk_log[0])
		result = logA + B * ekn_log

	elif ekn_log > xk_log_max :
		ekn_tmp = np.log10(R_to_Ekn(799.995, a_cr, mGeV, z_cr))
		y_tmp = inter.spline(xk_log, parIS, ekn_tmp, order=3, kind='smoothest')
		B = (parIS[5] - y_tmp) / (xk_log[5] - ekn_tmp)
		logA = parIS[5] - (B * xk_log[5])
		result = logA + B * ekn_log

	return np.power(10, result)

#--------------------------------

def Force_Field(Edata, par_IS, par_mod, name_exp_CR) :

	a_cr, z_cr, mGeV = CREntry(name_exp_CR)
	AZ = float(a_cr)/float(z_cr)
	flux_TOA = []

	for i in range(0, len(Edata)):
		flux_TOA.append((np.power(Edata[i],2)*a_cr + 2*mGeV*Edata[i]) / (np.power((Edata[i]+par_mod/AZ),2)*a_cr + 2*mGeV*(Edata[i]+par_mod/AZ)) * IS_Spline_flux(Edata[i]+par_mod/AZ, par_IS, name_exp_CR))
		#flux_IS2.append(IS_Spline_flux(Edata[i], par_IS, name_exp_CR))
		#flux_IS1.append(IS_Spline_flux(Edata[i]+z_cr*par_mod, par_IS, name_exp_CR))
	#for i in range(0, len(E_IS)):	
		#flux_IS2.append(IS_Spline_flux(E_IS[i], par_IS, name_exp_CR))
		#print "Etoa : ", Edata[i], "spline : ", IS_Spline_flux(Edata[i]+par_mod, par_IS, name_exp_CR)

	'''
	fig_all = plt.figure(figsize=(12,8))
	fig_all.set_facecolor('white')
	plot_IS1 = plt.plot(Edata+z_cr*par_mod, flux_IS1*np.power(Edata+z_cr*par_mod, 2.), c='black', lw = 2.5, label='Interstellar flux')
	#plot_IS2 = plt.plot(E_IS, flux_IS2*np.power(E_IS, 2.), c='black', lw = 2.5, label='Interstellar flux')
	#plot_IS2 = plt.plot(Edata, flux_IS2*np.power(Edata, 2.), c='black', lw = 2.5, label='Interstellar flux')
	plot_TOA = plt.plot(Edata, flux_TOA*np.power(Edata, 2.), c='red', lw = 2.5, label='TOA flux')
	plt.xscale('log')
	plt.yscale('log')
	plt.show()
	sys.exit()
	'''

	return np.array(flux_TOA)

#--------------------------------

def Momen(R, E0, AZ, par_IS, exp_CR):

	beta = R / sqrt(R*R + E0*E0*AZ*AZ)
  	KE   = beta*R          					# E-dependent diffusion coefficient 
  	T    = R/(beta * AZ) - E0      			# => T=Ekn, so that 'AZ'=A/Z

 	# dJ_IS/dEkn: GCR LIS from Webber and Lockwood [2001] in particles/m 2 .sr.s.MeV/n,
  	# with Ek=log(T) the kinetic energy in MeV/n.
  	#FLIS = 21.1*pow(T, -2.8) / (1 + 5.85*pow(T, -1.22) + 1.18*pow(T, -2.54))

  	FLIS = IS_Spline_flux(T, par_IS, exp_CR)
	return beta, KE, T, FLIS

#--------------------------------

def Model1D(Edata_exp, par_IS, par_K0, name_exp_CR):

	file_name = 'results/BestFit1D.txt'
	results = open(file_name, "w")

	# Cosmic ray properties
	a_cr, z_cr, mGeV = CREntry(name_exp_CR)
	E0 = mGeV / float(a_cr)
	AZ = float(a_cr)/float(z_cr)

	# Space grid
	N = 91
	dr = 1.0
	
	# Energy grid
	R 		= Ekn_to_R(np.amax(Edata_exp), a_cr, E0, z_cr)
	Rmin	= Ekn_to_R(np.amin(Edata_exp), a_cr, E0, z_cr)
	dlnR	= 0.02												# Optimized value dlnR   = 0.02

	# Lists and arrays initialization
	V0 	= 1.e3 * 400. / Convert_AUperyr_to_mpers()			# Solar wind V(r) = V0 [AU/yr]
	#Rd	= 29.
  	K, V, r, dKdr, dV, F, X, Xn, Y = (np.zeros(N+1) for _ in xrange(9)) 
	Rgrid, flux_TOA = ([] for _ in xrange(2))

	#---------------------------------------------
	
	beta, KE, T, spectrum = Momen(R, E0, AZ, par_IS, name_exp_CR)
	r[0] = -dr + 0.000001

  	# Fills r-dependent terms
  	for i in range(1,N+1):    
		F[i]  	= spectrum / (R*R*AZ)            	# Assuming (Eq.7 of Moraal) that f_p=(dJ/dEkn)/p^2, we have f_R=dJ/dEkn*A/(R^2*Z) (we recall 'AZ'=A/Z)
		r[i]  	= r[i-1] + dr                       # => Position r
		V[i] 	= V0								# V(r) = cst = V0
		dV[i]	= 2.0 * V[i] / r[i]					# 2V/r + dV/dr
		K[i] 	= par_K0 							# K(r, E) = K0*beta*R  => K(r) = K0
		dKdr[i] = 0. 
		#K[i] = K0 * exp((r[i]-1.)/Rd)
		#dKdr[i] = K[i]/Rd
  	  		
	# R = Rmax case	
  	Rgrid.append(R)
  	flux_TOA.append(spectrum)  	

	while R >= Rmin :
	
		R = R/exp(dlnR)	
		Rgrid.append(R/exp(dlnR/2.0))   # TEST
		F[0]	= F[2]                      	       	# Derivative at i_r=1 is 0 for all energies
		Xn[1] 	= -1.0                            		# [????] => all other XN[i] appear to be equal to 0
		beta, KE, T, spectrum = Momen(R/exp(dlnR/2.0), E0, AZ, par_IS, name_exp_CR)
			
  		# Loop on r: calculate CN coefficient and invert tridiagonal matrix on the fly
  		for i in range(1, N):
      			
			A = KE * (K[i]/(2.*dr*dr) - (2.0*K[i]/r[i]+dKdr[i])/(4.*dr)) + V[i]/(4.*dr)
			B = KE * K[i]/(dr*dr) - A
			C = - A - B - A*X[i-1] - dV[i]/(3. * dlnR)
      		#print "(%.3f,%.3f,%.3f) " % (A, B, C)
			
			X[i]=(B - A*Xn[i]) / C
			Y[i]=(-A*F[i-1] - (dV[i]/(3.*dlnR) - A - B) * F[i] - B*F[i+1] - A*Y[i-1]) / C	
  	
		F[N-1] = Y[N-1] - X[N-1] * spectrum / exp(2.0*log(R/exp(dlnR/2.0))) / AZ
		beta, KE, T, spectrum = Momen(R/exp(dlnR), E0, AZ, par_IS, name_exp_CR)
		F[N] = spectrum / exp(2.0*log(R/exp(dlnR))) / AZ
		for i in range(N-2,0, -1):
			F[i] = Y[i] - X[i]*F[i+1]
  		
		DP2 = exp(2.0*log(R/exp(dlnR)))
		
		flux_TOA.append(DP2 * AZ *F[2])	# r = 1 AU = r[2]
	
		# Write results if best fit
		#if (res == 1):
		#	results.write("%10.3e \t %10.3e \t %10.3e" % (R, T*1000, DP2 * AZ *F[2]))     	# R en GV, T en MeV
		#	for i in range(1,10) : 
		#		results.write("%10.3e" % (DP2*AZ*F[10*i+1]))      							# dJ/dT en #/m2/s/sr/MeV
		#	results.write("%10.3e \n" % (spectrum))
			
	Ekn_grid = R_to_Ekn(np.array(Rgrid), a_cr, E0, z_cr)
	flux_TOA = np.array(flux_TOA)
	#print "Ekn = \n", Ekn_grid, "flux_TOA (1AU) = \n", flux_TOA
	
	results.close()
	return Ekn_grid, flux_TOA

#--------------------------------

def Chi21D(init_pars, MODE, Edata_exp, flux_data_exp, sigma_exp, Nexp, Ndata, N_IS, list_CR, list_exp_CR):

	result = []    
	par_IS  = init_pars[0:N_IS]
	#print "par_IS : ", par_IS

	if (MODE == "FF"):
		par_mod = init_pars[N_IS:N_IS+Nexp]
	elif (MODE == "1D"):
		par_mod = np.power(10,init_pars[N_IS:N_IS+Nexp])
	
	ndof = 0
	Ncr = N_IS / len(list_CR)

	for i in range(0, Nexp):

		for k in range(0, len(list_exp_CR[i])):

			index = list_CR.index(list_exp_CR[i][k])
			par_IS_CR = par_IS[Ncr*index:Ncr+Ncr*index]
			ndof += Ndata[i][k]

			if (MODE == "FF") :
				flux_TOA = Force_Field(np.array(Edata_exp[i][k]), par_IS_CR, par_mod[i], list_exp_CR[i][k])
				result.append(np.sum(np.power((flux_TOA - flux_data_exp[i][k]) / sigma_exp[i][k] , 2)))

			elif (MODE == "1D"):
				Ekn_grid, flux_TOA = Model1D(np.array(Edata_exp[i][k]), par_IS_CR, par_mod[i], list_exp_CR[i][k])	
				Ekn_grid = np.flipud(Ekn_grid) 		# Revert array to have increasing value to make interpolation
				flux_TOA = np.flipud(flux_TOA)
				flux_TOA_interp = np.interp(Edata_exp[i][k], Ekn_grid, flux_TOA)
				result.append(np.sum(np.power((flux_TOA_interp - flux_data_exp[i][k]) / sigma_exp[i][k] , 2)))

	ndof -= len(init_pars)
	result = np.array(result)
	chi2_red = np.sum(result) / ndof
	return chi2_red

#--------------------------------

def Chi2_best1D(Best_IS, Best_mod, Edata_exp, flux_data_exp, sigma_exp, Nexp, Ndata, N_IS, list_CR , list_exp_CR):

	result = []
	chi2_red_exp = []   
	ndof = 0

	Ncr = N_IS / len(list_CR)
	
	for i in range(0, Nexp):
		chi2_exp_tmp = []

		for k in range(0, len(list_exp_CR[i])):

			index = list_CR.index(list_exp_CR[i][k])
			Best_IS_CR = Best_IS[Ncr*index:Ncr+Ncr*index]
			ndof += Ndata[i][k]

			if (MODE =="FF") :
				flux_TOA = Force_Field(np.array(Edata_exp[i][k]), Best_IS_CR, Best_mod[i], list_exp_CR[i][k])
				chi2_exp_tmp.append(np.sum(np.power((flux_TOA - flux_data_exp[i][k]) / sigma_exp[i][k] , 2)) / Ndata[i][k])
				result.append(np.sum(np.power((flux_TOA - flux_data_exp[i][k]) / sigma_exp[i][k] , 2)))

			elif (MODE == "1D"):
				Ekn_grid, flux_TOA = Model1D(np.array(Edata_exp[i][k]), Best_IS_CR, np.power(10,Best_mod[i]), list_exp_CR[i][k])
				Ekn_grid = np.flipud(Ekn_grid) 		# Revert array to have increasing value to make interpolation
				flux_TOA = np.flipud(flux_TOA)
				flux_TOA_interp = np.interp(Edata_exp[i][k], Ekn_grid, flux_TOA)
				chi2_exp_tmp.append(np.sum(np.power((flux_TOA_interp - flux_data_exp[i][k]) / sigma_exp[i][k] , 2)) / Ndata[i][k])
				result.append(np.sum(np.power((flux_TOA_interp - flux_data_exp[i][k]) / sigma_exp[i][k] , 2)))
			
		chi2_red_exp.append(chi2_exp_tmp)

	ndof -= (N_IS + Nexp)
	result = np.array(result)
	chi2_red = np.sum(result) / ndof
	print "Chi2_global = ", chi2_red
	return chi2_red, chi2_red_exp

#--------------------------------

def BestParams(MODE, Nexp, pars):

	Best_IS  = pars[0 : N_IS]
	Best_mod = pars[N_IS : N_IS+Nexp]

	if (MODE == "FF"):
		return Best_IS, Best_mod, Best_mod

	elif (MODE == "1D"):
		Best_phi = np.zeros(Nexp)
		for i in range(0,Nexp):
			Best_phi[i] = 1.e3 * 400. / Convert_AUperyr_to_mpers() / (3 * np.power(10,Best_mod[i])) * (90. - 1.)
		return Best_IS, Best_mod, Best_phi

#--------------------------------

def Best_Fit1D(MODE, E_IS, E_TOA, Best_IS, Best_mod, list_CR, list_exp_CR):

	flux_IS = []
	Ncr = len(Best_IS) / len(list_CR)
	print "Best_IS",Best_IS

	for i in range(0, len(list_CR)) :
		flux_IS_tmp = []
		Best_IS_CR = Best_IS[Ncr*i:Ncr+Ncr*i]
		print "Best_IS_CR", Best_IS_CR
		for k in range(0,len(E_IS)) :
			flux_IS_tmp.append(IS_Spline_flux(E_IS[k], Best_IS_CR, list_CR[i]))
		flux_IS.append(flux_IS_tmp)
	print E_IS, flux_IS

	flux_TOA = []   
	for i in range(0, Nexp) :
		flux_TOA.append([])

	all_TOA = []
	for i in range(0, len(list_CR)):
		all_TOA.append([])
		all_TOA[i].append(E_TOA.tolist())


	for i in range(0, Nexp): 
		flux_TOA_tmp2 = []

		for k in range(0, len(list_exp_CR[i])) :

			index = list_CR.index(list_exp_CR[i][k])
			Best_IS_CR = Best_IS[Ncr*index:Ncr+Ncr*index]

			if (MODE == "FF"):
				flux_TOA_FF = Force_Field(np.array(E_TOA), Best_IS_CR, Best_mod[i], list_exp_CR[i][k])
				flux_TOA[i].append(flux_TOA_FF)
				all_TOA[index].append(flux_TOA_FF)

			elif (MODE == "1D"):	
				Ekn_grid, flux_TOA_tmp = Model1D(np.array(E_TOA), Best_IS_CR, np.power(10,Best_mod[i]), list_exp_CR[i][k])
				print Ekn_grid, flux_TOA_tmp
				Ekn_grid = np.flipud(Ekn_grid)
				flux_TOA_tmp = np.flipud(flux_TOA_tmp)
				flux_TOA_interp = np.interp(E_TOA, Ekn_grid, flux_TOA_tmp)
				flux_TOA_interp = flux_TOA_interp.tolist()
				flux_TOA[i].append(flux_TOA_interp)
				all_TOA[index].append(flux_TOA_interp)

	#for i in range(0, len(list_CR)):
	#	flux_TOA[i] = np.array(flux_TOA[i])

	#print E_TOA, flux_TOA
	#np.savetxt('test.txt', (np.transpose(all_TOA)))#,fmt='%1.4e')
	return flux_IS, flux_TOA

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#						  	      								  CODE PRINCIPAL																#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

#--------------------------------
# Choose your modulation model : 
#		FF = force-field 
#		1D = 1D spherical model
#--------------------------------

MODE = "1D"

#--------------------------------
# Data extraction :  
# - ExtractData(list_exp) : From a list of experiments list_exp = ['AMS02','BESS00', ...] and a list of CR list_CR = ['H','He','Electron',...]
#                           Return - Nexp, Ndata, Edata_exp, flux_data_exp, sigma_exp, Edata, flux_data, sigma, date_list_mean, date_list_delta 
#--------------------------------

list_CR = ['He']
#list_exp = RefDataset()
list_exp = ['AMS02', 'BESS97', 'BESS00']
#list_exp = ['AMS01','AMS02', 'BESS97', 'BESS98', 'BESS99', 'BESS00', 'BESSPOLAR1', 'BESSPOLAR2', 'PAMELA0608']
Nexp, Ndata, Edata_exp, flux_data_exp, sigma_exp, Edata, flux_data, sigma, date_list_mean, date_list_delta, list_exp_CR = ExtractData(list_exp, list_CR)

# Voyager data
#Ndata_Voy, Edata_Voy, ydata_Voy, sigma_Voy = ExtractExp('VOYAGER1')

#--------------------------------------
# Define parameters for minimization :  
# - N_IS : nb of parameters for the interstellar flux
# - Nexp : nb of parameters for K0 (= nb of experiments, already define)
# - InitParameters(Nexp) : Fill init_pars list with [IS flux default parameters, ..., ..., and Nexp time phi0]
#--------------------------------------
'''

Emin_IS = 0.126
Emax_IS = 5.e3

#Best_IS = [3.8486011132, 1.93403494265, -0.471481132898, -1.3384898445, -3.02197021307, -3.86190526353] # H 
Best_IS = [4.28655743616, 1.04497401525, -0.961104403387, -1.75375209438, -3.34606895311, -4.14137773798]
Best_mod = [0.623936476681, 0.528795462052, 0.92673183789, 0.534770467312]
#Best_IS = [1.8,1.,-1.,-2.3,-3.8,-4.4] 

E_IS = np.logspace(log10(Emin_IS), log10(Emax_IS), num = 150)
E_TOA = np.logspace(log10(min(Edata)), log10(max(Edata)), num = 150)
flux_IS, flux_TOA = Best_Fit1D(MODE, E_IS, E_TOA, Best_IS, Best_mod, list_CR, list_exp_CR)
print flux_TOA

#fig_all = plt.figure(figsize=(12,8))
#fig_all.set_facecolor('white')
#plot_IS = plt.plot(E_IS, flux_IS[0]*np.power(E_IS, 2.), c='black', lw = 2.5, label='Interstellar flux')

for i in range(0, len(list_CR)):
	fig_all = plt.figure(figsize=(12,8))
	fig_all.set_facecolor('white')
	plot_IS = plt.plot(E_IS, flux_IS[0]*np.power(E_IS, 2.), c='black', lw = 2.5, label='Interstellar flux')

	for j in range(0,Nexp):   # Data and TOA fluxes best fits
		plot_exp = plt.errorbar(Edata_exp[j][0], flux_data_exp[j][0]*np.power(Edata_exp[j][0], 2.), xerr = None, yerr=sigma_exp[j][0]*np.power(Edata_exp[j][0], 2), fmt='o',ms=3., label=list_exp[j])
		col = plot_exp[0].get_color()
		#plot_TOA = plt.plot(E_TOA, flux_TOA[j][0]*np.power(E_TOA, 2.), color = col, ls ='-',ms=3.)

plt.xscale('log')
plt.yscale('log')
plt.show()

sys.exit()
'''

init_pars, N_IS = InitParameters1D(MODE, Nexp, list_CR) #par_IS, par_mod, 
bnds = SetBounds(MODE, N_IS, Nexp)
pars = fmin_slsqp(Chi21D, init_pars, bounds=bnds, args=(MODE, Edata_exp, flux_data_exp, sigma_exp, Nexp, Ndata, N_IS, list_CR, list_exp_CR), acc=1e-03, iprint=1, iter=1000, full_output=0, epsilon=1.49e-08) 
Best_IS, Best_mod, Best_phi = BestParams(MODE, Nexp, pars)

chi2_red, chi2_red_exp = Chi2_best1D(Best_IS, Best_mod, Edata_exp, flux_data_exp, sigma_exp, Nexp, Ndata, N_IS, list_CR, list_exp_CR)
stderr = np.zeros(N_IS+Nexp)
#stderr = Dispersion(list_exp)

# Print best fit results	    
PrintResults1Db(N_IS, Nexp, list_exp, list_CR, list_exp_CR, Best_IS, Best_phi, chi2_red_exp, stderr, chi2_red)

# Recover best fits values to plot results
# - For interstellar flux use : IS_Spline_flux(E, best_IS)
# - For toa fluxes use Best_Fit1D(E, best_IS, best_phi, Nexp) which returns flux_TOA "per block"


Emin_IS_tab, Emax_IS_tab =  ([] for _ in xrange(2))
for i in range(0, Nexp):
	for k in range(0, len(list_exp_CR[i])) :
		Emin_IS_tab.append(min(Edata_exp[i][k]+Best_phi[i]))
		Emax_IS_tab.append(max(Edata_exp[i][k]+Best_phi[i]))
#print min(Emin_IS_tab)
#Emax_IS = max(Emax_IS_tab)

#Emin_IS = 0.1
#Emax_IS = 1.5e3
Emin_IS = 0.1
#Emin_IS = min(Emin_IS_tab)
Emax_IS = 5.e3

E_IS = np.logspace(log10(Emin_IS), log10(Emax_IS), num = 150)
E_TOA = np.logspace(log10(min(Edata)), log10(max(Edata)), num = 150)

#Best_IS = [ 5., 1.6202, -0.8229, -1.6762, -3.3307, -4.1318]

flux_IS, flux_TOA = Best_Fit1D(MODE, E_IS, E_TOA, Best_IS, Best_mod, list_CR, list_exp_CR)

#---------------------
# Plot all results
#---------------------


# ONE PLOT / EXPERIMENT
#-----------------------
'''
for i in range(0, Nexp):
	fig = plt.figure(i, figsize=(12,8))
	fig.set_facecolor('white')
	plot_IS = plt.plot(E_IS, flux_IS*np.power(E_IS, 2.), c='black', lw = 2.5, label='Interstellar flux')    # Interstellar flux best fit  
	#plot_IS  = plt.errorbar(E_IS, flux_IS*np.power(E_IS, 2), xerr = None, yerr=np.power(10,stderr[0:7])*np.power(E_IS, 2),  c='black', lw = 2.5, label='Interstellar flux') 
	plot_exp = plt.errorbar(Edata_exp[i], flux_data_exp[i]*np.power(Edata_exp[i], 2.), xerr = None, yerr=sigma_exp[i]*np.power(Edata_exp[i], 2), fmt='bo',ms=3., label=list_exp[i] + ' data')
	plot_TOA = plt.plot(E_TOA, flux_TOA[i]*np.power(E_TOA, 2.), color='r', ls='-',ms=3., label=list_exp[i]+ ' TOA flux' + (', $\phi$ = %.2f' % best_phi[i]))

	SetAxis(r"E[GeV]", r"$\rm{Flux [GeV^{-1}.m^{-2}.s^{-1}.sr^{-1}]}$", min(Edata), 1e3, 0.1, 5.e5, "xylog")
	SetTitle(r"Global Fit")
	SetLegend()
	SaveFig("plots/all", list_exp[i], "png")
	#plt.close()
'''

# ALL EXPERIMENTS
#-----------------
'''
fig1 = plt.figure(30, figsize=(12,8))
fig1.set_facecolor('white')
plot_IS = plt.plot(E_IS, flux_IS*np.power(E_IS, 2.), c='black', lw = 2.5, label='Interstellar flux')
#plt.fill_between(E_IS, (flux_IS - np.power(10,stderr[0:7]))*np.power(E_IS, 2), (flux_IS + np.power(10,stderr[0:100]))*np.power(E_IS, 2), alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
#plot_IS  = plt.errorbar(E_IS, flux_IS*np.power(E_IS, 2), xerr = None, yerr=stderr[0:50]*np.power(E_IS, 2),  c='black', lw = 2.5, label='Interstellar flux')   
#plot_IS = plt.plot(E_IS, flux_IS*np.power(E_IS, 2), c='black', lw = 2.5, label='Interstellar flux')    # Interstellar flux best fit
#plt.fill_between(E_IS, flux_IS-(stderr[0:50]*np.power(E_IS, 2.8)), flux_IS+(stderr[0:50]*np.power(E_IS, 2.8)),alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')

for i in range(0,Nexp):   # Data and TOA fluxes best fits
	#if (i%3==0):
	plot_exp = plt.errorbar(Edata_exp[i], flux_data_exp[i]*np.power(Edata_exp[i], 2.), xerr = None, yerr=sigma_exp[i]*np.power(Edata_exp[i], 2), fmt='o',ms=3., label=list_exp[i] + ' data')
	col = plot_exp[0].get_color()
	plot_TOA = plt.plot(E_TOA, flux_TOA[i]*np.power(E_TOA, 2.), color = col, ls ='-',ms=3., label=list_exp[i]+ ' TOA flux' + (', $\phi$ = %.2f' % best_phi[i]) + 'GV')

	#plot_Voy = plt.errorbar(Edata_Voy, ydata_Voy*np.power(Edata_Voy, 2.), xerr = None, yerr = sigma_Voy*np.power(Edata_Voy, 2.), fmt='o',ms=3., label= 'Voyager1 data')


SetAxis(r"E[GeV]", r"$\rm{Flux [GeV^{-1}.m^{-2}.s^{-1}.sr^{-1}]}$", min(Edata), 1e4, 0.1, 5e5, "xylog")
SetTitle(r"Global Fit")
SetLegend(ncolumn = 2)
#SaveFig("plots/all", "all_exp", "png")
'''

#Un graphe par espece
for i in range(0, len(list_CR)):
	fig_all = plt.figure(figsize=(12,8))
	fig_all.set_facecolor('white')
	plot_IS = plt.plot(E_IS, flux_IS[i]*np.power(E_IS, 2.), c='black', lw = 2.5, label='Interstellar flux')

	for j in range(0,Nexp):   # Data and TOA fluxes best fits
		if list_CR[i] in list_exp_CR[j] :
			index = list_exp_CR[j].index(list_CR[i])
			plot_exp = plt.errorbar(Edata_exp[j][index], flux_data_exp[j][index]*np.power(Edata_exp[j][index], 2.), xerr = None, yerr=sigma_exp[j][index]*np.power(Edata_exp[j][index], 2), fmt='o',ms=3., label=list_exp[j] + (', $\phi$ = %.2f' % Best_phi[j]) + 'GV')
			col = plot_exp[0].get_color()
			plot_TOA = plt.plot(E_TOA, flux_TOA[j][index]*np.power(E_TOA, 2.), color = col, ls ='-',ms=3.)

	SetAxis(r"$\rm{E_{k/n}[GeV/n]}$", r"$\rm{J \times E^2_{k/n} [GeV/n.m^{-2}.s^{-1}.sr^{-1}]}$", min(Edata), 1e4, 0.1, 5e5, "xylog")
	SetTitle(r"Global Fit")
	SetLegend(ncolumn = 2)

plt.show()

# PHI EN FONCTION DU TEMPS
#---------------------------       
'''
fig2, ax = plt.subplots(1,1,figsize=(12,8))
fig2.set_facecolor('white')

for i in range(0, Nexp):
	if 'AMS' in list_exp[i]:
		col = 'navy'
	elif 'BESS' in list_exp[i]:
		col = 'indianred'
	elif 'PAMELA' in list_exp[i]:
		col = 'darkorange'
	else :
		col = 'steelblue'
	ax.errorbar(date_list_mean[i], best_phi[i], xerr=date_list_delta[i], fmt='o', color = col)   #, yerr=stderr[7+i:7+i+1]

for i in range(0,Nexp):
	ax.annotate(list_exp[i], (date_list_mean[i],best_phi[i]+0.05))

SetDateAxis(ax, "yearly")
SetAxis(r"Dates", r"$\phi \rm{[GV]}$", datetime.date(1992, 1, 1), datetime.date(2018, 1, 1), 0.2, 1.5, "lin")
SetTitle(r"Best $\phi$ values depending on time")
#SaveFig("plots/all", "phi_date", "png")


# PHI EN FONCTION DU CHI2       
fig3, ax = plt.subplots(1,1,figsize=(12,8))
fig3.set_facecolor('white')
plot_chi2 = plt.plot(chi2_red_exp, best_phi, c='black', marker='o', ls='')  
for i in range(0,Nexp):
	ax.annotate(list_exp[i], (chi2_red_exp[i],best_phi[i]+0.01))

SetAxis(r"$\chi^2$", r"$\phi \rm{[GV]}$", "None", "None", 0.2, 1.5, "lin")
SetTitle(r"Best $\phi$ values vs $\chi^2$")
#SaveFig("plots/all", "phi_chi2", "png")
	'''
	
	 

    
#---------------------

#start_time = time.time()
#main()
#print("--- %s seconds ---" % (time.time() - start_time))
#sys.exit() 

