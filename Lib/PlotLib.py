#!/usr/bin/python2.7
# -*-coding:Utf-8 -*

import matplotlib.pyplot as plt
from matplotlib.dates import YearLocator, YEARLY, MonthLocator, RRuleLocator, rrulewrapper, DateFormatter
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

#---------------------------------------------------------------
def SetFig ():
	fig = plt.figure(figsize=(10,8))
	fig.set_facecolor('white')
	return fig

#---------------------------------------------------------------
def SetTitle(title):
	plt.title(title,fontsize = 20)
	plt.tight_layout()

#---------------------------------------------------------------
def SetAxis(xlabel, ylabel, xmin, xmax, ymin, ymax, scale):

	plt.xlabel(xlabel,fontsize=30, labelpad=20)
	plt.ylabel(ylabel,fontsize=30, labelpad=20)
	plt.xticks(fontsize=26) ; plt.yticks(fontsize=26)

	if "x" in scale:
		plt.xscale('log')
	if "y" in scale :
		plt.yscale('log')

	if xmin and xmax != "None":
		plt.xlim(xmin, xmax) 
	if ymin and ymax != "None":
		plt.ylim(ymin, ymax) 

	plt.tight_layout()

#---------------------------------------------------------------
def SetAxisSub(ax1, ax2, xlabel, ylabel, xmin, xmax, ymin, ymax, sharexy, scale):

	ax2.set_xlabel(xlabel,fontsize=30)
	ax1.set_ylabel(ylabel,fontsize=30)  
	ax2.set_ylabel(ylabel,fontsize=30)  

	ax2.tick_params(axis='x', labelsize=26)
	ax1.tick_params(axis='y', labelsize=26)
	ax2.tick_params(axis='y', labelsize=26)

	if not "x" in sharexy :
		ax1.set_xlabel(xlabel,fontsize=30)
		ax1.tick_params(axis='x', labelsize=26)
	
	if "x" in scale:
		plt.xscale('log')
	if "y" in scale :
		plt.yscale('log')

	if xmin and xmax != "None":
		plt.xlim(xmin, xmax) 
	if ymin and ymax != "None":
		plt.ylim(ymin, ymax) 

	plt.tight_layout()

#---------------------------------------------------------------
def SetDateAxis(ax, form):
	if "year" in form :
		formatter = DateFormatter('%Y')
		rule = rrulewrapper(YEARLY, byeaster=1, interval=5)
		loc = RRuleLocator(rule)
	if "month" in form :
		formatter = DateFormatter('%Y-%m')
		loc = MonthLocator(interval = 4)
		
	ax.xaxis.set_major_locator(loc)
	ax.xaxis.set_major_formatter(formatter)
	plt.xticks(rotation=20)
	plt.tight_layout()

#---------------------------------------------------------------
def SetLegend(localisation = 1, ncolumn = 1):
	plt.legend(loc = localisation, ncol = ncolumn, numpoints = 1,frameon=False)
	leg = plt.gca().get_legend()
 	ltext  = leg.get_texts()            
	plt.setp(ltext, fontsize='xx-large')    # the legend text fontsize 
	plt.tight_layout()

#---------------------------------------------------------------
def SetLegendSub(ax, ncolumn = 1):
	ax.legend(ncol = ncolumn, numpoints = 1,frameon=False)
	plt.tight_layout()

#---------------------------------------------------------------
def SaveFig(dir, name, form):
	st = dir + "/" + name + "." + form
	st_eps = dir + "/eps/" + name + "." + "eps"
	plt.savefig(st, dpi=None, facecolor='w', edgecolor='w', format=form, transparent = True, bbox_inches="tight")
	plt.savefig(st_eps, dpi=None, facecolor='w', edgecolor='w', format='eps', transparent = True, bbox_inches="tight")

