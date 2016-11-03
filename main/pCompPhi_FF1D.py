import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, '/ams/aupetit/Documents/Pheno/Solmod/Lib/')
from Physics import *
from ExtractData import *
from PlotLib import *

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

#---------------------------------------------------------------------------------------------------

#list_exp = ['AMS01','AMS02', 'BESS97', 'BESS98']
list_exp = ['AMS01','AMS02', 'BESS97', 'BESS98', 'BESS99', 'BESS00', 'BESSPOLAR1', 'BESSPOLAR2', 'PAMELA0608', 'PAMELA2006', 'PAMELA2007', 'PAMELA2008', 'PAMELA2009']

npoint = 3
x = np.linspace(0., 2*len(list_exp), (len(list_exp)*npoint)+(len(list_exp)+1))

x1 = np.linspace(0.5, 0.5 + 2*(len(list_exp)-1), len(list_exp))
x2 = np.linspace(1. , 1.  + 2*(len(list_exp)-1), len(list_exp))
x3 = np.linspace(1.5, 1.5 + 2*(len(list_exp)-1), len(list_exp))

phi_FF_H = [0.564957654925, 0.726737856292, 0.498487216343, 0.596254826545, 0.673512788975, 1.31584562575 , 0.790936376939, 0.544054211101, 0.55730555504 , 0.633774773706, 0.589521076594, 0.519937648699, 0.47121552565 ]
phi_FF_He = [0., 1.10168215347 , 0.918966529971, 1.00815628958 , 1.09793092422 , 1.69758238841 , 1.13890900751 , 0.881098414174, 0., 0., 0., 0., 0.]
phi_FF_all = [0.628386706808, 0.798026332239, 0.585202908893, 0.680543921225, 0.761330967389, 1.38867466922, 0.851981713583, 0.601589512122, 0.619830507582, 0.697664946201, 0.652011475938, 0.582890601823, 0.533726782632]

phi_1D_H = [0.503081220333, 0.673363148176, 0.431587252929, 0.535710166793, 0.615901689143, 1.29697600794 , 0.742779956288, 0.479001239484, 0.494018254226, 0.567314281567, 0.521252117585, 0.449842612678, 0.399972026465]
phi_1D_He = [0., 0.96507993653, 0.76427892272, 0.86531566282, 0.96174097592, 1.6219526231 , 1.00656011965, 0.720803183058, 0., 0., 0., 0., 0.]
phi_1D_all = [0.330402612933, 0.507718143154, 0.285995984662, 0.384841582776, 0.469292294974, 1.12391792438 , 0.563540703536, 0.302858329392, 0.32248783652 , 0.3963113928  , 0.352259265436, 0.28392059071 , 0.237044657512]


# Three subplots sharing both x/y axes
f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True, figsize=(24,12))
f.set_facecolor('white')

ax1.errorbar(x1,phi_FF_H, xerr = None, yerr=0.050, fmt='bo',ms=8., label = "H")
ax1.errorbar(x2,phi_FF_He, xerr = None, yerr=0.050, fmt='ro',ms=8., label = "He")
ax1.errorbar(x3,phi_FF_all, xerr = None, yerr=0.050, fmt='ko',ms=8., label = "H + He")

ax2.errorbar(x1,phi_1D_H, xerr = None, yerr=0.050, fmt='bs',ms=8., label = "H")
ax2.errorbar(x2,phi_1D_He, xerr = None, yerr=0.050, fmt='rs',ms=8.,  label = "He")
ax2.errorbar(x3,phi_1D_all, xerr = None, yerr=0.050, fmt='ks',ms=8., label = "H + He")

my_xticks = ['']
for i in range (0, len(list_exp)):
	my_xticks.append('') ; my_xticks.append(list_exp[i]); 
	if i == (len(list_exp)-1):
		my_xticks.append('')
	else :
		my_xticks.append(''), my_xticks.append('')
plt.xticks(x, my_xticks)
plt.xticks(fontsize=24, rotation=50)


ax1.tick_params(axis='y', labelsize=24)
ax2.tick_params(axis='y', labelsize=24)

plt.xlim(0., 2*len(list_exp))
ax1.set_ylim([0.12, 1.8])
ax2.set_ylim([0.12, 1.8])
ax1.set_ylabel(r'$\phi$ [GV]',fontsize=25, labelpad=20)
ax2.set_ylabel(r'$\phi$ [GV]',fontsize=25, labelpad=20)

ax1.legend(ncol = 1, numpoints = 1,frameon=True, fontsize = 30) #prop={'size':30}
ax2.legend(ncol = 1, numpoints = 1,frameon=True, fontsize = 30)

ax1.annotate(r"\textbf{FF}", fontsize=30, xy=(0.025, 0.90), xycoords='axes fraction')
ax2.annotate(r"\textbf{1D}", fontsize=30, xy=(0.025, 0.90), xycoords='axes fraction')


xline = np.linspace(2.,2*(len(list_exp)-1), (len(list_exp)-1))
for i in range(0,xline.size) :
	ax1.axvline(x=xline[i], ymin=0., ymax = 2., ls='--', linewidth=1, color='k')
	ax2.axvline(x=xline[i], ymin=0., ymax = 2., ls='--', linewidth=1, color='k')


# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
plt.tight_layout()
f.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
plt.show()




