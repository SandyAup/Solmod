#!/usr/bin/python2.7
# -*-coding:Utf-8 -*

import os
import glob
import re
import numpy as np
from math import *
import datetime
from datetime import date
import string

import sys
import os.path


#--------------------------------
# Extraction des donnees 
#--------------------------------

def RefDataset():
    list_exp = ['AMS01','AMS02', 'BESS97', 'BESS98', 'BESS99', 'BESS00', 'BESSPOLAR1', 'BESSPOLAR2', 'PAMELA0608', 'PAMELA2006', 'PAMELA2007', 'PAMELA2008', 'PAMELA2009']
    return list_exp

#--------------------------------

def RefDatasetHe():
    list_exp = ['AMS02', 'BESS97', 'BESS98', 'BESS99', 'BESS00', 'BESSPOLAR1', 'BESSPOLAR2'] 
    return list_exp

#--------------------------------

def PrintInfos(list_CR, list_exp, MODE, powr):
    print "\n\t Performing analysis"
    print "\t ------------------- \n" 
    print "\t \t- Cosmic ray species \t\t:", list_CR
    print "\t \t- Experiments \t\t\t:", list_exp
    print "\t \t- Modulation model \t\t:", MODE, "model"
    if (MODE == "1D") : print "\t \t- Diffusion coefficient \t: K(r) = K0 * r ^ (", powr,")\n"   
    return 0 

#--------------------------------

def HeaderResultsFile(MODE, powr, list_CR, list_exp):

    head  = "----------------------------------------------------\n"
    head += "\t\t\tAnalysis summary (" + str(date.today()) + ") \n"
    head += "----------------------------------------------------\n" 
    head += " - Cosmic ray species \t: " + " ".join(str(x) for x in list_CR) + "\n"
    head += " - Experiments \t\t\t: " + ", ".join(str(x) for x in list_exp) + "\n"
    head += " - Modulation model \t\t: " + MODE + " model" + "\n"
    if (MODE == "1D") : head += " - Diffusion coefficient \t: K(r) = K0 * r ^ (" +  str(powr) + ")\n"
    head += "----------------------------------------------------\n"  

    return head

#--------------------------------

def HeaderResultsFileCR(MODE, powr, name_CR, list_CR, list_exp):

    list_exp_CR = [] 
    for name_exp in list_exp :
        st = '../Data/' + name_CR + '_data/data_' + name_exp + '.dat'
        if os.path.isfile(st)  :
            list_exp_CR.append(name_exp)

    head  = "----------------------------------------------------\n"
    head += "\t\t\tAnalysis summary (" + str(date.today()) + ") \n"
    head += "----------------------------------------------------\n" 
    head += " - Cosmic ray species \t: " + " ".join(str(x) for x in list_CR) + "\n"
    head += " - Experiments (%s data)\t: " % name_CR  + ", ".join(str(x) for x in list_exp_CR) + "\n"
    head += " - Modulation model \t\t: " + MODE + " model" + "\n"
    if (MODE == "1D") : head += " - Diffusion coefficient \t: K(r) = K0 * r ^ (" +  str(powr) + ")\n"
    head += "----------------------------------------------------\n"  

    return head

#--------------------------------

def DirResultsFile(MODE, list_CR, powr):

    dir = "results/" + MODE + "/" + "_".join(str(x) for x in list_CR) + "/"
    if MODE == "1D" and powr is not 0 :
        dir += "Kdiff/"
    return dir

#--------------------------------

def ExtractExp(name_exp, list_CR): 
# Extract data from one experiment file 

    st = []

    if "AMS" in name_exp and name_exp != "AMS02" and name_exp != "AMS01" :
        st.append('../Data/AMS/data_' + name_exp + '.dat')
    elif "PAMELA" in name_exp and not "PAMELA200" in name_exp and not "PAMELA0608" in name_exp :
        st.append('../Data/PAMELA/data_' + name_exp + '.dat')
    else :
        for i in range(0, len(list_CR)):
            st.append('../Data/' + list_CR[i] + '_data/data_' + name_exp + '.dat')

    Ncr = 0
    list_exp_CR =[]
    for i in range(0, len(list_CR)):

        if os.path.isfile(st[i])  :
            Ncr += 1
            list_exp_CR.append(list_CR[i])

            data_tmp      = np.loadtxt(st[i])
            Edata_tmp     = data_tmp[:,0]
            ydata_tmp     = data_tmp[:,3]
            sigma_tmp     = data_tmp[:,9]
            Ndata_tmp     = Edata_tmp.size

            if Ncr == 1 :
                Edata = [Edata_tmp.tolist()] ; ydata = [ydata_tmp.tolist()] ; sigma = [sigma_tmp.tolist()] ; Ndata = [Ndata_tmp] ;

            else :
                Edata.append(Edata_tmp.tolist()) ; ydata.append(ydata_tmp.tolist()) ; sigma.append(sigma_tmp.tolist()) ; Ndata.append(Ndata_tmp)

    return Ndata, Edata, ydata, sigma, list_exp_CR


#--------------------------------

def ExtractData(list_exp, list_CR, opt):
# Extract data from a list of experiment

    #print '\n\t Loading experiments : ', list_exp, "\n"
    Nexperiment, Ndata, EdataB, ydataB, sigmaB, date_list_mean, date_list_delta, list_exp_CR = ExtractDataBlock(list_exp, list_CR)   

    Edata, ydata, sigma = ([] for _ in xrange(3))
    for i in range(0, Nexperiment) :
        Edata_tmp = np.hstack(EdataB[i]) ; Edata.append(Edata_tmp) 
        ydata_tmp = np.hstack(ydataB[i]) ; ydata.append(ydata_tmp)
        sigma_tmp = np.hstack(sigmaB[i]) ; sigma.append(sigma_tmp)

    Edata = np.array(Edata) ; ydata = np.array(ydata) ; sigma = np.array(sigma)
    Edata = np.hstack(Edata) ; ydata = np.hstack(sigma) ;  sigma = np.hstack(sigma)

    if (opt == "all")   : return Nexperiment, Ndata, EdataB, ydataB, sigmaB, Edata, ydata, sigma, date_list_mean, date_list_delta, list_exp_CR
    if (opt == "light") : return Nexperiment, Ndata, EdataB, ydataB, sigmaB, Edata

#--------------------------------

def ExtractDataBlock(list_exp, list_CR): 
    # Extract data from a list of experiment and put it in an array 
    # Return an "array of arrays"
    
    Nexperiment = len(list_exp)
    Ndata = []

    count = 0    
    for name_exp in list_exp :
        
        if (count==0) : # First experiment
            Nexp, Eexp, yexp, sigmaexp, exp_CR = ExtractExp(name_exp, list_CR)
            Ndata.append(Nexp)
            Edata = [Eexp]
            ydata = [yexp]
            sigma = [sigmaexp]
            list_exp_CR = [exp_CR]
            count += 1
            
        else : # N-th experiment
            Nexp, Eexp, yexp, sigmaexp, exp_CR = ExtractExp(name_exp, list_CR)
            Ndata.append(Nexp)
            Edata.append(Eexp)
            ydata.append(yexp)
            sigma.append(sigmaexp)
            list_exp_CR.append(exp_CR)
            count += 1
            
    Edata = np.array(Edata)     
    ydata = np.array(ydata)
    sigma = np.array(sigma)
    
    date_list_mean, date_list_delta = ExtractDate(list_exp, list_CR)
     
    return Nexperiment, Ndata, Edata, ydata, sigma, date_list_mean, date_list_delta, list_exp_CR

#--------------------------------

def AddMonthlyAMS(list_exp):
    for i in range(1,29):
        list_exp.append("AMS%i" % i)
    return 

#--------------------------------

def AddMonthlyPAMELA(list_exp):
    for i in range(1,46):
        list_exp.append("PAMELA%i" % i)
    return 

#--------------------------------

def PrintResults(N_IS, Nexp, list_exp, best_IS, best_phi, chi2_red, std_error):
    print '\n\t Minimization results : \n'
    
    print '\t\t Interstellar flux parameters : '
    for i in range(0,N_IS):
            #print  '\t\t c%i = ' % i, pow(10,best_IS[i])   , '+/- ', pow(10,std_error[i]), 'diff = ', (pow(10,best_IS[i]) - pow(10,std_error[i]))
            print  '\t\t c%i = ' % i, best_IS[i]   , '+/- ', std_error[i]
    for i in range(0,Nexp):
        if (i<Nexp-1):
            print  '\t\t phi(', list_exp[i], ') = ', best_phi[i], '+/- ', std_error[N_IS + i], ', Chi2 =', chi2_red[i]
        elif (i == Nexp-1) : 
            print '\t\t phi(', list_exp[i], ') = ', best_phi[i], '+/- ', std_error[N_IS + i], ',Chi2 =', chi2_red[i], '\n'
    return 
  
#--------------------------------

def SetBounds(MODE, N_IS, Nexp, min_val, max_val):
    bnds = []
    for i in range(0, N_IS):
        bnds.append((-5, 5))

    if (MODE == "FF"):
        for i in range(0, Nexp):
            bnds.append((min_val, max_val))

    elif (MODE == "1D"):
        for i in range(0, Nexp):
            bnds.append((min_val, max_val))

    return bnds

#--------------------------------

def PrintResults1D(N_IS, Nexp, list_exp, best_IS, best_phi, chi2_red, std_error, chi2):
    print '\n\t Minimization results : \n'
    
    print '\t\t Interstellar flux parameters : '
    for i in range(0,N_IS):
            #print  '\t\t c%i = ' % i, pow(10,best_IS[i])   , '+/- ', pow(10,std_error[i]), 'diff = ', (pow(10,best_IS[i]) - pow(10,std_error[i]))
            print  '\t\t c%i = ' % i, best_IS[i]   , '+/- ', std_error[i]
    for i in range(0,Nexp):
        if (i<Nexp-1):
            print  '\t\t phi(', list_exp[i], ') = ', best_phi[i], '+/- ', std_error[N_IS + i], ', Chi2 =', chi2_red[i]
        elif (i == Nexp-1) : 
            print '\t\t phi(', list_exp[i], ') = ', best_phi[i], '+/- ', std_error[N_IS + i], ',Chi2 =', chi2_red[i], '\n'
    print "Chi2 global = ", chi2
    return 

#--------------------------------

def PrintResults1Db(N_IS, Nexp, list_exp, list_CR, list_exp_CR, best_IS, best_phi, chi2_red, std_error, chi2):
    print '\n\t Minimization results :'
    print '\t -------------------- \n'
    print '\t\t Interstellar flux parameters : '
    for k in range (0, len(list_CR)):
            print '\n\t\t\t ', list_CR[k], 'flux :'   
            for i in range(6*k,6*k+6):
                #print  '\t\t c%i = ' % i, pow(10,best_IS[i])   , '+/- ', pow(10,std_error[i]), 'diff = ', (pow(10,best_IS[i]) - pow(10,std_error[i]))
                print  '\t\t\t c%i = ' % i, best_IS[i]   , '+/- ', std_error[i]

    print "\n"
    for i in range(0,Nexp):
        if (i<Nexp-1):
            print  '\t\t phi(', list_exp[i], ') = ', best_phi[i], '+/- ', std_error[N_IS + i], ', Chi2 =', (chi2_red[i][j] for j in range(0,len(list_CR)))
        elif (i == Nexp-1) : 
            print '\t\t phi(', list_exp[i], ') = ', best_phi[i], '+/- ', std_error[N_IS + i], ',Chi2 =', (chi2_red[i][j] for j in range(0,len(list_CR))), '\n'
    print "Chi2 global = ", chi2
    return 
    
#--------------------------------

def ExtractDate(name_exp_list, list_CR):
 
    date_list_mean = []
    date_list_delta = []

    for name_exp in name_exp_list :

        st = []

        if "AMS" in name_exp and name_exp != "AMS02" and name_exp != "AMS01" :
            st.append('../Data/AMS/data_' + name_exp + '.dat')
        elif "PAMELA" in name_exp and not "PAMELA200" in name_exp and not "PAMELA0608" in name_exp :
            st.append('../Data/PAMELA/data_' + name_exp + '.dat')
        else :
            for i in range(0, len(list_CR)):
                st.append('../Data/' + list_CR[i] + '_data/data_' + name_exp + '.dat')

        Ncr = 0
        for i in range(0, len(list_CR)):

            if os.path.isfile(st[i])  :

                Ncr += 1
                date_beg, date_end = ReadDates(st[i])
                date_mean = date_beg + (date_end - date_beg)/2
                date_delta = (date_end - date_beg) / 2

                if Ncr == 1 :
                    date_mean_tmp = [date_mean] ; date_delta_tmp = [date_delta] ;

                else :
                    date_mean_tmp.append(date_mean) ; date_delta_tmp.append(date_delta) ;

        date_list_mean.append(date_mean_tmp)
        date_list_delta.append(date_delta_tmp)

    date_list_mean = np.array(date_list_mean)  
    date_list_delta = np.array(date_list_delta)
    
    return date_list_mean, date_list_delta
    
#--------------------------------

def ReadDates(file_name):

    extract_header  = open(file_name)
    header          = extract_header.readline()
    dates           = re.search('\((.*?)\)',header).group(1)
              
    if ('-' in header and re.search('\((.*?)\-',header)) :
        date_beg = re.search('\((.*?)\-',header).group(1)
        date_end = re.search(date_beg+'-(.*?)\)',header).group(1)
    else :
        date_beg = date_end = re.search('\((.*?)\)',header).group(1)
                
    date_beg = string.replace(date_beg,'/','-')
    date_end = string.replace(date_end,'/','-')

    date_beg = datetime.datetime.strptime(date_beg, "%Y-%m").date()
    date_end = datetime.datetime.strptime(date_end, "%Y-%m").date()

    return date_beg, date_end


#---------------------------------------------------------------

def ExtractExpDict(name_exp, list_CR): 
    # Extract data from one experiment file 

    st = []

    if "AMS" in name_exp and name_exp != "AMS02" and name_exp != "AMS01" :
        st.append('../Data/AMS/data_' + name_exp + '.dat')
    elif "PAMELA" in name_exp and not "PAMELA200" in name_exp and not "PAMELA0608" in name_exp :
        st.append('../Data/PAMELA/data_' + name_exp + '.dat')
    else :
        for i in range(0, len(list_CR)):
            st.append('../Data/' + list_CR[i] + '_data/data_' + name_exp + '.dat')

    Ncr = 0
    list_exp_CR =[]
    dict_Data = {}

    for i in range(0, len(list_CR)):

        if os.path.isfile(st[i])  :
            Ncr += 1
            list_exp_CR.append(list_CR[i])

            data_tmp      = np.loadtxt(st[i])
            dict_Data["Edata_%s_%s" % (list_CR[i],name_exp)] = data_tmp[:,0]
            dict_Data["ydata_%s_%s" % (list_CR[i],name_exp)] = data_tmp[:,3]
            dict_Data["sigma_%s_%s" % (list_CR[i],name_exp)] = data_tmp[:,9]
            Ndata_tmp = dict_Data["Edata_%s_%s" % (list_CR[i],name_exp)].size

            print dict_Data
            print Ndata_tmp

            if Ncr == 1 :
                Ndata = [Ndata_tmp]
            else :
                Ndata.append(Ndata_tmp)

    return Ndata, dict_Data, list_exp_CR
    