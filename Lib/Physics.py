#!/usr/bin/python2.7
# -*-coding:Utf-8 -*

#
#  Conversions utiles
#

import numpy as np
from math import *

#--------------------------------
# Cosmic ray entry definition   
#--------------------------------

def CREntry(st):

	# CR properties [A, |Z|, mGeV] 
	if (st == "H"):
		cr = [1,1,mp_GeV()]
	elif (st == "He"):
		cr = [4,2,m4He_GeV()]
	elif (st == "Electron"):
		cr = [0,1,me_GeV()]
	return cr[0],cr[1],cr[2]

def GetA(cr):
	return cr[0]

def GetZ(cr):
	return cr[1]

def GetmGeV(cr):
	return cr[2]

#--------------------------------
# Conversion de distances (PDG)   
#--------------------------------

def Convert_cm_to_kpc():                
	return 3.2407793e-22

def Convert_m_to_kpc():                 
     	return 3.2407793e-20

def Convert_km_to_kpc():                
     	return 3.2407793e-17

def Convert_kpc_to_km():                
     	return 3.0856775e+16

def Convert_au_to_m():                  
	return 149597870700.


#--------------------------------
# Conversion de temps (PDG)   
#--------------------------------

# Sideral year = 1/31557600. s (from PDG 2005)

def Convert_s_to_Myr():             
	return 3.168808781e-14

def Convert_s_to_yr():             
	return 3.168808781e-8

def Convert_Myr_to_s():             
	return 3.15576e13

def Convert_h_to_Myr():             
	return 8.802246613e-18

def Convert_d_to_Myr():             
	return 3.667602755e-19

def Convert_yr_to_Myr():      
	return 1.e-6

def Convert_kyr_to_Myr():             
	return 1.e-3


#--------------------------------
# Conversion de vitesse
#--------------------------------

def Convert_kmpers_to_kpcperMyr():      
	return 1.022712e-3

def Convert_cmpers_to_kpcperMyr():      
	return 1.022712e-8

def Convert_mpers_to_kpcperMyr():       
	return 1.022712e-6

def Convert_AUperyr_to_mpers():         
	return 4.74047046293062522e+03

# Vitesse de la lumiere

def C_cms():
	return 2.99792458e10

def C_ms():
	return 2.99792458e8

def C_cmMyr():
	return 9.46073047e23


#--------------------------------
# Conversion coefficient de diffusion   
#--------------------------------

#Diffusion kpc/Myr -> cm2/s 

def Convert_kpc2perMyr_to_cm2pers():
	return 3.01715e29

def Convert_AU2peryr_to_cm2pers():      
	return 7.09164287370664960e+18


#--------------------------------
# Constantes   
#--------------------------------

#Electron mass (GeV)

def me_GeV():
	return 0.510998902e-3

#Proton mass (GeV)

def mp_GeV():
	return 0.938271998

#3He mass (GeV)

def m3He_GeV(): 
	return 2.8094132

#4He mass (GeV)

def m4He_GeV():
	return 3.728401


#-----------------------------------------------------------------   
#   Conversions entre les energies Etot, Ek, Ekn, R, p, beta... 
#   
#     - m is the mass (in amu or in GeV)
#     - E (total energy), Ek (kinetic energy) and p (momentum) are in GeV
#     - Ekn (kinetic energy per nucleon) is in GeV/n
#     - R (rigidity) is in GV
#-----------------------------------------------------------------   

def Beta_pE(p_gev, e_gev):

   #--- Return beta of the particle (= v/c = p/E).
   #  p                 Momentum of the particle [GeV]
   #  E                 Total energy of the particle [GeV]

   	return p_gev / e_gev 

#--------------------------------
def Beta_mEk(m_gev, ek_gev):

   #--- Returns beta of the particle (= v/c = p/E).
   #  m_gev             Mass of the particle [GeV]
   #  Ek                Kinetic energy (Ek=E+m) [GeV]

   	return sqrt(ek_gev * (ek_gev + 2. * m_gev)) / (ek_gev + m_gev) 


#--------------------------------
def Beta_mAEkn(m_gev, a_cr, ekn_gevn):

   #--- Returns beta of the particle (= v/c = p/E).
   #  m_gev             Mass of the particle [GeV]
   #  a_cr              Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)
   #  Ekn               Kinetic energy per nucleon (Ekn=Ek/n) of the particle [GeV/n]

   	if (a_cr == 0):
			a_cr = 1
   	ek_gev = ekn_gevn * float(a_cr)
   	return np.sqrt(ek_gev * (ek_gev + 2. * m_gev)) / (ek_gev + m_gev)

#--------------------------------
def Ekn_to_E(ekn_gevn, a_cr, m_gev):

   #--- Returns E from Ekn.
   #  E                 Total energy of the particle [GeV]
   #  Ekn               Kinetic energy per nucleon (Ekn=Ek/n) of the particle [GeV/n]
   #  A                 Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)
   #  m_gev             Mass of the particle [GeV]

   	if (a_cr == 0):
			a_cr = 1 
   	return ekn_gevn * float(a_cr) + m_gev 


#--------------------------------
def Ekn_to_Ek(ekn_gevn, a_cr):

   #--- Returns Ek from Ekn.
   #  Ek                Kinetic energy (Ek=E-m) of the particle [GeV]
   #  Ekn               Kinetic energy per nucleon (Ekn=Ek/n) of the particle [GeV/n]
   #  A                 Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)

   	if (a_cr == 0):
			a_cr = 1 
   	return ekn_gevn * float(a_cr) 


#--------------------------------
def Ekn_to_p(ekn_gevn, a_cr, m_gev):

   #--- Returns p from Ekn.
   #  p                 Momentum of the particle [GeV]
   #  Ekn               Kinetic energy per nucleon (Ekn=Ek/n) of the particle [GeV/n]
   #  A                 Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)
   #  m_gev             Mass of the particle [GeV]

   	if (a_cr == 0):	
			a_cr = 1 
   	ek_gev = ekn_gevn * float(a_cr) 
   	return sqrt(ek_gev * (ek_gev +  2. * m_gev)) 


#--------------------------------
def Ekn_to_R(ekn_gevn, a_cr, m_gev, z_cr):

   #--- Returns R from Ekn.
   #  R                 Rigidity (R=pc/Ze) of the particle [GV]
   #  Ekn               Kinetic energy per nucleon (Ekn=Ek/n) of the particle [GeV/n]
   #  A                 Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)
   #  Z                 Charge of the species

   	if (a_cr == 0):
			a_cr = 1 
   	ek_gev = ekn_gevn * float(a_cr) 
   	return np.sqrt(ek_gev * (ek_gev +  2. * m_gev)) / (float)(abs(z_cr)) 


#--------------------------------
def Ek_to_E(ek_gev, m_gev):

   #--- Returns E from Ek.
   #  E                 Total energy of the particle [GeV]
   #  Ek                Kinetic energy (Ek=E-m) of the particle [GeV]
   #  m_gev             Mass of the particle [GeV]

   	return ek_gev + m_gev 


#--------------------------------
def Ek_to_Ekn(ek_gev, a_cr):

   #--- Returns Ekn from Ek.
   #  Ekn               Kinetic energy per nucleon (Ekn=Ek/n) of the particle [GeV/n]
   #  Ek                Kinetic energy (Ek=E-m) of the particle [GeV]
   #  A                 Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)

   	if (a_cr == 0):
			a_cr = 1 
   	return ek_gev / float(a_cr) 

#--------------------------------
def Ek_to_p(ek_gev, m_gev):

   #--- Returns p from Ek.
   #  p                 Momentum of the particle [GeV]
   #  Ek                Kinetic energy (Ek=E-m) of the particle [GeV]
   #  m_gev             Mass of the particle [GeV]

	return sqrt(ek_gev * (ek_gev +  2. * m_gev)) 

#--------------------------------
def Ek_to_R(ek_gev, m_gev, z_cr):

   #--- Returns R from Ek.
   #  R  	       Rigidity (R=pc/Ze) of the particle [GV]
   #  Ek 	       Kinetic energy (Ek=E-m) of the particle [GeV]
   #  m_gev	       Mass of the particle [GeV]
   #  Z  	       Charge of the species

	return sqrt(ek_gev * (ek_gev +  2. * m_gev)) / float(abs(z_cr)) 

#--------------------------------
def E_to_Ek(e_gev, m_gev):

   #--- Returns Ek from E.
   #  E  	       Total energy of the particle [GeV]
   #  Ek 	       Kinetic energy (Ek=E-m) of the particle [GeV]
   #  m_gev	       Mass of the particle [GeV]

	return e_gev - m_gev 

#--------------------------------
def E_to_Ekn(e_gev, a_cr, m_gev):

   #--- Returns Ekn from E.
   #  E  	       Total energy of the particle [GeV]
   #  A  	       Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)
   #  Ekn	       Kinetic energy per nucleon (Ekn=Ek/n) of the particle [GeV/n]
   #  m_gev	       Mass of the particle [GeV]

	if (a_cr == 0):
		a_cr = 1 
	return (e_gev -  m_gev) / float(a_cr) 

#--------------------------------
def E_to_p(e_gev, m_gev):

   #--- Returns p from E.
   #  E  	       Total energy of the particle [GeV]
   #  p  	       Momentum of the particle [GeV]
   #  m_gev	       Mass of the particle [GeV]

	return sqrt(e_gev * e_gev - m_gev * m_gev) 

#--------------------------------
def E_to_R(e_gev, m_gev, z_cr):

   #--- Returns R from E.
   #  E  	       Total energy of the particle [GeV]
   #  m_gev	       Mass of the particle [GeV]
   #  Z  	       Charge of the species

	return sqrt(e_gev * e_gev - m_gev * m_gev) / float(abs(z_cr)) 

#--------------------------------
def Gamma_Em(e_gev, m_gev):

   #--- Returns gamma lorentz of the particle (= 1/(1-beta^2) = E/m).
   #  E  	       Total energy of the particle [GeV]
   #  m_gev	       Mass of the particle [GeV]

	return e_gev / m_gev 

#--------------------------------
def Gamma_mAZR(m_gev, a_cr, z_cr, r_gv):

   #--- Returns gamma lorentz of the particle (= 1/(1-beta^2) = E/m).
   #  m_gev	       Mass of the particle [GeV]
   #  A  	       Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)
   #  Z  	       Charge of the species
   #  R  	       Rigidity (R=pc/Ze) of the particle [GV]

   	if (a_cr == 0):
			a_cr = 1 
   	e = R_to_E(r_gv, m_gev, z_cr) 
   	return e / m_gev 

#--------------------------------
def Gamma_mAEkn(m_gev, a_cr, ekn_gevn):

   #--- Returns gamma lorentz of the particle (= 1/(1-beta^2) = E/m).
   #  m_gev	       Mass of the particle [GeV]
   #  A  	       Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)
   #  Ekn	       Kinetic energy per nucleon (Ekn=Ek/n) of the particle [GeV/n]
   #  E  	       Total energy of the particle [GeV]

   	if (a_cr == 0):
			a_cr = 1 
   	return (ekn_gevn * float(a_cr) + m_gev) / m_gev 

#--------------------------------
def p_to_Ek(p_gev, m_gev):

   #--- Returns Ek from p: Ek=E-m is the kinetic energy of the particle [GeV].
   #  p  	       Momentum of the particle [GeV]
   #  m_gev	       Mass of the particle [GeV]

   	return sqrt(p_gev * p_gev + m_gev * m_gev) - m_gev 

#--------------------------------
def p_to_Ekn(p_gev, a_cr, m_gev):

   #--- Returns Ekn from p.
   #  Ekn	       Kinetic energy per nucleon (Ekn=Ek/n) of the particle [GeV/n]
   #  p  	       Momentum of the particle [GeV]
   #  A  	       Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)
   #  m_gev	       Mass of the particle [GeV]

   	if (a_cr == 0):
			a_cr = 1 
   	return (sqrt(p_gev * p_gev + m_gev * m_gev) - m_gev) / float(a_cr) 

#--------------------------------
def p_to_E(p_gev, m_gev):

   #--- Returns E from p.
   #  E  	       Total energy of the particle [GeV]
   #  p  	       Momentum of the particle [GeV]
   # m_gev	       Mass of the particle [GeV]

   	return sqrt(p_gev * p_gev + m_gev * m_gev) 

#--------------------------------
def p_to_R(p_gev, z_cr):

   #--- Returns R from p.
   #  R  	       Rigidity (R=pc/Ze) of the particle [GV]
   #  p  	       Momentum of the particle [GeV]
   #  Z  	       Charge of the species

   	return p_gev / float(abs(z_cr)) 
	
#--------------------------------
def R_to_Ek(r_gv, m_gev, z_cr):

   #--- Returns Ek from R: Ek=E-m is the kinetic energy of the particle [GeV].
   #  R  	       Rigidity (R=pc/Ze) of the particle [GV]
   #  m_gev	       Mass of the particle [GeV]
   #  Z  	       Charge of the species

   	p_gev = r_gv * float(abs(z_cr)) 
   	return sqrt(p_gev * p_gev + m_gev * m_gev) - m_gev 

#--------------------------------
def R_to_Ekn(r_gv, a_cr, m_gev, z_cr):

   #--- Returns Ekn from R.
   #  Ekn	       Kinetic energy per nucleon (Ekn=Ek/n) of the particle [GeV/n]
   #  R  	       Rigidity (R=pc/Ze) of the particle [GV]
   #  A  	       Atomic number of the species (A=1 for leptons so that 'Ekn'(=Ek) exists)
   #  m_gev	       Mass of the particle [GeV]
   #  Z  	       Charge of the species

		if (a_cr == 0):
			a_cr = 1 
		p_gev = r_gv * float(abs(z_cr)) 
		return (np.sqrt(p_gev * p_gev + m_gev * m_gev) - m_gev) / float(a_cr) 

#--------------------------------
def R_to_E(r_gv, m_gev, z_cr):

   #--- Returns E from R.
   #  E  	       Total energy of the particle [GeV]
   #  R  	       Rigidity (R=pc/Ze) of the particle [GV]
   #  m_gev	       Mass of the particle [GeV]
   #  Z  	       Charge of the species

   	p_gev = r_gv * float(abs(z_cr)) 
   	return sqrt(p_gev * p_gev + m_gev * m_gev) 

#--------------------------------
def R_to_p(r_gv, z_cr):

   #--- Returns p from R.
   #  p  	       Momentum of the particle [GeV]
   #  R  	       Rigidity (R=pc/Ze) of the particle [GV]
   #  Z  	       Charge of the species

   	return r_gv * float(abs(z_cr)) 






