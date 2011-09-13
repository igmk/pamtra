# -*- coding: utf-8 -*-
# met formeln fuer python
#
# ----------------------------------------------------------
# (c) Jan Schween 2005 (gnuplot) -> Mario Mech 2009 (python)
# converted to sSI units
# vaphet, pseudoAdiabLapseRate, rh2a, moist_rho_rh added by Max Maahn 2011
# ----------------------------------------------------------
#
#from math import *
#import sys
import numpy as np


Grav = 9.80665  # m/s^2 der Wert fuer mittlere Breiten
Rair = 287.04  # J/kg/K
Rvapor = 461.5  # J/kg/K
Cp = 1005.0	 # J/kg/K specific heat capacity
Gamma = -Grav/Cp  # =-0.0097..K/m  adiabatic temperature gradient
Lv = 2.5e6  # J/kg  bei 0C Lv heat of vaporization
Mwml = 0.622	# dimlos, Molmassenverhaeltnis
Tnull = -273.15 # K Absoluter Nullpunkt
Kadiab = Rair/Cp # dimlos  Adiabatenexponent

def moist_rho_rh(p,T,rh):
	'''
	Input:
	p is in Pa
	T is in K
	rh is in Pa/Pa
	
	Output:
	
	
	'''
	if np.any(rh > 1.5): raise TypeError("rh must not be in %")
	return p/(Rair*T_virt_rh(T,rh,p))

def moist_rho_q(p,T,q):
	'''
	Input p is in Pa
	T is in K
	q is in kg/kg
	
	'''
	return p/(Rair*T_virt_q(T,q))

def T_virt_rh(T,rh,p):
	'''
	Calculate the virtual temperature from air temperature,
	pressure, and relative humidity.

	Input:
	T is in K
	rh is in Pa/Pa
	p is in Pa
	
	Output:
	T_virt in K
	'''
	if np.any(rh > 1.5): raise TypeError("rh must not be in %")
	return T_virt_q(T,rh2q(rh,T,p))

def T_virt_q(T,q):
	  '''
	  Calculate the virtual temperature from air temperature and specific humidity.

	  Input:
	  T is in K
	  q is in kg/kg
	  
	  Output:
	  T_virt in K
	  '''
	  return T + T * 0.6078 * q


def rh2q(rh,T,p):
	'''
	Calculate the specific humidity from relative humidity, air temperature,
	and pressure. 

	Input:
	T is in K
	rh is in Pa/Pa
	p is in Pa
	
	Output
	q in kg/kg
	'''
	if np.any(rh > 1.5): raise TypeError("rh must not be in %")
	e = rh*e_sat_gg_water(T)
	q = Mwml*e/(p-(1-Mwml)*e)
	
	return q
	
	
def e_sat_gg_water(T):
	'''
	Calculates the saturation pressure over water after Goff and Gratch (1946).
	It is the most accurate that you can get for a temperture range from -90°C to +80°C.
	Source: Smithsonian Tables 1984, after Goff and Gratch 1946
	http://cires.colorado.edu/~voemel/vp.html
	http://hurri.kean.edu/~yoh/calculations/satvap/satvap.html

	Input:
	T in Kelvin.
	Output:
	e_sat_gg_water in Pa.
	'''
	return 100 * 1013.246 * 10**( -7.90298*(373.16/T-1) + 5.02808*np.log10(373.16/T) - 1.3816e-7*(10**(11.344*(1-T/373.16))-1) + 8.1328e-3 * (10**(-3.49149*(373.16/T-1))-1) )

def rh_to_iwv(relhum_lev,temp_lev,press_lev,hgt_lev):
	'''
	Calculate the integrated water vapour

	Input:
	T is in K
	rh is in Pa/Pa
	p is in Pa
	z is in m
	
	Output
	iwv in kg/m^2
	'''
	dz = np.diff(hgt_lev,axis=-1)
	relhum = (relhum_lev[...,0:-1] + relhum_lev[...,1:])/2.
	temp = (temp_lev[...,0:-1] + temp_lev[...,1:])/2.

	xp = -1.*np.log(press_lev[...,1:]/press_lev[...,0:-1])/dz
	press = -1.*press_lev[...,0:-1]/xp*(np.exp(-xp*dz)-1.)/dz

	q = meteoSI.rh2q(relhum,temp,press)
	rho_moist = meteoSI.moist_rho_q(press,temp,q)

	return np.sum(q*rho_moist*dz)


