# -*- coding: utf-8 -*-

"""
Helper functions to get SSRGA parameter from normalized rime mass M follwing Maherndl et al.[2022, in prep.]
"""

from __future__ import division, print_function
import numpy as np

def ssrga_parameter(M):
    '''
    Calculate ssrga parameter kappa, beta, gamma, zeta1 and alpha_eff from the normalized rime mass M as defined in Seifert et al. (2019), 
    here: M = mass of rime / mass of size equivalent graupel assuming graupel density of 700kg/m³

    Input:
    M (unitless)

    Output
    kappa, beta, gamma, zeta1, alpha_eff
    '''

    p0 = 0.514
    p1 = np.array([0.16, -0.1, 4.06, -1.27, 0.127])
    p2 = np.array([0.187, 0.068, -7.45, 1.79, -0.091])
    p3 = np.array([0.575, 0.194, 5.42, 2.76, 0.067])

    alpha_eff = p1[0]*M**(2*p0)+p2[0]*M**p0+p3[0]
    kappa = 	p1[1]*M**(2*p0)+p2[1]*M**p0+p3[1]
    beta = 		p1[2]*M**(2*p0)+p2[2]*M**p0+p3[2]
    gamma = 	p1[3]*M**(2*p0)+p2[3]*M**p0+p3[3]
    zeta1 = 	p1[4]*M**(2*p0)+p2[4]*M**p0+p3[4]

    return kappa, beta, gamma, zeta1, alpha_eff

def scattering_name(M):
    try:
        para = ssrga_parameter(M)
        stringname = 'ss-rayleigh-gans_' + str(np.round(ssrga_parameter(M), 3)[0]) + '_' + str(np.round(ssrga_parameter(M), 3)[1]) + '_' + str(np.round(ssrga_parameter(M), 3)[2]) + '_' + str(np.round(ssrga_parameter(M), 3)[3]) 
        
            
        if len(stringname) > 40:
            print('enter single M value as float, not array!')
    except:
        print('enter single M value as float!')
        
    return stringname