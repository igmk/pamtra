# -*- coding: utf-8 -*-
# met formeln fuer python
#
# ----------------------------------------------------------
# (c) Jan Schween 2005 (gnuplot) -> Mario Mech 2009 (python)
# converted to sSI units
# vaphet, pseudoAdiabLapseRate, rh2a, moist_rho2rh added by Max Maahn 2011
# 2026: most formulas here now delegate to the standalone meteo_si package
# (https://github.com/maahn/meteo_si), which carries the same formulas
# under independent tests/CI -- see AI.md for which functions still live
# here only (PAMTRA-specific adiabatic-LWC cloud physics) and why.
# ----------------------------------------------------------
#
from __future__ import division, print_function
import numpy as np
import meteo_si
try:
  from collections.abc import Iterable
except ImportError:
  from collections import Iterable

# Re-exported from meteo_si.constants as the single source of truth --
# `g` and `missingNumber` have no meteo_si.constants equivalent (the
# former is just a PAMTRA-local alias of Grav, the latter is a PAMTRA
# sentinel value, not a physical constant).
Grav = meteo_si.constants.Grav
Rair = meteo_si.constants.Rair
Rvapor = meteo_si.constants.Rvapor
Cp = meteo_si.constants.Cp
Gamma = meteo_si.constants.Gamma
Lv = meteo_si.constants.Lv
Mwml = meteo_si.constants.Mwml
Tnull = meteo_si.constants.Tnull
Kadiab = meteo_si.constants.Kadiab
g = Grav  # gravitational acceleration

missingNumber = -9999


# e_sat_gg_water is, despite the name, the Goff & Gratch (1946) formula --
# not to be confused with meteo_si.humidity.e_sat_gg_water, which (despite
# sharing that name) is a *different* formula (WMO CIMO Guide 2008). This
# alias keeps PAMTRA's own e_sat_gg_water name and Goff-Gratch behavior
# intact -- consistency with src/e_sat_gg_water.f90 (the Fortran core uses
# the same formula) is what matters here, not the meteo_si name it
# happens to share. See meteo_si's AI.md for the same caveat from its
# side.
e_sat_gg_water = meteo_si.humidity.e_sat_goffgratch_water


def moist_rho_rh(p, T, rh, *qm):
    '''
    Input:
    p is in Pa
    T is in K
    rh is in Pa/Pa
    Optional, several possible:
    qm is in kg/kg other species which contribute to the air mass! (ice, snow, cloud etc.)

    Output:
    density of moist air [kg/m^3]

    Example:
    moist_rho_rh(p,T,rh,q_ice,q_snow,q_rain,q_cloud,q_graupel,q_hail)


    '''
    if np.any(rh > 5):
        raise TypeError("rh must not be in %")

    q = rh2q(rh, T, p)

    return moist_rho_q(p, T, q, *qm)


def moist_rho_q(p, T, q, *qm):
    '''
    Input p is in Pa
    T is in K
    q is in kg/kg
    Optional, several possible:
    qm is in kg/kg other species which contribute to the air mass! (ice, snow, cloud etc.)

    Output:
    density of moist air [kg/m^3]

    Example:
    moist_rho_q(p,T,q,q_ice,q_snow,q_rain,q_cloud,q_graupel,q_hail)
    '''

    if len(qm) > 0:
        # get rid of masked data!
        qm = np.ma.array(qm).filled(0)
        qm[qm < 0] = 0
        qm = np.sum(qm, axis=0)
    else:
        qm = 0

    moist_rho_q = meteo_si.density.moist_rho_q(p, T, q, qm)

    if np.any(moist_rho_q < 0):
        if np.any(moist_rho_q < -0.001):
            raise ValueError(
                "meteoSI.moist_rho_q calculated negative densities!")
        else:
            try:
                moist_rho_q[moist_rho_q < 0] = 0
            except:
                moist_rho_q = 0

    return moist_rho_q


def T_virt_rh(T, rh, p):
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
    if np.any(rh > 5):
        raise TypeError("rh must not be in %")
    return meteo_si.temperature.T_virt_rh(T, rh, p)


def T_virt_q(T, q):
    '''
    Calculate the virtual temperature from air temperature and specific humidity.

    Input:
    T is in K
    q is in kg/kg

    Output:
    T_virt in K
    '''
    return meteo_si.temperature.T_virt_q(T, q)


def e2q(e, p):
    '''
    Calculate the specific humidity from water vapour pressure and air pressure.

    Input:
    e is in Pa
    p is in Pa

    Output
    q in kg/kg
    '''
    return meteo_si.humidity.e2q(e, p)


def q2e(q, p):
    '''
    Calculate water vapour pressure from the specific humidity and air pressure.

    Input:
    q in kg/kg
    p is in Pa

    Output
    e is in Pa
    '''
    return meteo_si.humidity.q2e(q, p)


def rh2q(rh, T, p):
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
    return meteo_si.humidity.rh2q(rh, T, p, e_sat_func=e_sat_gg_water)


def rh2a(rh, T):
    '''
    Calculate the absolute humidity from relative humidity, air temperature,
    and pressure.

    Input T is in K
    rh is in Pa/Pa
    p is in Pa
    Output
    a in kg/m^3

    Source: Kraus: Chapter 8.1.2
    '''
    return meteo_si.humidity.rh2a(rh, T, e_sat_func=e_sat_gg_water)


def a2rh(a, T):
    '''
    Calculate the relative from absolute humidity and air temperature.

    Input
    T is in K
    a in kg/kg
    Output
    rh in Pa/Pa

    Source: Kraus: Chapter 8.1.2
    '''
    return meteo_si.humidity.a2rh(a, T, e_sat_func=e_sat_gg_water)


def q2rh(q, T, p):
    '''
    Calculate relative humidity from specific humidity

    Input:
    T is in K
    p is in Pa
    q in kg/kg

    Output:
    rh is in Pa/Pa
    '''
    return meteo_si.humidity.q2rh(q, T, p, e_sat_func=e_sat_gg_water)


def rh_to_iwv(relhum_lev, temp_lev, press_lev, hgt_lev):
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
    return meteo_si.humidity.rh_to_iwv(
        relhum_lev, temp_lev, press_lev, hgt_lev, e_sat_func=e_sat_gg_water)


def detect_liq_cloud(z, t, rh):  # , rh_thres, t_thres):

    # UL NOV 2007
    # tranlated to python by mx 2011
    # ***********
    # INPUT
    # z: height grid
    # T: temperature on z
    # rh: relative humidty on z
    # rh_thres: relative humidity threshold for the detection on liquid clouds on z
    # T_thres: don not detect liquid water clouds below this value (scalar)
    # ***********
    # OUTPUT
    # z_top: array of cloud tops
    # z_base: array of cloud bases
    # z_cloud: array of cloudy height levels
    # ***********

    rh_thres = 0.95  # 1
    t_thres = 253.15  # K
    #import pdb; pdb.set_trace()
    n = len(z)
    # print "!",n
    # ***determine cloud boundaries
    # --> layers where mean rh GT rh_thres

    cloud_bound_ind = np.zeros(n, dtype=int)
    for i in np.arange(0, (n - 1)):
        #print ((rh[i + 1] + rh[i]) / 2. > rh_thres)
        #print ((t[i + 1] + t[i]) / 2. > t_thres)
        if ((rh[i + 1] + rh[i]) / 2. > rh_thres) and ((t[i + 1] + t[i]) / 2. > t_thres):
            cloud_bound_ind[i] = np.bitwise_or(1, cloud_bound_ind[i])
            cloud_bound_ind[i + 1] = np.bitwise_or(2, cloud_bound_ind[i + 1])
        # end if
    # end for
    #import pdb; pdb.set_trace()
    i_cloud = np.where(cloud_bound_ind != 0)[0]

    # ***determine z_base & z_top arrays
    #z_top = -99.
    #z_base = -99.
    #z_cloud = -99

    i_top = []
    i_base = []

    if len(i_cloud) != 0:
        #z_cloud = z[i_cloud]
        i_base = np.where(cloud_bound_ind == 1)[0]
        i_top = np.where(cloud_bound_ind == 2)[0]
        n_base = len(i_base)
        n_top = len(i_top)
        if n_top != n_base:
            print('something wrong, number of bases NE number of cloud tops!')
            return [], [], []
        # end if
    # end if
    #z_top = z[i_top]
    #z_base = z[i_base]

    return i_top, i_base, i_cloud
# end def detect_liq_cloud


def adiab(i, T, P, z):
    """
    Adiabtic liquid water content assuming pseudoadiabatic lapse rate
    throughout the whole cloud layer. Thus the assumed temperature
    profile is differnt from the measured one
    Input:
    i no of levels
    T is in K
    p is in Pa
    z is in m
    Output:
    LWC in
    translated to Python from adiab.pro by mx.
    """

    #   Set actual cloud base temperature to the measured one
    #   Initialize Liquid Water Content (LWC)
    #   Compute adiabatic LWC by integration from cloud base to level I

    TCL = T[0]
    LWC = 0.0

    for j in range(1, i+1):

        deltaz = z[j]-z[j-1]

        #   Compute actual cloud temperature

        #   1. Compute air density
        #   2. Compute water vapor density of saturated air
        #   3. Compute mixing ratio of saturated air
        #   4. Compute pseudoadiabatic lapse rate
        #   5. Compute actual cloud temperature

        R = moist_rho_rh(P[j], T[j], 1.)
        RWV = rh2a(1., T[j])
        WS = RWV/(R-RWV)
        DTPS = pseudoAdiabLapseRate(T[j], WS)
        TCL = TCL+DTPS*(deltaz)

        #   Compute adiabatic LWC

        #   1. Compute air density
        #   2. Compute water vapor density of saturated air
        #   3. Compute mixing ratio of saturated air
        #   4. Compute specific heat of vaporisation
        #   5. Compute adiabatic LWC

        R = moist_rho_rh(P[j], TCL, 1.)
        RWV = rh2a(1., TCL)
        WS = RWV/(R-RWV)
        L = vaphet(TCL)

        LWC = LWC+(R*Cp/L*((Grav/Cp)-pseudoAdiabLapseRate(TCL, WS))*(deltaz))
    # end for

    return LWC
# end def adiab


def mod_ad(T_cloud, p_cloud, z_cloud, fak):

    # ;IN: T_cloud, p_cloud, z_cloud
    # ;OUT: lwc
    # ;Einheiten: SI!
    # translated to Python from mod_ad.pro by mx.

    n_level = len(T_cloud)
    lwc = np.zeros(n_level-1)

    thick = 0.

    for jj in range(n_level-1):
        deltaz = z_cloud[jj+1] - z_cloud[jj]
        thick = deltaz + thick
        lwc[jj] = adiab(jj+1, T_cloud, p_cloud, z_cloud)
        lwc[jj] = lwc[jj]*(-0.144779*np.log(thick/fak) + 1.239387)
    # end for
    return lwc
# end def mod_ad


def pseudoAdiabLapseRate(T, Ws):
    """
    Pseudoadiabatic lapse rate
    Input: T   [K]  thermodynamic temperature
    Ws   [1]  mixing ratio of saturation
    Output: PSEUDO [K/m] pseudoadiabatic lapse rate
    Constants: Grav   [m/s2]     : constant of acceleration
        CP  [J/(kg K)]    : specific heat cap. at const. press
        Rair  [J/(kg K)]  : gas constant of dry air
        Rvapor [J/(kg K)] : gas constant of water vapor
    Source: Rogers and Yau, 1989: A Short Course in Cloud Physics
    (III.Ed.), Pergamon Press, 293p. Page 32
    translated to Python from pseudo1.pro by mx
    """

    # Compute specific humidity of vaporisation
    L = vaphet(T)

    # Compute pseudo-adiabatic temperature gradient
    x = (Grav/Cp) * (1+(L*Ws/Rair/T)) / (1+(Ws*L**2/Cp/Rvapor/T**2))

    return x


def vaphet(T):
    """
    Compute specific heat of vaporisation
    Input  : T      [K]      thermodynamic temperature
    Output : VAPHET [J/kg]   specific heat of vaporisation
    Source: Liljequist, G.H. und K. Cehak, 1984: Allgemeine
        Meteorologie (III.Ed.). Vieweg, 396p. Page 95
    translated to Python from vaphet.pro by mx
    """

    x = (2500.8-2.372*(T-273.15)) * 1000.0

    return x


usStandard = meteo_si.atmosphere.usStandard
