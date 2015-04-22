

Descriptor File Structure
==========================


This file contains microphysical specifications for all hydrometeors. It includes central parameters for PSD, m-D, A-D relation and selection of scattering model for frozen hydrometeors. All values (e.g. a and b parameter for m-D relation) have to be provided in SI units.

The complete list of options which can be used to define the PSD is given in table below. The table should be used as follows: given a distribution dist_name and a set of the parameters p_1 to p_4, moment_in describes the available choices for the moments of the PSD that have to be included in the input\_file. 1 means the total number concentration, 2 the effective radius, 3 the mass concentration, 13 both number and mass concentration and so on. d_1 and d_2 describe the diameter of the smallest and largest particle considered to calculate absorption and scattering properties for a given class of hydrometeors.

=============== =============== =============== =============== =============== ============== =============== =============== 
**dist\_name**  **p\_1**        **p\_2**        **p\_3**        **p\_4**        **moment\_in** **d\_1**        **d\_2**
=============== =============== =============== =============== =============== ============== =============== =============== 
mono            -99             -99             -99              -99            1 or 3         d\_1            -99
mono             N_T            -99             -99              -99            2 or 3         -99             -99
mono\_cosmo     -99             -99             -99              -99            3              -99             -99
exp             N_T             -99             -99              -99            2 or 3         d\_1            d\_2
exp             -99             R_eff           -99              -99            1 or 3         d\_1            d\_2
exp             -99             -99             N_0              -99            1, 2, or 3     d\_1            d\_2
exp             -99             -99             -99              -99            12, 13, or 23  d\_1            d\_2
exp\_field\_t   -99             -99             -99              -99            2 or 3         d\_1            d\_2
exp\_field\_tq  -99             -99             -99              -99            3              d\_1            d\_2
exp\_ryan       -99             -99             -99              -99            1 or 3         d\_1            d\_2
logn            N_T             -99             sigma            -99            2 or 3         d\_1            d\_2
logn            N_T             R_eff           -99              -99            3              d\_1            d\_2
logn            -99             R_eff           sigma            -99            1 or 3         d\_1            d\_2
logn            -99             -99             sigma            -99            12, 13, or 23  d\_1            d\_2
mgamma          N_T             -99             mu               gamma          2 or 3         d\_1            d\_2
mgamma          -99             R_eff           mu               gamma          1 or 3         d\_1            d\_2
mgamma          -99             -99             mu               gamma          12, 13, or 23  d\_1            d\_2
=============== =============== =============== =============== =============== ============== =============== =============== 

An exemplary descriptor\_file for a simulation with two hydrometeors is given in table \ref{tab:descriptor}. For the first hydrometeor category named cwc_q, the log-normal distribution have been used with a sigma of 0.38, fixed for the whole simulation via the parameter p\_3. Since the hydrometeor is in liquid phase, the density and the a and b parameters of the mass-size relation will be ignored. The aspect ratio and canting angle has been set to missing value (-99) because single scattering properties are calculated with Mie theory assuming spherical particles. The second hydrometeor swc_q is used to simulate snow. In this case an exponential distribution is used with a temperature dependent intercept parameter.
The thermodynamic state of the atmospheric columns included in the simulation and the surface properties are provided via the input\_file (IF), which can be in ascii or net-cdf format. The minimum set of parameters that has to be included for clear sky simulations are: surface temperature, longitude, latitude, date and profiles of temperature, pressure, height and relative humidity. In case of cloudy profiles, the moments of the PSD specified in the descriptor\_file need to be included. The PAMTRA model gives as output TB for 2 observational heights, one is fixed at ground and the other can be specified in the IF for each profile.

========= =============== ============== ============= =========== =========== =========== ========== ================ ========== ================ ========== ========== ========== ========== ============ ========== ================ ========================= ==================
**name**   **as\_ratio**   **liq\_ice**   **rho\_ms**   **a\_ms**   **b\_ms**   **alpha**   **beta**   **moment\_in**   **nbin**   **dist\_name**   **p\_1**   **p\_2**   **p\_3**   **p\_4**   **d\_1**     **d\_2**   **scat\_name**   **vel\_size\_mod**        **canting**
========= =============== ============== ============= =========== =========== =========== ========== ================ ========== ================ ========== ========== ========== ========== ============ ========== ================ ========================= ==================
cwc\_q     -99.            1              -99.          -99.        -99.        -99.        -99.       23               100       logn              -99.       -99.       0.38       -99.       1.e-12       1 .e-2      mie-sphere       khvorostyanov01\_drops    -99.
swc\_q     -99.            -1             -99.          0.038       2.0         0.3971      1.88       3                100        exp\_field\_t    -99.       -99.       -99.       -99.       0.51e-10     2 .e-2      mie-sphere       heymsfield10\_particles   -99.
========= =============== ============== ============= =========== =========== =========== ========== ================ ========== ================ ========== ========== ========== ========== ============ ========== ================ ========================= ==================


