# test the adaption of McSnow

import numpy as np
import sys
import pyPamtra


# Initialize PyPAMTRA instance
pam = pyPamtra.pyPamtra()
#define hydrometeor properties
n_bins = 100 #bins for PAMTRA calculation
n_heights = 100 #heigts for PAMTRA calculation
pam.df.addHydrometeor(('ice_nonspher',  -999,           -1,        -99.,      -99.,    -99.,   -99.,     -99.,   13,           n_bins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'ss-rayleigh-gans',  'dummy',        0.))


'''
set namelist parameter
'''
#turn off passive calculations
pam.nmlSet["passive"] = False  # Activate this for Microwave radiometer
pam.nmlSet["active"] = True    # Activate this for Cloud radar
pam.nmlSet["radar_mode"] = "simple"
pam.nmlSet["save_psd"] = False    # save particle size distribution
pam.nmlSet["radar_attenuation"] = "disabled" #"bottom-up"
pam.nmlSet["hydro_adaptive_grid"] = False    # uncomment this line before changing max and min diameter of size distribution
pam.nmlSet["conserve_mass_rescale_dsd"] = False    #necessary because of separating into m-D-relationship regions
pam.nmlSet["hydro_fullspec"] = True #use full-spectra as input
pam.nmlSet["radar_allow_negative_dD_dU"] = True #allow negative dU dD which can happen at the threshold between different particle species

#show some messages
debugging = False
if debugging:
    print 'debugging'
    pam.set["verbose"] = 5
else:
    pam.set["verbose"] = 0
pam.set["pyVerbose"] = 0



#load asci with all timesteps and all SPs
allSPalltimesteps = np.loadtxt("/home/mkarrer/Dokumente/McSnow/MCSNOW/experiments/1d_xi100000_nz5000_lwc20_ncl0_dtc5_nrp30_rm10_rt2_vt2_h10-20_ba500/mass2fr.dat")

#TODO: select timesteps
allSP = allSPalltimesteps
#read individual properties of the SPs
m_tot = allSP[:,0] #read mass of all particles
sp_height = allSP[:,2]
diam = allSP[:,9] #diameter
multipl = allSP[:,-1] #multiplicity


#seperate by height
model_top = 5000 #top of model / m
heightvec = np.linspace(0,5000,n_heights+1)



#define arrays with diam
d_bound_ds = np.logspace(-12,0,n_bins+1)
d_ds = d_bound_ds[:-1] + 0.5*np.diff(d_bound_ds)
d_counts = np.zeros([n_bins,n_heights])
#d_counts_no_mult = np.zeros(n_bins)
#get the number of particles (SP*multiplicity) at each bin (binned by height and diameter)
for i in range(0,n_heights):
    for j in range(0,n_bins):
                d_counts[j,i] = np.sum(np.where(np.logical_and(
                                np.logical_and(d_bound_ds[j]<=diam,diam<d_bound_ds[j+1]),
                                np.logical_and(heightvec[i]<=sp_height,sp_height<heightvec[i+1]),
                                ),multipl,0))
                #d_counts_no_mult[i] = np.sum(np.where(np.logical_and(d_bound_ds[i]<diam,diam<d_bound_ds[i+1]),1,0))


#output SP and RP counts per bin
print 'diameter       counts:  RP, SP'
for j in range(0,n_heights):
    print heightvec[j],d_counts[j,:]#,d_counts_no_mult[i]


'''
set up PAMTRA with some default values
'''
# Generate PAMTRA data dictonary
pamData = dict()
vec_shape = [1,n_heights]
## Copy data to PAMTRA dictonary
pamData["press"]  =  np.ones(vec_shape)*80000 # Pressure Pa #TODO: not randomly set here 800hPa
pamData["relhum"] =  np.ones(vec_shape)*80  # Relative Humidity in  #TODO: not randomly set here 80%
pamData["timestamp"] =  0 #unixtime #TODO: not randomly set here 80%
#pamData["lat"] = 0
#pamData["lon"] = 0
#TODO: set lfrac u,v,T_g
#pamData["lfrac"] = icon_frland[timeframe]
#pamData["wind10u"] = iconData["u"][0,1]
#pamData["wind10v"] = iconData["v"][0,1]
#pamData["groundtemp"] = iconData["t_g"][0]

pamData["hgt"] = np.ones(vec_shape)*heightvec[:-1] #np.arange(0.,12000.,(12000.-0.)/(vec_shape[1])) #height in m  #TODO: not randomly set here 80%
pamData["temp"] = np.ones(vec_shape)*263.15 #T in K   #TODO: not randomly set here -10C
pamData["hydro_q"] = np.zeros([1,n_heights,1]) #TODO: not just one category (last number)
pamData["hydro_n"] = np.zeros([1,n_heights,1]) #TODO: not just one category (last number)

# Add them to pamtra object and create profile
pam.createProfile(**pamData)
#set hydrometeor properties
pam.df.dataFullSpec

pam.nmlSet["hydro_fullspec"] = True
pam.df.addFullSpectra()


#print pam.df.dataFullSpec["d_bound_ds"].shape #dimensions are x,y,z,category,bins
#generate diameter arrays
pam.df.dataFullSpec["d_bound_ds"][0,0,:,0,:],dum =  np.meshgrid(d_bound_ds,np.arange(0,n_heights))#2D-grid dimension:(height,bins); matrix with diameters which is repeted N_height times
pam.df.dataFullSpec["d_ds"][0,0,:,0,:] = pam.df.dataFullSpec["d_bound_ds"][0,0,:,0,:-1] + 0.5 * np.diff(pam.df.dataFullSpec["d_bound_ds"][0,0,:,0,:])#center of the bins defined by d_bound_ds
################################
#'feed' PAMTRA with hydrometeors
################################
#number per bin
pam.df.dataFullSpec["n_ds"][0,0,:,0,:] = d_counts
#mass at middle of bin
b_agg=2.1;a_agg = 2.8*10**(2*b_agg-6)
pam.df.dataFullSpec["mass_ds"][0,0,:,0,:] = a_agg*pam.df.dataFullSpec["d_ds"][0,0,:,0,:]**b_agg #TODO: quick and dirty: m-D coefficients for aggregates
#area of middle of bin
d_agg=1.88;a_agg = 2.285*10**(2*d_agg-5)
pam.df.dataFullSpec["area_ds"][0,0,:,0,:] = pam.df.dataFullSpec["d_ds"][0,0,:,0,:]**b_agg #TODO: quick and dirty: m-D coefficients for aggregates
#aspect ratio
pam.df.dataFullSpec["as_ratio"][0,0,:,0,:] = 0.6


#run PAMTRA
pam.runPamtra(35.5)
# Write output to NetCDF4 file
pam.writeResultsToNetCDF('output/McSnow_test.nc')


