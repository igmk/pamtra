'''
adaption of McSnow to PAMTRA
'''
#import modules
import numpy as np
import sys #for some debugging
import pyPamtra
import re #for selecting numbers from file-string
import subprocess #for using shell commands with subprocess.call()
import os #import environment variables

#self written
import __postprocess_McSnow

#read variables passed by shell script
tstep = os.environ["tstep"]
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)

# Initialize PyPAMTRA instance
pam = pyPamtra.pyPamtra()
#define hydrometeor properties
n_bins = 100 #bins for PAMTRA calculation
n_heights = 51 #heights for PAMTRA calculation
pam.df.addHydrometeor(('ice_nonspher',  -999,           -1,        -99.,      -99.,    -99.,   -99.,     -99.,   13,           n_bins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'ss-rayleigh-gans',  'khvorostyanov01_particles',        0.))

'''
set namelist parameter
'''
#turn off passive calculations
pam.nmlSet["passive"] = False  # Activate this for Microwave radiometer
pam.nmlSet["active"] = True    # Activate this for Cloud radar
pam.nmlSet["radar_mode"] = "simple"
pam.nmlSet["save_psd"] = False    # save particle size distribution
pam.nmlSet["radar_attenuation"] = "disabled" #"bottom-up"
pam.nmlSet["hydro_adaptive_grid"] = False    # uncomment this line before changing max and min sp_diameter of size distribution
pam.nmlSet["conserve_mass_rescale_dsd"] = False    #necessary because of separating into m-D-relationship regions
pam.nmlSet["hydro_fullspec"] = True #use full-spectra as input
#pam.nmlSet["radar_allow_negative_dD_dU"] = True #allow negative dU dD which can happen at the threshold between different particle species

#show some messages
debugging = False
if debugging:
    print 'debugging'
    pam.set["verbose"] = 8
else:
    pam.set["verbose"] = 0
pam.set["pyVerbose"] = 0

#directory of experiments
directory = "/home/mkarrer/Dokumente/McSnow/MCSNOW/experiments/"
#experiment name (this also contains a lot of information about the run
#experiment="1d_xi100000_nz5000_lwc20_ncl0_dtc5_nrp30_rm10_rt0_vt2_h10-20_ba500"
#choose file (including timestep)
filestring = directory + experiment + "/mass2fr_" + tstep + ".dat"

#read mass2fr.dat file and get SP-dictionary
SP = __postprocess_McSnow.read_mass2frdat(experiment,filestring)

from IPython.core.debugger import Tracer ; Tracer()()

'''
seperate by height
'''
#create height vector
model_top = 5000. #top of model / m
heightvec = np.linspace(0,model_top,n_heights) #start with 0+z_res and go n_heigts step up to model_top

#calculate values binned to here defined h-D bins
binned_val,heightvec,d_bound_ds,d_ds,zres = __postprocess_McSnow.seperate_by_height_and_diam(SP,nbins=100,diamrange=[-9,0],nheights=51,model_top=5000)

#calculate volume of box
Vbox = __postprocess_McSnow.calculate_Vbox(experiment,zres)
#divide by box volume to acchieve [#]->[#/m3]
binned_val["d_counts"] = __postprocess_McSnow.conv_num2numpm3(binned_val["d_counts"],Vbox)
binned_val["d_counts_no_mult"] = __postprocess_McSnow.conv_num2numpm3(binned_val["d_counts"],Vbox)


#output RP/m3 counts per bin
#print 'sp_diameter       counts:  RP/m3 SP/m3'
#print 'Volume of bin: V=' + str(Vbox)
#for j in range(0,n_heights-1):
#    print heightvec[j],binned_val["d_counts"][j,:]#,binned_val["d_counts_no_mult"]_no_mult[j,:]


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

pamData["hgt"] = np.ones(vec_shape)*heightvec #np.arange(0.,12000.,(12000.-0.)/(vec_shape[1])) #height in m  #TODO: not randomly set here 80%
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
#generate sp_diameter arrays
pam.df.dataFullSpec["d_bound_ds"][0,0,:,0,:],dum =  np.meshgrid(d_bound_ds,np.arange(0,n_heights))#2D-grid dimension:(height,bins); matrix with sp_diameters which is repeted N_height times
pam.df.dataFullSpec["d_ds"][0,0,:,0,:] = pam.df.dataFullSpec["d_bound_ds"][0,0,:,0,:-1] + 0.5 * np.diff(pam.df.dataFullSpec["d_bound_ds"][0,0,:,0,:])#center of the bins defined by d_bound_ds
################################
#'feed' PAMTRA with hydrometeors
################################
#number per bin
pam.df.dataFullSpec["n_ds"][0,0,:,0,:] = binned_val["d_counts"]
#mass at middle of bin
b_agg=2.1;a_agg = 2.8*10**(2*b_agg-6)
pam.df.dataFullSpec["mass_ds"][0,0,:,0,:] = a_agg*pam.df.dataFullSpec["d_ds"][0,0,:,0,:]**b_agg #TODO: quick and dirty: m-D coefficients for aggregates
#area of middle of bin
d_agg=1.88;c_agg = 2.285*10**(2*d_agg-5)
pam.df.dataFullSpec["area_ds"][0,0,:,0,:] = c_agg*pam.df.dataFullSpec["d_ds"][0,0,:,0,:]**d_agg #TODO: quick and dirty: m-D coefficients for aggregates
#aspect ratio
pam.df.dataFullSpec["as_ratio"][0,0,:,0,:] = 0.6


#run PAMTRA
pam.runPamtra([9.6,35.5,95])
# Write output to NetCDF4 file
pam.writeResultsToNetCDF("output/adaptv1_" + experiment + "_t" + tstep + ".nc")
subprocess.call(["cp","output/adaptv1_" + experiment + "_t" + tstep + ".nc",directory + experiment + "/" + "adaptv1" + "_t" + tstep + ".nc"])
print "check results at: " + directory + experiment + "/" + "adaptv1" + "_t" + tstep + ".nc"
