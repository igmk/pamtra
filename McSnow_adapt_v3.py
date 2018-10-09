'''
adaption of McSnow to PAMTRA
'''
#from IPython.core.debugger import Tracer ; Tracer()() #insert this line somewhere to debug

#import modules
import numpy as np
import sys #for some debugging
import pyPamtra
import re #for selecting numbers from file-string
import subprocess #for using shell commands with subprocess.call()
import os #import environment variables
from netCDF4 import Dataset

#self written #these files are linked here from the pythoncode/functions directory
import __postprocess_McSnow
import __general_utilities

#read variables passed by shell script
tstep = int(os.environ["tstep"])
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"] #"more readable" string of the experiment specifications
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output
MC_dir = os.environ["MC"]
# Initialize PyPAMTRA instance
pam = pyPamtra.pyPamtra() #load pyPamtra class (in core.py)

#define some general values
n_heights = 51  #heights for PAMTRA calculation
debug = False #True->debug this code; False-> no special output

'''
set namelist parameter
'''
#turn off passive calculations
pam.nmlSet["passive"] = False  # Activate this for Microwave radiometer
pam.nmlSet["active"] = True    # Activate this for Cloud radar
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet["save_psd"] = False    # save particle size distribution
pam.nmlSet["radar_attenuation"] = "disabled" #"bottom-up"
pam.nmlSet["hydro_fullspec"] = True #use full-spectra as input
pam.nmlSet["radar_allow_negative_dD_dU"] = True #allow negative dU dD which can happen at the threshold between different particle species
#pam.nmlSet["conserve_mass_rescale_dsd"] = False
#pam.nmlSet["radar_nfft"] = 8192 #1024
#pam.nmlSet["radar_max_V"] = 3.
#pam.nmlSet["radar_min_V"] = -3.

#show some messages
if debug:
    print 'debugging'
    pam.set["verbose"] = 8
else:
    pam.set["verbose"] = 0
pam.set["pyVerbose"] = 0

#deactivate particle types:
deact='1100' #set values to zero to deactivate particle types; order: small ice, unrimed aggregates, partially rimed, graupel

#directory of experiments
directory = MC_dir + "/experiments/"

#load file with superparticles (SP)
SP_file = Dataset(directory + experiment + '/mass2fr_' + str(tstep).zfill(4) + 'min_avtstep_' + str(av_tstep) + '.ncdf',mode='r')
#create dictionary for all variables from PAMTRA
SP = dict()
#print pam_file.variables.keys() #if needed show keys (variables) of SPs

#read PAMTRA variables to pamData dictionary
for var in SP_file.variables:#read files and write it with different names in Data
    SP[var] = np.squeeze(SP_file.variables[var])

#read atmospheric variables
filestring_atmo = directory + experiment + "/atmo.dat"
atmo = __postprocess_McSnow.read_atmo(experiment,filestring_atmo)

#create height vector
#get Sp with maximum height for upper limit
model_top = np.nanmax(SP['height'])
#model_top = 5000. #top of model / m #TODO: flexible input for model_top
heightvec_bound = np.linspace(0,model_top,n_heights) #start with 0+z_res and go n_heigts step up to model_top
heightvec_center = heightvec_bound[:-1] + 0.5 * np.diff(heightvec_bound) #center of the height-bins
zres = heightvec_bound[1]-heightvec_bound[0]
#interpolate atmospheric variables to heightvec
atmo_interpolated = __postprocess_McSnow.interpolate_atmo(atmo,heightvec_center)

#seperate SP dictionary by different height-bins of heightvec
for var in SP_file.variables:#read files and write it with different names in Data
    for i_height,height_now in enumerate(heightvec_bound[:-1]):
        condition_in_height = np.logical_and(heightvec_bound[i_height]<=SP["height"],SP["height"]<heightvec_bound[i_height+1])
        SP[var + "_heightbin_" + str(i_height)] = SP[var][condition_in_height]
        
    
#calculate volume of box
Vbox = __postprocess_McSnow.calculate_Vbox(experiment,zres)

# Generate PAMTRA data dictonary
pamData = dict()
#vec_shape = [1,n_heights-1]
## Copy data to PAMTRA dictonary
#TODO: interpolate atmo if not a multiple of 5m should be used for vertical spacing
pamData["press"]  =  atmo_interpolated["p"] # Pressure Pa
pamData["relhum"] =  atmo_interpolated["rh"]  # Relative Humidity in  
pamData["timestamp"] =  0 #unixtime
#pamData["lat"] = 0
#pamData["lon"] = 0
#TODO: set lfrac u,v,T_g
#pamData["lfrac"] = icon_frland[timeframe]
#pamData["wind10u"] = iconData["u"][0,1]
#pamData["wind10v"] = iconData["v"][0,1]
#pamData["groundtemp"] = iconData["t_g"][0]

pamData["hgt"] = atmo_interpolated["z"] #np.arange(0.,12000.,(12000.-0.)/(vec_shape[1])) #height in m 
pamData["temp"] = atmo_interpolated["T"] #T in K 

#determine number of SP to get number of necessary categories
number_ofSP = SP['m_tot'].shape[0]

#get necessary parameter of m-D and A-D relationship
mth,unr_alf,unr_bet,rhoi,rhol,Dth = __postprocess_McSnow.return_parameter_mD_AD_rel()
#selecting fallspeed model from testcase string
if "HW" in testcase:
    fallsp_model='heymsfield10_particles'
else:
    fallsp_model='khvorostyanov01_particles'

###
#handling the categories (TODO: until now just one)
###
nbins = 100
pam.df.addHydrometeor(("McSnowsmallice",  1.0,           -1,      -99,        -99,      -99,    -99,   -99,     13,              nbins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'mie-sphere',  fallsp_model,        0.)) #add hydrometeors: see
pam.df.addHydrometeor(("McSnowaggreg",  0.6,           -1,      -99,        -99,      -99,    -99,   -99,     13,              nbins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'ss-rayleigh-gans',  fallsp_model,        0.)) #add hydrometeors: see
N_cat = 2; #i_cat = 0 

#give some properties already here (f.e. scattering model can not be given on the run)
pamData["hydro_q"] = np.zeros([1,n_heights-1,N_cat]) #n_heights-1 because we are loosing one height due to binning
pamData["hydro_n"] = np.zeros([1,n_heights-1,N_cat]) 
# Add them to pamtra object and create profile
pam.createProfile(**pamData)
#set hydrometeor properties
pam.df.dataFullSpec
#initialize FullSpectra 
pam.nmlSet["hydro_fullspec"] = True
pam.df.addFullSpectra()

#initialize diameter arrays; dimensions are: (x,y,z,hydrometeor,bin)
for i_cat in range(0,N_cat):
    pam.df.dataFullSpec["d_bound_ds"][0,0,:,i_cat,:],dum =  np.meshgrid(10**np.linspace(-9,0,nbins+1),np.arange(0,n_heights-1))#2D-grid dimension:(height,bins); matrix with sp_diameters which is repeted N_height times
    pam.df.dataFullSpec["d_ds"][0,0,:,i_cat,:] = pam.df.dataFullSpec["d_bound_ds"][0,0,:,i_cat,:-1] + 0.5 * np.diff(pam.df.dataFullSpec["d_bound_ds"][0,0,:,i_cat,:])#center of the bins defined by d_bound_ds

#loop over all heights to perform KDE at each height
for idx_height,height_now in enumerate(heightvec_center):
    #constant for bandwidth from Shima 2009
    sigma0         = 0.62     #! Shima 2009, Sec. 5.1.4 
    #determine number of SP
    number_ofSP_in_heightbin = SP["m_tot_heightbin_" + str(idx_height)].shape[0] #; print "number_ofSP",number_ofSP
    # sigma for kernel estimate, sigma = sigma0/N_s**(1/5), see Shima Sec 5.1.4
    sigmai = (sigma0 / number_ofSP_in_heightbin**0.2) #ATTENTION:  this is defined **-1 in McSnow's mo_output.f90 
    
    #transform to logspace and set target array
    x = SP["diam_heightbin_" + str(idx_height)] #np.log(SP["diam_heightbin_" + str(idx_height)]) #logarithmate (base e)
    x_grid = pam.df.dataFullSpec["d_ds"][0,0,idx_height,0,:] #this must be equidistant (in log space) to calculate x_grid_logdiff!!
    x_grid_log = np.log(x_grid)
    x_grid_logdiff=x_grid_log[1]-x_grid_log[0]
    #x_grid[0]=-100
    
    #perform kde self-written routine (take multiplicity into account)
    pdf_fixed_bandwith_self_written_weighted = __postprocess_McSnow.kernel_estimate(x,x_grid,sigmai,weight=SP["xi_heightbin_" + str(idx_height)])
    pdf_fixed_bandwith_self_written_weighted= pdf_fixed_bandwith_self_written_weighted*x_grid_logdiff #normalize it so the integral is one (this is the real pdf)
    #TODO: convert pdf to [#/m-3] by multipling with absolute number of RP and dividing by the volume of the box
    n_ds = pdf_fixed_bandwith_self_written_weighted*sum(SP["xi_heightbin_" + str(idx_height)])/Vbox
    #print "n_ds",n_ds,"height",height_now
    ################################
    #'feed' PAMTRA with hydrometeors
    ################################
    
    #separate between small ice and unrimed aggregates
    #Dth=1e-3
    i_th = np.argmax(pam.df.dataFullSpec["d_ds"][0,0,idx_height,0,:]>Dth) #index corresponds to first diameter with exceeds Dth (threshold between spherical small ice and aggregates)
    
    #pass values for small ice
    #'''
    if not all(n_ds[:i_th])==0 or deact[0]=='1':#see switch deact at beginning #pass values for small spherical ice
        i_cat = 0
        #print "i_th",i_th,n_ds[:i_th-1],n_ds[i_th:]
        #print height_now,"n_ds",n_ds[:i_th]
        pam.df.dataFullSpec["n_ds"][0,0,idx_height,i_cat,:i_th] = n_ds[:i_th] #number density [#/m3] for small spherical ice
        pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i_cat,:i_th] = np.pi/6.*rhoi*pam.df.dataFullSpec["d_ds"][0,0,idx_height,0,:i_th]**3 #number density [#/m3] for small spherical ice
        pam.df.dataFullSpec["rho_ds"][0,0,idx_height,i_cat,:i_th] = rhoi #*pam.df.dataFullSpec["d_ds"][0,0,idx_height,0,i_th:]**3 #mass at the middle of the bin for small spherical ice
        pam.df.dataFullSpec["area_ds"][0,0,idx_height,i_cat,:i_th] = np.pi/4.*pam.df.dataFullSpec["d_ds"][0,0,idx_height,i_cat,:i_th]**2. #area at the middle of the bin for small spherical ice
        pam.df.dataFullSpec["as_ratio"][0,0,idx_height,i_cat,:i_th] = 1.0
        #'''
    if not all(n_ds[i_th:])==0 or deact[1]=='1':#see switch deact at beginning :
        #pass values for aggregates
        i_cat = 1
        #from IPython.core.debugger import Tracer ; Tracer()()
        pam.df.dataFullSpec["n_ds"][0,0,idx_height,i_cat,i_th:] = n_ds[i_th:] #number density [#/m3] for small spherical ice
        #pam.df.dataFullSpec["n_ds"][0,0,idx_height,i_cat,:] = n_ds[:] #number density [#/m3] for small spherical ice
        #mass at middle of bin
        b_agg=2.1;a_agg = 2.8*10**(2*b_agg-6)
        pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i_cat,:] = a_agg*pam.df.dataFullSpec["d_ds"][0,0,idx_height,i_cat,:]**b_agg #mass at the middle of the bin for aggregates
        #area of middle of bin
        d_agg=1.88;c_agg = 2.285*10**(2*d_agg-5)
        pam.df.dataFullSpec["area_ds"][0,0,idx_height,i_cat,:] = c_agg*pam.df.dataFullSpec["d_ds"][0,0,idx_height,i_cat,:]**d_agg #area at the middle of the bin for aggregates
        pam.df.dataFullSpec["as_ratio"][0,0,idx_height,i_cat,:] = 0.6

'''
import matplotlib.pyplot as plt
fig, axes = plt.subplots(nrows=4, ncols=1,)
axes[0].plot(pam.df.dataFullSpec["n_ds"][0,0,0,0,:])
axes[1].semilogy(pam.df.dataFullSpec["mass_ds"][0,0,0,0,:])
axes[2].semilogy(pam.df.dataFullSpec["area_ds"][0,0,0,0,:])
axes[3].plot(pam.df.dataFullSpec["as_ratio"][0,0,0,0,:])
plt.show()
'''
#from IPython.core.debugger import Tracer ; Tracer()()

#run PAMTRA
pam.runParallelPamtra([9.6,35.5,95], pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers='auto') #pam.runPamtra([9.6,35.5,95])
# Write output to NetCDF4 file
out_filename = "adaptv3_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
pam.writeResultsToNetCDF("output/" + out_filename + ".nc")
#copy netcdf to corresponding MCSNOW/experiment folder and show the user some output where it is
subprocess.call(["cp","output/" + out_filename + ".nc",directory + experiment + '/' + out_filename + ".nc"])
print "check results at: " + directory + experiment + "/" + out_filename + ".nc"
