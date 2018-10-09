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
n_heights = 50  #heights for PAMTRA calculation
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
pam.nmlSet["conserve_mass_rescale_dsd"] = False
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
deact='1111' #set values to zero to deactivate particle types; order: small ice, unrimed aggregates, partially rimed, graupel

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
heightvec = np.linspace(model_top/n_heights,model_top,n_heights) #start with 0+z_res and go n_heigts step up to model_top
zres = heightvec[1]-heightvec[0]
#interpolate atmospheric variables to heightvec
atmo_interpolated = __postprocess_McSnow.interpolate_atmo(atmo,heightvec)

#calculate volume of box
Vbox = __postprocess_McSnow.calculate_Vbox(experiment,zres)

# Generate PAMTRA data dictonary
pamData = dict()
vec_shape = [1,n_heights]
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

###
#handling the categories (as many as SP)
###
#generate a list of string with enough entries (more than SP>>10^6)
if number_ofSP<238328:#count_len=3->2.3*10^5 diff. strings
    count_len = 3
elif 238328<number_ofSP<14776336: #count_len=4->1.4*10^7 diff. strings
    count_len = 4
elif 14776336<number_ofSP:
    print "check number of SP (exceeding anticipated number in this adaption)"
    sys.exit(0)
#get the string (3-4 characters) which contains no duplicate from a function
no_dupl_str = __general_utilities.gen_shortest_strings_without_duplicate(count_len=count_len)

#get necessary parameter of m-D and A-D relationship
mth,unr_alf,unr_bet,rhoi,rhol,Dth = __postprocess_McSnow.return_parameter_mD_AD_rel()
#selecting fallspeed model from testcase string
if "HW" in testcase:
    fallsp_model='heymsfield10_particles'
else:
    fallsp_model='khvorostyanov01_particles'

nbins=2 #we give always just one SP per category, but for PAMTRA-internal reasons (calculating dD) we have to pass 2 bins
#Initialize array with information about particle type [0: not defined /or deactivated; 1: small ice; 2: aggregates; 3: partially rimed; 4: graupel]
particle_type = np.zeros(number_ofSP)
i_active = 0 #counts all active particles (some could have been deactivated by deact[i]=0)
for i in range(0,number_ofSP): #number_ofSP):
    ##calculate m_i and save to SP dictionary
    SP["m_i"] = SP["m_tot"]-SP["m_rime"]-SP["m_wat"]
    SP["V_r"] = SP["m_rime"]/SP["rhor"]
    #define some epsilon defining distance to second artificial diameter bin
    diam_eps=1e-8
    #descriptor file structure
    #name 	as_ratio 	liq_ice 	rho_ms 	a_ms 	b_ms 	alpha 	beta 	moment_in 	nbin 	dist_name 	p_1 	p_2 	p_3 	p_4 	d_1 	d_2 	scat_name 	vel_size_mod 	canting
    #seperate between different m-D,A-D cases as in McSnows mo_mass2diam.f90 m2d_fillin()
    if SP["m_i"][i]<mth: #small ice
        if debug:
            print "smallice", SP["Frim"][i],SP["V_r"][i], SP["diam"][i],SP["proj_A"][i],SP["m_tot"][i]
        if deact[0]=='1':#see switch deact at beginning
            pam.df.addHydrometeor((no_dupl_str[i_active] + '_' + str(int(SP["Frim"][i]*10000)).zfill(5) + '_' + str(int(SP["vt"][i]*1000)).zfill(5),  1.0,           -1,      -99,        -99,      -99,    -99,   -99,     13,              nbins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'mie-sphere',  fallsp_model,        0.)) #add hydrometeors: see descriptorFile.py to see what is done here
            particle_type[i]=1
        else:
            particle_type[i]=0
    else: #all other than small ice
        #calculate the maximum dimension of the unrimed aggregate
        d_unrimed = (SP["m_i"][i]/unr_alf)**(1./unr_bet) #volume of the unrimed aggregate
        #calculate critical volume (if there is more rime mass that fit in this volume Fr is > 1 and the particle is considered as infilled)
        v_crit = np.pi/6*d_unrimed**3-SP["m_i"][i]*1./rhoi
        Fr = np.fmax(SP["V_r"][i],SP["m_rime"][i]*1./rhoi+SP["m_wat"][i]*1./rhol)/v_crit #this is not the same as SP["Frim"], rather it shows the relation between the actual rime volume and the critical rime volume, when the ice crystal is filled in
        if Fr==0.: #aggregates
            if debug:
                print "aggregates",Fr, SP["Frim"][i],SP["V_r"][i], SP["diam"][i],SP["proj_A"][i],SP["m_tot"][i]
            if deact[1]=='1':#see switch deact at beginning
                pam.df.addHydrometeor((no_dupl_str[i_active] + '_' + str(int(SP["Frim"][i]*10000)).zfill(5) + '_' + str(int(SP["vt"][i]*1000)).zfill(5),  0.6,           -1,      -99,        -99,      -99,     -99,   -99,     13,              nbins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'ss-rayleigh-gans',  fallsp_model,        0.))
                particle_type[i]=2
            else:
                particle_type[i]=0
        elif 0<Fr<1: #partially rimed (for now just use same properties as for aggregates)
            if debug:
                print "partially-rimed",Fr, SP["Frim"][i],SP["V_r"][i], SP["diam"][i],SP["proj_A"][i],SP["m_tot"][i]
            if deact[2]=='1':#see switch deact at beginning
                pam.df.addHydrometeor((no_dupl_str[i_active] + '_' + str(int(SP["Frim"][i]*10000)).zfill(5) + '_' + str(int(SP["vt"][i]*1000)).zfill(5),  -99.,           -1,      -99,        -99,      -99,     -99,   -99,     13,              nbins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'ss-rayleigh-gans',  fallsp_model,        0.))
                particle_type[i]=3
            else:
                particle_type[i]=0            
        elif 1<=Fr: #graupel
            if debug:
                print "graupel",Fr, SP["Frim"][i],SP["V_r"][i], SP["diam"][i],SP["proj_A"][i],SP["m_tot"][i]
            if deact[2]=='1':#see switch deact at beginning
                pam.df.addHydrometeor((no_dupl_str[i_active] + '_' + str(int(SP["Frim"][i]*10000)).zfill(5) + '_' + str(int(SP["vt"][i]*1000)).zfill(5),  1.0,           -1,      -99,        -99,      -99,    -99.,   -99.,     13,              nbins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'mie-sphere',  fallsp_model,        0.))
                particle_type[i]=4
            else:
                particle_type[i]=0   
        else:
            print "something wrong with the seperation of the categories?"
            print Fr, SP["Frim"][i],SP["V_r"][i], SP["diam"][i],np.maximum(SP["V_r"][i],SP["m_rime"][i]*1./rhoi+SP["m_wat"][i]*1./rhol),v_crit
            sys.exit(1)
    if not particle_type[i]==0: #increase counter if particle type where the current particle belongs is not deactivated
        i_active+=1

#give some properties already here (f.e. scattering model can not be given on the run)
pamData["hydro_q"] = np.zeros([1,n_heights,i_active]) 
pamData["hydro_n"] = np.zeros([1,n_heights,i_active]) 
# Add them to pamtra object and create profile
pam.createProfile(**pamData)
#set hydrometeor properties
pam.df.dataFullSpec
#initialize FullSpectra 
pam.nmlSet["hydro_fullspec"] = True
pam.df.addFullSpectra()


#print pam.df.dataFullSpec["d_bound_ds"].shape #dimensions are x,y,z,category,bins
##generate sp_diameter arrays
i = 0
i_active = 0 #counts all active particles (some could have been deactivated by deact[i]=0)
while i<number_ofSP: #number_ofSP):
    if particle_type[i]==0:
        i+=1; continue
    #select heightbin
    idx_height = np.argmax(SP["height"][i]<=heightvec)

    #dimensions are: (x,y,z,hydrometeor,bin)
    #PAMTRA needs two bins to calculate del_ds; only the first bin is really used
    pam.df.dataFullSpec["d_bound_ds"][0,0,idx_height,i_active,:] =  np.array([SP["diam"][i]-diam_eps,SP["diam"][i]+diam_eps,SP["diam"][i]+2.0*diam_eps]) 
    pam.df.dataFullSpec["d_ds"][0,0,idx_height,i_active,:] = np.array([SP["diam"][i],SP["diam"][i]+1.5*diam_eps])

    ################################
    #'feed' PAMTRA with hydrometeors
    ################################
    #number per bin
    if SP["diam"][i]<1.0: #set upper threshold [in m]
        pam.df.dataFullSpec["n_ds"][0,0,idx_height,i_active,0] = SP["xi"][i]/Vbox #just the first of the two bins is filled; nevertheless mass_ds and area_ds is needed at both bins to calculate dU_dD in rescale_spectrum.f90 correctly
    else:
        pam.df.dataFullSpec["n_ds"][0,0,idx_height,i_active,0] = 0
    #area of middle of bin
    pam.df.dataFullSpec["area_ds"][0,0,idx_height,i_active,:] = SP["proj_A"][i]
    if particle_type[i]==1: #small ice
        pam.df.dataFullSpec["as_ratio"][0,0,idx_height,i_active,:] = 1.0 #":,:]" defines the aspect ratio also for other heights and the second diameter bin (allthough there are no particle), otherwise PAMTRA will throw an error: "nan or negative as_ratio"
        pam.df.dataFullSpec["rho_ds"][0,0,idx_height,i_active,:] = rhoi #6./np.pi/SP["diam"][i]**3*SP["m_tot"][i]
        pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i_active,:] = SP["m_tot"][i] #this is not used by mie-sphere, but needed by rescale_spectrum.f90
    elif particle_type[i]==2: #aggregates
        pam.df.dataFullSpec["as_ratio"][0,0,idx_height,i_active,:] = 0.6
        pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i_active,:] = SP["m_i"][i]
    elif particle_type[i]==3: #partially rimed (for now just use same properties as for aggregates)
        pam.df.dataFullSpec["as_ratio"][0,0,idx_height,i_active,:] = 0.6                
        pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i_active,:] = SP["m_i"][i]  #this is not used by mie-sphere, but needed by rescale_spectrum.f90
    elif particle_type[i]==4: #graupel
        pam.df.dataFullSpec["as_ratio"][0,0,idx_height,i_active,:] = 1.0
        pam.df.dataFullSpec["rho_ds"][0,0,idx_height,i_active,:] = SP["m_tot"][i]/(SP["V_r"][i]+SP["m_i"][i]/rhoi) #TODO: must be rho=f(rhoi,rhor)
        pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i_active,:] = pam.df.dataFullSpec["rho_ds"][0,0,idx_height,i_active,:]* pam.df.dataFullSpec["d_ds"][0,0,idx_height,i_active,:]**3
    if not particle_type[i]==0: #increase counter if particle type where the current particle belongs is not deactivated
        i+=1
        i_active+=1 #counts all active particles (some could have been deactivated by deact[i]=0)

#run PAMTRA
pam.runParallelPamtra([9.6,35.5,95], pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers='auto') #pam.runPamtra([9.6,35.5,95])
# Write output to NetCDF4 file
out_filename = "adaptv2_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
pam.writeResultsToNetCDF("output/" + out_filename + ".nc")
#copy netcdf to corresponding MCSNOW/experiment folder and show the user some output where it is
subprocess.call(["cp","output/" + out_filename + ".nc",directory + experiment + '/' + out_filename + ".nc"])
print "check results at: " + directory + experiment + "/" + out_filename + ".nc"
