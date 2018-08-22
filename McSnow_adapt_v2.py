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
tstep = os.environ["tstep"]
experiment = os.environ["experiment"] #experiment name (this also contains a lot of information about the run)
testcase = os.environ["testcase"] #"more readable" string of the experiment specifications
av_tstep = int(os.environ["av_tstep"]) #average window for the McSnow output

# Initialize PyPAMTRA instance
pam = pyPamtra.pyPamtra()

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
#pam.nmlSet["hydro_adaptive_grid"] = True    # uncomment this line before changing max and min sp_diameter of size distribution
pam.nmlSet["conserve_mass_rescale_dsd"] = False    #necessary because of separating into m-D-relationship regions
pam.nmlSet["hydro_fullspec"] = True #use full-spectra as input
pam.nmlSet["radar_allow_negative_dD_dU"] = True #allow negative dU dD which can happen at the threshold between different particle species
#pam.nmlSet["radar_nfft"] = 4 #4096
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
directory = "/home/mkarrer/Dokumente/McSnow/MCSNOW/experiments/"

#load file
SP_file = Dataset(directory + experiment + '/mass2fr_' + tstep + '_avtstep_' + str(av_tstep) + '.ncdf',mode='r')
#create dictionary for all variables from PAMTRA
SP = dict()
#print pam_file.variables.keys() 

#if necessary change name of variables
varlist = SP_file.variables
#read PAMTRA variables to pamData dictionary
for var in varlist:#read files and write it with different names in Data
    SP[var] = np.squeeze(SP_file.variables[var])

#choose file (including timestep)
filestring_mass2fr = directory + experiment + "/mass2fr_" + tstep + ".dat"

#read mass2fr.dat file and get SP-dictionary
#SP = __postprocess_McSnow.read_mass2frdat(experiment,filestring_mass2fr)
atmo = __postprocess_McSnow.read_atmo(experiment)
'''
seperate by height
'''
#create height vector
model_top = 5000. #top of model / m
heightvec = np.linspace(model_top/n_heights,model_top,n_heights) #start with 0+z_res and go n_heigts step up to model_top
zres = heightvec[1]-heightvec[0]
#interpolate athmospheric variables to heightvec
atmo_interpolated = __postprocess_McSnow.interpolate_atmo(atmo,heightvec)

#calculate volume of box
Vbox = __postprocess_McSnow.calculate_Vbox(experiment,zres)

'''
set up PAMTRA with some arbitrary values
'''
# Generate PAMTRA data dictonary #TODO: input of test case and radiosoundes
pamData = dict()
vec_shape = [1,n_heights]
## Copy data to PAMTRA dictonary
#TODO: interpolate atmo if not a multiple of 5m should be used for vertical spacing
pamData["press"]  =  atmo_interpolated["p"] # Pressure Pa #TODO: not randomly; set here to 800hPa
pamData["relhum"] =  atmo_interpolated["rh"]  # Relative Humidity in  #TODO: not randomly; set here to 80%
pamData["timestamp"] =  0 #unixtime #TODO: not randomly set here
#pamData["lat"] = 0
#pamData["lon"] = 0
#TODO: set lfrac u,v,T_g
#pamData["lfrac"] = icon_frland[timeframe]
#pamData["wind10u"] = iconData["u"][0,1]
#pamData["wind10v"] = iconData["v"][0,1]
#pamData["groundtemp"] = iconData["t_g"][0]

pamData["hgt"] = atmo_interpolated["z"] #np.arange(0.,12000.,(12000.-0.)/(vec_shape[1])) #height in m  #TODO: not randomly set here 80%
pamData["temp"] = atmo_interpolated["T"] #T in K   #TODO: not randomly set here -10C

#determine number of SP to get number of necessary categories
number_ofSP = SP['m_tot'].shape[0]
#give some properties already here (f.e. scattering model can not be given on the run)
pamData["hydro_q"] = np.zeros([1,n_heights,number_ofSP]) #TODO: not just one category (last number)
pamData["hydro_n"] = np.zeros([1,n_heights,number_ofSP]) #TODO: not just one category (last number)

'''
handling the categories (as many as SP)
'''

#generate a list of string with enough entries (more than SP>>10^6)
if number_ofSP<238328:#count_len=3->2.3*10^5 diff. strings
    count_len = 3
elif 238328<number_ofSP<14776336: #count_len=4->1.4*10^7 diff. strings
    count_len = 4
elif 14776336<number_ofSP:
    print "check number of SP (exceeding anticipated number in this adaption)"
    sys.exit(0)
#get the string finally from a function
no_dupl_str = __general_utilities.gen_shortest_strings_without_duplicate(count_len=count_len)

#get necessary parameter of m-D and A-D relationship
mth,unr_alf,unr_bet,rhoi,rhol = __postprocess_McSnow.return_parameter_mD_AD_rel()
#print SP.keys()

#selecting fallspeed model from testcase string
if "HW" in testcase:
    fallsp_model='heymsfield10_particles'
else:
    fallsp_model='khvorostyanov01_particles'

nbins=2 #we give always just one SP per category, but for PAMTRA-internal reasons (calculating dD) we have to pass 2 bins
for i in range(0,number_ofSP): #number_ofSP):
    ##calculate m_i and save to SP dictionary
    SP["m_i"] = SP["m_tot"]-SP["m_rime"]-SP["m_wat"]
    SP["V_r"] = SP["m_rime"]/SP["rhor"]
    #define some epsilon for diameter bin
    diam_eps=1e-8
    #descriptor file structure
    #name 	as_ratio 	liq_ice 	rho_ms 	a_ms 	b_ms 	alpha 	beta 	moment_in 	nbin 	dist_name 	p_1 	p_2 	p_3 	p_4 	d_1 	d_2 	scat_name 	vel_size_mod 	canting
    #seperate between different m-D,A-D cases as in McSnows mo_mass2diam.f90 m2d_fillin()
    if SP["m_i"][i]<mth: #small ice
        if debug:
            print "smallice", SP["Frim"][i],SP["V_r"][i], SP["diam"][i],SP["proj_A"][i],SP["m_tot"][i]
        pam.df.addHydrometeor((no_dupl_str[i] + '_' + str(int(SP["Frim"][i]*10000)).zfill(5) + '_' + str(int(SP["vt"][i]*1000)).zfill(5),  1.0,           -1,      -99,        -99,      -99,    -99,   -99,     13,              nbins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'mie-sphere',  fallsp_model,        0.))
    else: #all other than small ice
        d_unrimed = (SP["m_i"][i]/unr_alf)**(1./unr_bet) #volume of the unrimed aggregate
        
        v_crit = np.pi/6*d_unrimed**3-SP["m_i"][i]*1./rhoi
        
        
        Fr = np.fmax(SP["V_r"][i],SP["m_rime"][i]*1./rhoi+SP["m_wat"][i]*1./rhol)/v_crit #this is not the same as SP["Frim"], rather it shows the relation between the actual rime volume and the critical rime volume, when the ice crystal is filled in
        if Fr==0.: #aggregates
            if debug:
                print "aggregates",Fr, SP["Frim"][i],SP["V_r"][i], SP["diam"][i],SP["proj_A"][i],SP["m_tot"][i]
            pam.df.addHydrometeor((no_dupl_str[i] + '_' + str(int(SP["Frim"][i]*10000)).zfill(5) + '_' + str(int(SP["vt"][i]*1000)).zfill(5),  0.6,           -1,      -99,        -99,      -99,     -99,   -99,     13,              nbins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'ss-rayleigh-gans',  fallsp_model,        0.))
        elif 0<Fr<1: #partially rimed (for now just use same properties as for aggregates)
            if debug:
                print "partially-rimed",Fr, SP["Frim"][i],SP["V_r"][i], SP["diam"][i],SP["proj_A"][i],SP["m_tot"][i]
            pam.df.addHydrometeor((no_dupl_str[i] + '_' + str(int(SP["Frim"][i]*10000)).zfill(5) + '_' + str(int(SP["vt"][i]*1000)).zfill(5),  -99.,           -1,      -99,        -99,      -99,     -99,   -99,     13,              nbins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'ss-rayleigh-gans',  fallsp_model,        0.))

        elif 1<=Fr: #graupel
            if debug:
                print "graupel",Fr, SP["Frim"][i],SP["V_r"][i], SP["diam"][i],SP["proj_A"][i],SP["m_tot"][i]
            pam.df.addHydrometeor((no_dupl_str[i] + '_' + str(int(SP["Frim"][i]*10000)).zfill(5) + '_' + str(int(SP["vt"][i]*1000)).zfill(5),  1.0,           -1,      -99,        -99,      -99,    -99.,   -99.,     13,              nbins,       'dummy',          -99.,    -99.,     -99.,   -99.,  -99.,    -99.        ,'mie-sphere',  fallsp_model,        0.))
        else:
            print "something wrong with the seperation of the categories?"
            print Fr, SP["Frim"][i],SP["V_r"][i], SP["diam"][i],np.maximum(SP["V_r"][i],SP["m_rime"][i]*1./rhoi+SP["m_wat"][i]*1./rhol),v_crit
            sys.exit(1)

# Add them to pamtra object and create profile
pam.createProfile(**pamData)
#set hydrometeor properties
pam.df.dataFullSpec


pam.nmlSet["hydro_fullspec"] = True
pam.df.addFullSpectra()


#print pam.df.dataFullSpec["d_bound_ds"].shape #dimensions are x,y,z,category,bins
##generate sp_diameter arrays
#d_bound_ds = SP["diam"]+diam_eps
d_bound_ds = np.logspace(-12,1,nbins+1) #array from 1nm to 1m
for i in range(0,number_ofSP): #number_ofSP):
    #select heightbin
    idx_height = np.argmax(SP["height"][i]<=heightvec)

    #dimensions are: (x,y,z,hydrometeor,bin)
    #PAMTRA needs two bins to calculate del_ds; only the first bin is really used
    pam.df.dataFullSpec["d_bound_ds"][0,0,idx_height,i,:] =  np.array([SP["diam"][i]-1.0*diam_eps,SP["diam"][i],SP["diam"][i]+1.0*diam_eps]) 
    pam.df.dataFullSpec["d_ds"][0,0,idx_height,i,:] = np.array([SP["diam"][i]-0.5*diam_eps,SP["diam"][i]+0.5*diam_eps])

    ################################
    #'feed' PAMTRA with hydrometeors
    ################################
    #number per bin
    pam.df.dataFullSpec["n_ds"][0,0,idx_height,i,:] = SP["xi"][i]/Vbox/2 #*diam_eps #/Vbox
    #mass at middle of bin
    #pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i,0:1] = SP["m_i"][i] #this is not used by mie-sphere
    #pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i,:] = (1.+SP["Frim"][i]/(1.-SP["Frim"][i]))*unr_alf*SP["diam"][i]**unr_bet #
    #area of middle of bin
    pam.df.dataFullSpec["area_ds"][0,0,idx_height,i,:] = SP["proj_A"][i]
    if SP["m_i"][i]<mth: #small ice
        if deact[0]=='1':#see switch deact at beginning
            pam.df.dataFullSpec["as_ratio"][0,0,idx_height,:,:] = 1.0 #":,:]" defines the aspect ratio also for other heights and the second diameter bin (allthough there are no particle), otherwise PAMTRA will throw an error: "nan or negative as_ratio"
            pam.df.dataFullSpec["rho_ds"][0,0,idx_height,i,:] = rhoi #6./np.pi/SP["diam"][i]**3*SP["m_tot"][i]
            pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i,:] = SP["m_tot"][i] 
        else:
            pam.df.dataFullSpec["n_ds"][0,0,idx_height,i,:] = 0
    else: #all other than small ice
        d_unrimed = (SP["m_i"][i]/unr_alf)**(1./unr_bet) #volume of the unrimed aggregate
        
        v_crit = np.pi/6*d_unrimed**3-SP["m_i"][i]*1./rhoi
        
        
        Fr = np.fmax(SP["V_r"][i],SP["m_rime"][i]*1./rhoi+SP["m_wat"][i]*1./rhol)/v_crit #this is not(?) the same as SP["Frim"]
        
        if Fr==0.: #aggregates
            if deact[1]=='1':#see switch deact at beginning
                pam.df.dataFullSpec["as_ratio"][0,0,idx_height,:,:] = 0.6
                pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i,:] = SP["m_i"][i] #this is not used by mie-sphere
            else:
                pam.df.dataFullSpec["n_ds"][0,0,idx_height,i,:] = 0
        elif 0<Fr<1: #partially rimed (for now just use same properties as for aggregates)
            if deact[2]=='1':#see switch deact at beginning
                pam.df.dataFullSpec["as_ratio"][0,0,idx_height,i,:] = 0.6                
                pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i,:] = SP["m_i"][i] #this is not used by mie-sphere

            else:
                pam.df.dataFullSpec["n_ds"][0,0,idx_height,i,:] = 0
        elif 1<=Fr: #graupel
            if deact[3]=='1':#see switch deact at beginning
                pam.df.dataFullSpec["as_ratio"][0,0,idx_height,i,:] = 1.0
                pam.df.dataFullSpec["rho_ds"][0,0,idx_height,i,:] = SP["m_tot"][i]/(SP["V_r"][i]+SP["m_i"][i]/rhoi) #TODO: must be rho=f(rhoi,rhor)
                pam.df.dataFullSpec["mass_ds"][0,0,idx_height,i,:] = pam.df.dataFullSpec["rho_ds"][0,0,idx_height,i,:]* pam.df.dataFullSpec["d_ds"][0,0,idx_height,i,:]**3
                #print SP["m_tot"][i]/(SP["V_r"][i]+SP["m_i"][i]/rhoi),SP["m_rime"][i]/SP["V_r"][i],SP["rhor"][i]
            else:
                pam.df.dataFullSpec["n_ds"][0,0,idx_height,i,:] = 0
        else:
            print "something wrong with the seperation of the categories?"
            print Fr, SP["Frim"][i],SP["V_r"][i], SP["diam"][i],np.maximum(SP["V_r"][i],SP["m_rime"][i]*1./rhoi+SP["m_wat"][i]*1./rhol),v_crit
            sys.exit(1)        
#run PAMTRA
pam.runPamtra([9.6,35.5,95])
# Write output to NetCDF4 file
out_filename = "adaptv2_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + tstep
pam.writeResultsToNetCDF("output/" + out_filename + ".nc")
subprocess.call(["cp","output/" + out_filename + ".nc",directory + experiment + '/' + out_filename + ".nc"])
print "check results at: " + directory + experiment + "/" + out_filename + ".nc"
