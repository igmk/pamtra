'''
connect the embedded 2moment scheme (from McSnow) to PAMTRA and run it 
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

# Initialize PyPAMTRA instance
pam = pyPamtra.pyPamtra()

#define some general values
n_heights = 50  #heights for PAMTRA calculation
debug = False #True->debug this code; False-> no special output

'''
set namelist parameter
'''
# Load and Modify Namelist file
pam.nmlSet
#turn off passive calculations
pam.nmlSet["passive"] = False  # Activate this for Microwave radiometer
pam.nmlSet["active"] = True    # Activate this for Cloud radar
pam.nmlSet["radar_mode"] = "spectrum"
pam.nmlSet["save_psd"] = False    # save particle size distribution
pam.nmlSet["radar_attenuation"] = "disabled" #"bottom-up"
#pam.nmlSet["hydro_adaptive_grid"] = True    # uncomment this line before changing max and min sp_diameter of size distribution
pam.nmlSet["conserve_mass_rescale_dsd"] = True    #necessary because of separating into m-D-relationship regions
#pam.nmlSet["hydro_fullspec"] = True #use full-spectra as input
#pam.nmlSet["radar_allow_negative_dD_dU"] = True #allow negative dU dD which can happen at the threshold between different particle species
#pam.nmlSet["radar_nfft"] = 4 #4096
#pam.nmlSet["radar_max_V"] = 3.
#pam.nmlSet["radar_min_V"] = -3.

#choose descriptor file
#desc_pamtra_file = 'descriptorfiles/descriptor_file_2m_crystal_ssrg.txt'
desc_pamtra_file = 'descriptorfiles/descriptor_file_2m_crystal_ssrg.txt'
#desc_pamtra_file = 'descriptorfiles/descriptor_file_2m_crystal-uli_modified.txt'
#desc_pamtra_file = 'descriptorfiles/descriptor_file_2m_icon.txt'
#load descriptor file
pam.df.readFile(desc_pamtra_file)

#show some messages
if debug:
    print 'debugging'
    pam.set["verbose"] = 8
else:
    pam.set["verbose"] = 0
pam.set["pyVerbose"] = 0

#deactivate particle types:
if "deact" in os.environ:
    deact = os.environ["deact"]
else:
    deact='111111' #set values to zero to deactivate particle types; order: c,i,r,s,g,h
print "using the following categories (c,i,r,s,g,h): ", deact

#directory of experiments
directory = "/home/mkarrer/Dokumente/McSnow/MCSNOW/experiments/"

#load file
twomom_file = Dataset(directory + experiment + '/twomom_d.ncdf',mode='r')
#create dictionary for all variables from PAMTRA
twomom = dict()
#print pam_file.variables.keys() 

#if necessary change name of variables
varlist = twomom_file.variables
#read PAMTRA variables to pamData dictionary
i_timestep=(tstep/30)-1 #there is no output for t=0min after that there are output steps in 30 minute steps (this could vary)
for var in varlist:#read files and write it with different names in Data
    if var in ("times","heights"):
        twomom[var] = np.squeeze(twomom_file.variables[var])
    else:
        twomom[var] = np.squeeze(twomom_file.variables[var][i_timestep,:])


#read atmospheric variables
atmo = __postprocess_McSnow.read_atmo(experiment)
'''
seperate by height
'''
#create height vector
model_top = 5000. #top of model / m
heightvec = np.linspace(model_top/n_heights,model_top,n_heights) #[::-1] #start with 0+z_res and go n_heigts step up to model_top
zres = heightvec[1]-heightvec[0]

#interpolate hydrometeor variables to heightvec    
twomom = __postprocess_McSnow.interpolate_2height(twomom,heightvec,twomom["heights"])
#interpolate atmospheric variables to heightvec
atmo_interpolated = __postprocess_McSnow.interpolate_atmo(atmo,heightvec)
# Generate PAMTRA data dictonary
pamData = dict()
vec_shape = [1,n_heights]
## Copy data to PAMTRA dictonary
#TODO: interpolate atmo if not a multiple of 5m should be used for vertical spacing
pamData["press"]  =  atmo_interpolated["p"] # Pressure Pa
pamData["relhum"] =  atmo_interpolated["rh"]  # Relative Humidity in  
pamData["timestamp"] =  0 #unixtime #TODO: not randomly set here
#pamData["lat"] = 0
#pamData["lon"] = 0
#TODO: set lfrac u,v,T_g
#pamData["lfrac"] = icon_frland[timeframe]
#pamData["wind10u"] = iconData["u"][0,1]
#pamData["wind10v"] = iconData["v"][0,1]
#pamData["groundtemp"] = iconData["t_g"][0]

pamData["hgt"] = atmo_interpolated["z"] #np.arange(0.,12000.,(12000.-0.)/(vec_shape[1])) #height in m 
pamData["temp"] = atmo_interpolated["T"] #T in K 


#write hydrometeors into one variable named hydro_cmpl for mass mixing ratios and hydro_num_cmpl for number
hydro_cmpl 	= 	np.zeros([1,1,n_heights,6])#allocate
hydro_num_cmpl 	= 	np.zeros([1,1,n_heights,6])#allocate

for i_hydro,hydromet in enumerate(["qc","qi","qr","qs","qg","qh"]):
    i_hydro
    if deact[i_hydro]=='1':
        hydro_cmpl[0,0,:,i_hydro]		=	twomom[hydromet]
        print hydromet,twomom[hydromet]
    elif deact[i_hydro]=='0':
        hydro_cmpl[0,0,:,i_hydro]		=	np.zeros(twomom[hydromet].shape)
for i_hydro,hydromet in enumerate(["qnc","qni","qnr","qns","qng","qnh"]):
    if deact[i_hydro]=='1':
        hydro_num_cmpl[0,0,:,i_hydro]		=	twomom[hydromet]
        print hydromet,twomom[hydromet]

    elif deact[i_hydro]=='0':
        hydro_num_cmpl[0,0,:,i_hydro]		=	np.zeros(twomom[hydromet].shape)

pamData["hydro_q"] = hydro_cmpl[:,:]
pamData["hydro_n"] = hydro_num_cmpl[:,:]
# Add them to pamtra object and create profile
pam.createProfile(**pamData)

#set hydrometeor properties
pam.df.data
# Execute PAMTRA for MIRA 35 GHz
pam.runParallelPamtra([9.6,35.5,95], pp_deltaX=1, pp_deltaY=1, pp_deltaF=1, pp_local_workers='auto') #pam.runPamtra([9.6,35.5,95])
# Write output to NetCDF4 file
if not deact=='111111': #add descriptor of used categories if not all categories are used
    testcase = testcase + '_cat' + deact

out_filename = "PAMTRA_2mom_" + testcase + '_av_' + str(av_tstep) + '_' + experiment + "_t" + str(tstep).zfill(4) + 'min'
pam.writeResultsToNetCDF("output/" + out_filename + ".nc")
subprocess.call(["cp","output/" + out_filename + ".nc",directory + experiment + '/' + out_filename + ".nc"])
