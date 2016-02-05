# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import csv
import cPickle as pickle
import time,calendar,datetime
import sys
import os
import random
import string
import itertools
from copy import deepcopy
from matplotlib import mlab
import multiprocessing
import logging
import glob
import warnings

from .libWrapper import PamtraFortranWrapper, parallelPamtraFortranWrapper
from .descriptorFile import pamDescriptorFile
from .tools import sftp2Cluster, formatExceptionInfo
from .meteoSI import detect_liq_cloud, mod_ad, moist_rho_rh,rh2q
from .fortranNamelist import Namelist

missingNumber=-9999.
missingIntNumber=int(missingNumber)

#logging.basicConfig(filename='example.log',level=logging.DEBUG)


    
class pyPamtra(object):

  '''
  Class for pamtra calculations. Initalisations fill dictonary 'set' with default values and 
  also fills 'dimensions' and 'units'
  
  '''
   
  
  def __init__(self):
    """
    set setting default values
    """
    
    self.default_p_vars = ["timestamp","lat","lon","lfrac","wind10u","wind10v","hgt","press","temp","relhum","hgt_lev","press_lev","temp_lev","relhum_lev","q","hydro_q","hydro_n","hydro_reff","wind10u","wind10v","obs_height", "ngridy","ngridx","max_nlyrs","nlyrs","model_i","model_j","unixtime","airturb","radar_prop","groundtemp","wind_w"]
    self.nmlSet = dict() #:settings which are required for the nml file. keeping the order is important for fortran
    #keys MUST be lowercase for f2py!
    self.nmlSet["hydro_threshold"]=  1.e-10   # [kg/kg] 
    #set namelist defaults#
    self.nmlSet["write_nc"]= True
    self.nmlSet["data_path"]= 'data/'
    self.nmlSet["save_psd"]= False #also saves the PSDs used for radiative transfer
    self.nmlSet["save_ssp"]= False #also saves the single scattering properties used for radiative transfer
    #self.nmlSet["noutlevels"]= 2 # output levels are given from top to bottom in meters
    self.nmlSet["add_obs_height_to_layer"] = False # add observation levels to atmospheric levels by interpolation 
    self.nmlSet["outpol"]= 'VH'
    self.nmlSet["file_desc"]= ''
    self.nmlSet["creator"]= 'Pamtrauser'
    self.nmlSet["active"]= True
    self.nmlSet["passive"]= True
    self.nmlSet["radar_mode"]= "simple" #"splitted"|"moments"|"spectrum"
    self.nmlSet["randomseed"] = 0 #0 is real noise, other value gives always the same random numbers
    self.nmlSet["ground_type"]= 'L'
    self.nmlSet["salinity"]= 33.0
    self.nmlSet["emissivity"]= 0.6
    self.nmlSet["lgas_extinction"]= True
    self.nmlSet["gas_mod"]= 'R98'
    self.nmlSet["lhyd_absorption"]= True
    self.nmlSet["lhyd_scattering"]= True
    self.nmlSet["lhyd_emission"]= True
    self.nmlSet["liq_mod"]= "Ell"
    self.nmlSet["hydro_includehydroinrhoair"] = True
    self.nmlSet["hydro_fullspec"] = False
    self.nmlSet["hydro_limit_density_area"] = True
    self.nmlSet["hydro_softsphere_min_density"] = 10.    
    self.nmlSet["hydro_adaptive_grid"] = True
    self.nmlSet["tmatrix_db"] = "none" 
    self.nmlSet["tmatrix_db_path"] = "database/" 
    #: number of FFT points in the Doppler spectrum [typically 256 or 512]
    self.nmlSet["radar_nfft"]= 256
    #: number of average spectra for noise variance reduction, typical range [1 150]
    self.nmlSet["radar_no_ave"]= 150
    #: MinimumNyquistVelocity in m/sec
    self.nmlSet["radar_max_v"]= 7.885
    #: MaximumNyquistVelocity in m/sec
    self.nmlSet["radar_min_v"]= -7.885
    #: radar noise in 1km in same unit as Ze 10*log10(mm⁶/m³). noise is calculated with noise"]=  radar_pnoise0 + 20*log10(range/1000)
    self.nmlSet["radar_pnoise0"]= -32.23 # mean value for BArrow MMCR during iSDAC
    self.nmlSet["radar_airmotion"]=  False
    self.nmlSet["radar_airmotion_model"]=  "step" #: "constant","linear","step"
    self.nmlSet["radar_airmotion_vmin"]=  -4.e0
    self.nmlSet["radar_airmotion_vmax"]=  +4.e0
    self.nmlSet["radar_airmotion_linear_steps"]=  30
    self.nmlSet["radar_airmotion_step_vmin"]=  0.5e0

    self.nmlSet["radar_aliasing_nyquist_interv"]=  1
    self.nmlSet["radar_save_noise_corrected_spectra"]=  False
    self.nmlSet["radar_use_hildebrand"]=  False
    self.nmlSet["radar_min_spectral_snr"]=  1.2#threshold for peak detection. if radar_no_Ave >> 150, it can be set to 1.1
    self.nmlSet["radar_convolution_fft"]=  True #use fft for convolution of spectrum. is alomst 10 times faster, but can introduce aretfacts for radars with *extremely* low noise levels or if noise is turned off at all.
    self.nmlSet["radar_smooth_spectrum"]=  True #smooth spectrum before moment estimation
    self.nmlSet["radar_k2"]=  0.93 # dielectric constant |K|² (always for liquid water by convention) for the radar equation
    self.nmlSet["radar_npeaks"] = 1
    self.nmlSet["radar_noise_distance_factor"]=  2.0
    self.nmlSet["radar_receiver_uncertainty_std"]=  0.e0 #dB
    self.nmlSet["radar_receiver_miscalibration"]=  0.e0 #dB
    self.nmlSet["radar_attenuation"]=  "disabled" #! "bottom-up" or "top-down"
    self.nmlSet["radar_polarisation"]=  "NN" #! comma separated
    self.nmlSet["radar_use_wider_peak"]=  False #
    self.nmlSet["liblapack"]=  True # use liblapack for matrix inversion
    
    
    
    #all settings which do not go into the nml file go here:
    self.set = dict()
    self.set["pyVerbose"] = 0
    self.set["verbose"] = 0
    self.set["freqs"] = []
    #self.set["nfreqs"] = 0
    self.set["namelist_file"] = "TMPFILE"

    self._nmlDefaultSet = deepcopy(self.nmlSet)
        
    self.dimensions = dict()
    
    self.dimensions["unixtime"] = ["ngridx","ngridy"]
    
    self.dimensions["ngridx"] = []
    self.dimensions["ngridy"] = []
    self.dimensions["nlyrs"] = ["ngridx","ngridy"]
    
    self.dimensions["model_i"] = ["ngridx","ngridy"]
    self.dimensions["model_j"] = ["ngridx","ngridy"]
    self.dimensions["lat"] = ["ngridx","ngridy"]
    self.dimensions["lon"] = ["ngridx","ngridy"]
    self.dimensions["lfrac"] = ["ngridx","ngridy"]
    self.dimensions["wind10u"] = ["ngridx","ngridy"]
    self.dimensions["wind10v"] = ["ngridx","ngridy"]
    self.dimensions["wind_w"] = ["ngridx","ngridy","max_nlyrs"]
    self.dimensions["obs_height"] = ["ngridx","ngridy","noutlevels"]

    self.dimensions["iwv"] = ["ngridx","ngridy"]

    self.dimensions["hgt_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
    self.dimensions["temp_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
    self.dimensions["p_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
    self.dimensions["relhum_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
    
    self.dimensions["hydro_q"] = ["ngridx","ngridy","max_nlyrs","nhydro"]
    self.dimensions["hydro_n"] = ["ngridx","ngridy","max_nlyrs","nhydro"]
    self.dimensions["hydro_reff"] = ["ngridx","ngridy","max_nlyrs","nhydro"]
    
    self.dimensions["radar_hgt"] = ["ngridx","ngridy","max_nlyrs"]
    self.dimensions["radar_prop"] = ["ngridx","ngridy","2"]

    self.dimensions["Ze"] = ["gridx","gridy","lyr","frequency","radar_npol","radar_npeaks"]
    self.dimensions["Att_hydro"] = ["gridx","gridy","lyr","frequency","att_npol"]
    self.dimensions["Att_atmo"] = ["gridx","gridy","lyr","frequency"]
    self.dimensions["tb"] = ["gridx","gridy","outlevels","angles","frequency","passive_npol"]
    
    self.units = dict()
    
    self.units["unixtime"] = "seconds since 1970-01-01 00:00"
    
    self.units["ngridx"] = "-"
    self.units["ngridy"] = "-"
    self.units["nlyrs"] = "-"
    self.units["noutlevels"] = "-"
    
    self.units["model_i"] = "-"
    self.units["model_j"] = "-"
    self.units["lat"] = "deg.dec"
    self.units["lon"] = "deg.dec"
    self.units["lfrac"] = "-"
    self.units["wind10u"] = "m/s"
    self.units["wind10v"] = "m/s"
    self.units["wind_w"] = "m/s"
    self.units["obs_height"] = "m"

    self.units["iwv"] = "kg/m^2"
    
    self.units["hgt_lev"] = "m"
    self.units["temp_lev"] = "K"
    self.units["p_lev"] = "Pa"
    self.units["relhum_lev"] = "1"
    self.units["rdar_prop"] = "dBz"

    self.units["hydro_q"] = "kg/kg"
    self.units["hydro_n"] = "-"
    self.units["hydro_reff"] = "m"
    
    self.units["radar_hgt"] = "m"
    self.units["Ze"] = "dBz"
    self.units["Att_hydros"] = "dB"
    self.units["Att_atmo"] = "dB"
    self.units["tb"] = "K"  
      
    self._nstokes = 2
    self._nangles = 16
    
    self.df = pamDescriptorFile(self)
    
    self.p = dict() #:  test
    self.r = dict()
    
    self.p["ngridx"] = 0
    self.p["ngridy"] = 0
    self.p["max_nlyrs"] = 0
  
    return
  
  
  def writeNmlFile(self,nmlFile):
    raise NotImplementedError("not yet implemented again in v1.0")
    
    #for keygr in self.nmlSet.keys():
      #for key in self.nmlSet[keygr].keys():
        #if key not in self._nmlDefaultKeys:
          #warnings.warn("Warning can not parse setting: "+str(key))
          
          
    #f = open(nmlFile,"w")
    #for keygr in self.nmlSet.keys():
      #f.write("&%s\n\r"%keygr)
      #for key in self.nmlSet[keygr].keys():
        #if self.set["pyVerbose"] > 1: print "write: ", keygr, key
        #if type(self._nmlDefaultValues[keygr][key])==bool:
          #value = str(self.nmlSet[keygr][key]).lower()
          #f.write("%s=.%s.\n\r"%(key,value,))
        #elif type(self._nmlDefaultValues[keygr][key]) in [int,np.int32,np.int64]:
          #value = int(self.nmlSet[keygr][key])
          #f.write("%s=%i\n\r"%(key,value,))  
        #elif type(self._nmlDefaultValues[keygr][key]) in [float,np.float32,np.float64]:
          #value = np.float64(self.nmlSet[keygr][key])
          #f.write("%s=%f\n\r"%(key,value,))  
        #elif type(self._nmlDefaultValues[keygr][key]) in [str]:
          #value = str(self.nmlSet[keygr][key])
          #f.write("%s='%s'\n\r"%(key,value,))
        #else:
          #raise ValueError("cannot determine type of nml key "+ key)
      #f.write("/\n\r")
    #f.close()

  def readNmlFile(self,inputFile):
    """
    read classical Pamtra Namelist File from inputFile
    """

    nmlFile = Namelist(inputFile)
    if nmlFile == {}:
      raise IOError("file not found: "+inputFile)

    for key in nmlFile.keys():
      for subkey in nmlFile[key]["par"][0].keys():
        if subkey.lower() in self._nmlDefaultSet.keys():
          if nmlFile[key]["par"][0][subkey][0] == ".true.": value = True
          elif nmlFile[key]["par"][0][subkey][0] == ".false.": value = False
          else: value = nmlFile[key]["par"][0][subkey][0]
          if type(value) != bool:
	    try: value = int(value)
	    except: 
	      try: value = float(value.replace("d","e").replace("D","e"))
	      except: pass
          self.nmlSet[subkey.lower()] = value
          if self.set["pyVerbose"] > 1: print subkey.lower(), ":", value
        elif self.set["pyVerbose"] > 0:
          print "Setting '"+ subkey.lower() +"' from '"+ inputFile +"' skipped."
    if self.set["pyVerbose"] > 1: print "reading nml file done: ", inputFile
    return
    
    
  def readPamtraProfile(self,inputFile):
    #make sure that a descriptor file was defined
    assert self.df.nhydro > 0
    
    f = open(inputFile,"r")
    g = csv.reader(f,delimiter=" ",skipinitialspace=True)
    levLay = inputFile.split(".")[-1]

    # check whether it is the new style. otherwise read the old style
    if levLay not in ["lay","lev"]:
      f.close()
      self.readClassicPamtraProfile(inputFile)
      return
      
    self.p["ngridx"], self.p["ngridy"], self.p["max_nlyrs"], self.nmlSet["noutlevels"] = [int(el) for el in g.next()[:4]]

    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape3Dout = (self.p["ngridx"],self.p["ngridy"],self.p["noutlevels"],)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,1)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,0)

    self.p["model_i"] = np.ones(self._shape2D,dtype=int) *missingIntNumber
    self.p["model_j"] = np.ones(self._shape2D,dtype=int)* missingIntNumber
    self.p["lat"] = np.ones(self._shape2D) * np.nan
    self.p["lon"] = np.ones(self._shape2D) * np.nan
    self.p["lfrac"] = np.ones(self._shape2D) * np.nan
    self.p["wind10u"] = np.ones(self._shape2D) * np.nan
    self.p["wind10v"] = np.ones(self._shape2D) * np.nan
    self.p["groundtemp"] = np.ones(self._shape2D) * np.nan
    self.p["unixtime"] = np.ones(self._shape2D,dtype=int) * missingIntNumber
    self.p["obs_height"] = np.ones(self._shape3Dout) * np.nan
    self.p["nlyrs"] = np.ones(self._shape2D,dtype=int) * missingIntNumber
    self.p["iwv"] = np.ones(self._shape2D) * np.nan
    self.p["radar_prop"] = np.ones(self._shape2D+tuple([2])) * np.nan

    self.p["hydro_q"] = np.ones(self._shape4D)  * np.nan
    self.p["hydro_n"] = np.ones(self._shape4D) * np.nan
    self.p["hydro_reff"] = np.ones(self._shape4D) * np.nan
    if (self.nmlSet["active"] and (self.nmlSet["radar_mode"] in ["moments","spectrum"])):
      self.p["airturb"] = np.ones(self._shape4D) * np.nan
    
    self.p["hgt_lev"] = np.ones(self._shape3Dplus) * np.nan
    self.p["hgt"] = np.ones(self._shape3D) * np.nan
    self.p["temp"] = np.ones(self._shape3D) * np.nan
    self.p["press"] = np.ones(self._shape3D) * np.nan
    self.p["relhum"] = np.ones(self._shape3D) * np.nan
    self.p["temp_lev"] = np.ones(self._shape3Dplus) * np.nan
    self.p["press_lev"] = np.ones(self._shape3Dplus) * np.nan
    self.p["relhum_lev"] = np.ones(self._shape3Dplus) * np.nan

    
    for xx in xrange(self._shape2D[0]):
      for yy in xrange(self._shape2D[1]):
        
        #first all the stuff without heigth dimension
        year,month,day,time, self.p["nlyrs"][xx,yy], self.p["model_i"][xx,yy], self.p["model_j"][xx,yy] = g.next()[:7]
        self.p["unixtime"][xx,yy] = calendar.timegm(datetime.datetime(year = int(year), month = int(month), day = int(day), hour = int(time[0:2]), minute = int(time[2:4]), second = 0).timetuple())
        self.p["obs_height"][xx,yy,:] = np.array(np.array(g.next()[:int(self.p["noutlevels"])]),dtype=float)
        self.p["lat"][xx,yy], self.p["lon"][xx,yy], self.p["lfrac"][xx,yy],self.p["wind10u"][xx,yy],self.p["wind10v"][xx,yy],self.p["groundtemp"][xx,yy],self.p["hgt_lev"][xx,yy,0]  = np.array(np.array(g.next()[:7]),dtype=float)
 
        self.p["iwv"][xx,yy] = np.array(np.array(g.next()[0]),dtype=float)
        
        #if levels are provided we have one line more:
        if levLay == "lev":
          dataLine = g.next()
          #in case we have spaces after the last value...
          try: dataLine = dataLine.remove("")
          except: pass
          dataLine = np.array(np.array(dataLine),dtype=float)
          self.p["hgt_lev"][xx,yy,0],self.p["press_lev"][xx,yy,0],self.p["temp_lev"][xx,yy,0],self.p["relhum_lev"][xx,yy,0] = dataLine

          #in the file its actually nlevels, so:
          self.p["nlyrs"][xx,yy] = self.p["nlyrs"][xx,yy] - 1
          
        for zz in xrange(self.p["nlyrs"][xx,yy]):
          dataLine =g.next()
          dataLineCOPY = deepcopy(dataLine)
          try: dataLine = dataLine.remove("")
          except: pass
          dataLine = map(float,dataLine)
          #do we have layer or levels
          hgt,press,temp,relhum = dataLine[:4]
          if levLay == "lay":
            self.p["hgt"][xx,yy,zz] = dataLine.pop(0)
            self.p["press"][xx,yy,zz] = dataLine.pop(0)
            self.p["temp"][xx,yy,zz] = dataLine.pop(0)
            self.p["relhum"][xx,yy,zz] = dataLine.pop(0)
          elif levLay == "lev":
            self.p["hgt_lev"][xx,yy,zz+1] = dataLine.pop(0)
            self.p["press_lev"][xx,yy,zz+1] = dataLine.pop(0)
            self.p["temp_lev"][xx,yy,zz+1] = dataLine.pop(0)
            self.p["relhum_lev"][xx,yy,zz+1] = dataLine.pop(0)
          else:
            raise IOError("Did not understand lay/lev: "+layLev)
    
          #for hydrometeors it's always layers:
          for hh in xrange(self.df.nhydro):
            if self.df.data["moment_in"][hh] == 0:
              pass #nothing to read...
            elif self.df.data["moment_in"][hh] == 1:
              self.p["hydro_n"][xx,yy,zz,hh] = dataLine.pop(0)
            elif self.df.data["moment_in"][hh] == 2:
              self.p["hydro_reff"][xx,yy,zz,hh] = dataLine.pop(0)             
            elif self.df.data["moment_in"][hh] == 3:
              self.p["hydro_q"][xx,yy,zz,hh] = dataLine.pop(0)          
            elif self.df.data["moment_in"][hh] == 12:
              self.p["hydro_n"][xx,yy,zz,hh] = dataLine.pop(0)
              self.p["hydro_reff"][xx,yy,zz,hh] = dataLine.pop(0)             
            elif self.df.data["moment_in"][hh] == 13:
              self.p["hydro_n"][xx,yy,zz,hh] = dataLine.pop(0)
              self.p["hydro_q"][xx,yy,zz,hh] = dataLine.pop(0)      
            elif self.df.data["moment_in"][hh] == 23:
              self.p["hydro_reff"][xx,yy,zz,hh] = dataLine.pop(0)
              self.p["hydro_q"][xx,yy,zz,hh] = dataLine.pop(0)   
            else:
              raise IOError ('Did not understand df.data["moment_in"]')
            if "airturb" in self.p.keys():
              self.p["airturb"][xx,yy,zz,hh] = dataLine.pop(0)   
              
          #make sure we used all the data!  
          assert len(dataLine)  == 0
                    
    f.close()
    # if layer input is used, include hgt_lev. This is required when inserting observation heights in runPamtra orr runParallelPamtra
    if levLay == "lay":
      self.p["hgt_lev"][...,1:-1] = 0.5 * (self.p["hgt"][...,:-1] + self.p["hgt"][...,1:])
      self.p["hgt_lev"][...,-1] = self.p["hgt"][...,-1] + (self.p["hgt"][...,-1] - self.p["hgt"][...,-2])*0.5
    
    return 
    
  def readClassicPamtraProfile(self,inputFile,n_moments=1):
    """
    read classical pamtra profile from file
    """
    
    f = open(inputFile,"r")
    g = csv.reader(f,delimiter=" ",skipinitialspace=True)
    year,month,day,time, self.p["ngridx"], self.p["ngridy"], self.p["nlyrs"], deltax,deltay = g.next()
    
    self.p["ngridx"], self.p["ngridy"], self.p["nlyrs"],  = int(self.p["ngridx"]), int(self.p["ngridy"]), int(self.p["nlyrs"]), 
    
    self.p["unixtime"] = calendar.timegm(datetime.datetime(year = int(year), month = int(month), day = int(day), hour = int(time[0:2]), minute = int(time[2:4]), second = 0).timetuple())
    
    self.p["ngridx"] = int(self.p["ngridx"])
    self.p["ngridy"] = int(self.p["ngridy"])
    self.p["nlyrs"] = int(self.p["nlyrs"])
    self.p["max_nlyrs"] = deepcopy(self.p["nlyrs"])
    
    
    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],4+n_moments)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],4+n_moments,1)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],4+n_moments,0)

    self.p["model_i"] = np.zeros(self._shape2D)
    self.p["model_j"] = np.zeros(self._shape2D)
    self.p["lat"] = np.zeros(self._shape2D)
    self.p["lon"] = np.zeros(self._shape2D)
    self.p["lfrac"] = np.zeros(self._shape2D)
    self.p["wind10u"] = np.zeros(self._shape2D)
    self.p["wind10v"] = np.zeros(self._shape2D)
    
    self.p["iwv"] = np.zeros(self._shape2D)
    self.p["cwp"] = np.zeros(self._shape2D)
    self.p["iwp"] = np.zeros(self._shape2D)
    self.p["rwp"] = np.zeros(self._shape2D)
    self.p["swp"] = np.zeros(self._shape2D)
    self.p["gwp"] = np.zeros(self._shape2D)
    self.p["hwp"] = np.zeros(self._shape2D)
    self.p["radar_prop"] = np.ones(self._shape2D+tuple([2])) * np.nan

    self.p["hgt_lev"] = np.zeros(self._shape3Dplus)
    self.p["temp_lev"] = np.zeros(self._shape3Dplus)
    self.p["press_lev"] = np.zeros(self._shape3Dplus)
    self.p["relhum_lev"] = np.zeros(self._shape3Dplus)
    
    self.p["cwc_q"] = np.zeros(self._shape3D)
    self.p["iwc_q"] = np.zeros(self._shape3D)
    self.p["rwc_q"] = np.zeros(self._shape3D)
    self.p["swc_q"] = np.zeros(self._shape3D)
    self.p["gwc_q"] = np.zeros(self._shape3D)
    
    if n_moments == 1:
      self.p["hwc_q"] = np.ones(self._shape3D)*missingNumber
      self.p["cwc_n"] = np.ones(self._shape3D)*missingNumber
      self.p["iwc_n"] = np.ones(self._shape3D)*missingNumber
      self.p["rwc_n"] = np.ones(self._shape3D)*missingNumber
      self.p["swc_n"] = np.ones(self._shape3D)*missingNumber
      self.p["gwc_n"] = np.ones(self._shape3D)*missingNumber
      self.p["hwc_n"] = np.ones(self._shape3D)*missingNumber
      
      for x in np.arange(self.p["ngridx"]):
        for y in np.arange(self.p["ngridy"]):
          self.p["model_i"][x,y], self.p["model_j"][x,y] = np.array(np.array(g.next()),dtype=int)
          self.p["lat"][x,y], self.p["lon"][x,y], self.p["lfrac"][x,y],self.p["wind10u"][x,y],self.p["wind10v"][x,y]  = np.array(np.array(g.next()),dtype=float)
          self.p["iwv"][x,y],self.p["cwp"][x,y],self.p["iwp"][x,y],self.p["rwp"][x,y],self.p["swp"][x,y],self.p["gwp"][x,y] = np.array(np.array(g.next()),dtype=float)
          self.p["hgt_lev"][x,y,0],self.p["press_lev"][x,y,0],self.p["temp_lev"][x,y,0],self.p["relhum_lev"][x,y,0] = np.array(np.array(g.next()),dtype=float)
          for z in np.arange(self.p["nlyrs"]):
            self.p["hgt_lev"][x,y,z+1],self.p["press_lev"][x,y,z+1],self.p["temp_lev"][x,y,z+1],self.p["relhum_lev"][x,y,z+1],self.p["cwc_q"][x,y,z],self.p["iwc_q"][x,y,z],self.p["rwc_q"][x,y,z],self.p["swc_q"][x,y,z],self.p["gwc_q"][x,y,z] = np.array(np.array(g.next()),dtype=float)

    elif n_moments == 2:
      self.p["hwc_q"] = np.zeros(self._shape3D)
      self.p["cwc_n"] = np.zeros(self._shape3D)
      self.p["iwc_n"] = np.zeros(self._shape3D)
      self.p["rwc_n"] = np.zeros(self._shape3D)
      self.p["swc_n"] = np.zeros(self._shape3D)
      self.p["gwc_n"] = np.zeros(self._shape3D)
      self.p["hwc_n"] = np.zeros(self._shape3D)
      
      for x in np.arange(self.p["ngridx"]):
        for y in np.arange(self.p["ngridy"]):
          self.p["model_i"][x,y], self.p["model_j"][x,y] = np.array(np.array(g.next()),dtype=int)
          self.p["lat"][x,y], self.p["lon"][x,y], self.p["lfrac"][x,y],self.p["wind10u"][x,y],self.p["wind10v"][x,y]  = np.array(np.array(g.next()),dtype=float)
          self.p["iwv"][x,y],self.p["cwp"][x,y],self.p["iwp"][x,y],self.p["rwp"][x,y],self.p["swp"][x,y],self.p["gwp"][x,y],self.p["hwp"][x,y] = np.array(np.array(g.next()),dtype=float)
          self.p["hgt_lev"][x,y,0],self.p["press_lev"][x,y,0],self.p["temp_lev"][x,y,0],self.p["relhum_lev"][x,y,0] = np.array(np.array(g.next()),dtype=float)
          for z in np.arange(self.p["nlyrs"]):
            self.p["hgt_lev"][x,y,z+1],self.p["press_lev"][x,y,z+1],self.p["temp_lev"][x,y,z+1],self.p["relhum_lev"][x,y,z+1],self.p["cwc_q"][x,y,z],self.p["iwc_q"][x,y,z],self.p["rwc_q"][x,y,z],self.p["swc_q"][x,y,z],self.p["gwc_q"][x,y,z],self.p["hwc_q"][x,y,z],self.p["cwc_n"][x,y,z],self.p["iwc_n"][x,y,z],self.p["rwc_n"][x,y,z],self.p["swc_n"][x,y,z],self.p["gwc_n"][x,y,z],self.p["hwc_n"][x,y,z] = np.array(np.array(g.next()),dtype=float)
    else:
      raise IOError("n_moments mus be 1 or 2")


    #in PyPamtra I want relhum not in %
    self.p["relhum_lev"] = self.p["relhum_lev"]/100.

    #finally make an array from nlyrs and unixtime
    self.p["nlyrs"] = np.ones(self._shape2D,dtype=int)*self.p["nlyrs"]
    self.p["unixtime"] = np.ones(self._shape2D,dtype=int)*self.p["unixtime"]

    f.close()

    self.p["hydro_q"] = np.zeros(self._shape4D) 
    self.p["hydro_n"] = np.zeros(self._shape4D)
    self.p["hydro_reff"] = np.zeros(self._shape4D)
    
    self.p["hydro_q"][:,:,:,0] = self.p["cwc_q"]
    self.p["hydro_q"][:,:,:,1] = self.p["iwc_q"]
    self.p["hydro_q"][:,:,:,2] = self.p["rwc_q"]
    self.p["hydro_q"][:,:,:,3] = self.p["swc_q"]
    self.p["hydro_q"][:,:,:,4] = self.p["gwc_q"]
    if n_moments == 2:
      self.p["hydro_q"][:,:,:,5] = self.p["hwc_q"]
      self.p["hydro_n"][:,:,:,0] = self.p["cwc_n"]
      self.p["hydro_n"][:,:,:,1] = self.p["iwc_n"]
      self.p["hydro_n"][:,:,:,2] = self.p["rwc_n"]
      self.p["hydro_n"][:,:,:,3] = self.p["swc_n"]
      self.p["hydro_n"][:,:,:,4] = self.p["gwc_n"]
      self.p["hydro_n"][:,:,:,5] = self.p["hwc_n"]

      
    self.p["airturb"] = np.zeros(self._shape3D)
    self.p["wind_w"] = np.zeros(self._shape3D) + np.nan
    for key in ["cwc_q","iwc_q","rwc_q","swc_q","gwc_q","hwc_q","cwc_n","iwc_n","rwc_n","swc_n","gwc_n","hwc_n"] :
      del self.p[key]
      
    return
    
  def writePamtraProfile(self,profileFile):
  
    levLay = profileFile.split(".")[-1]
    
    firstTime = datetime.datetime.utcfromtimestamp(self.p["unixtime"][0,0])
    year=str(firstTime.year)
    mon=str(firstTime.month).zfill(2)
    day=str(firstTime.day).zfill(2)
    hhmm=datetime.datetime.strftime(firstTime,"%H%M")
    
    if "iwv" not in self.p.keys():
      self.p['iwv'] = np.ones((self._shape2D[0],self._shape2D[1]))*-9999.
      self.p['hydro_wp'] = np.ones((self._shape2D[0],self._shape2D[1],self.df.nhydro))*-9999.
      self.p['hydro_tn'] = np.ones((self._shape2D[0],self._shape2D[1],self.df.nhydro))*-9999.
      #self.addIntegratedValues()

    s = str(self._shape2D[0])+" "+str(self._shape2D[1])+" "+str(self._shape3D[2])+" "+str(self._shape3Dout[2])+"\n"
    
    for xx in range(self._shape2D[0]):
      for yy in range(self._shape2D[1]):
	s += year+" "+mon+" "+day+" "+hhmm+" "+str(self._shape3D[2])+" "+str(xx+1)+" "+str(yy+1)+"\n"
	s += ' '.join(['%9e'%height for height in self.p['obs_height'][xx,yy,:]])+"\n"
	s += '%3.2f'%self.p["lat"][xx,yy]+" "+'%3.2f'%self.p["lon"][xx,yy]+" "+str(self.p["lfrac"][xx,yy])+" "+str(self.p["wind10u"][xx,yy])+" "+str(self.p["wind10v"][xx,yy])+" "+str(self.p['groundtemp'][xx,yy])+" "+str(self.p['hgt_lev'][xx,yy,0])+"\n"
	s += str(self.p["iwv"][xx,yy])
	for ihyd in range(self.df.nhydro):
	  if self.df.data['moment_in'][ihyd] == 1:
	    s +=" "+'%9e'%self.p["hydro_wp"][xx,yy,ihyd]
	  if self.df.data['moment_in'][ihyd] == 13:
	    s +=" "+'%9e'%self.p["hydro_wp"][xx,yy,ihyd]+" "+'%9e'%self.p["hydro_tn"][xx,yy,ihyd]
	s += "\n"
	if levLay == 'lev':
	  s += '%6.1f'%self.p["hgt_lev"][xx,yy,0]+" "+'%6.1f'%self.p["press_lev"][xx,yy,0]+" "+'%3.2f'%self.p["temp_lev"][xx,yy,0]+" "+'%1.4f'%(self.p["relhum_lev"][xx,yy,0])+"\n"
	  for zz in range(1,self._shape3D[2]+1):
	    s += '%6.1f'%self.p["hgt_lev"][xx,yy,zz]+" "+'%6.1f'%self.p["press_lev"][xx,yy,zz]+" "+'%3.2f'%self.p["temp_lev"][xx,yy,zz]+" "+'%1.4f'%(self.p["relhum_lev"][xx,yy,zz])+" "
	    for ihyd in range(self.df.nhydro):
	      if self.df.data['moment_in'][ihyd] == 1:
		s +=str('%9e'%self.p["hydro_q"][xx,yy,zz,ihyd])+" "
	      if self.df.data['moment_in'][ihyd] == 13:
		s +=str('%9e'%self.p["hydro_n"][xx,yy,zz,ihyd])+" "+str('%9e'%self.p["hydro_q"][xx,yy,zz,ihyd])+" "
	    s += "\n"
	elif levLay == 'lay':
	  for zz in range(0,self._shape3D[2]):
	    s += '%6.1f'%self.p["hgt"][xx,yy,zz]+" "+'%6.1f'%self.p["press"][xx,yy,zz]+" "+'%3.2f'%self.p["temp"][xx,yy,zz]+" "+'%1.4f'%(self.p["relhum"][xx,yy,zz])+" "
	    for ihyd in range(self.df.nhydro):
	      if self.df.data['moment_in'][ihyd] == 1:
		s +=str('%9e'%self.p["hydro_q"][xx,yy,zz,ihyd])+" "
	      if self.df.data['moment_in'][ihyd] == 13:
		s +=str('%9e'%self.p["hydro_n"][xx,yy,zz,ihyd])+" "+str('%9e'%self.p["hydro_q"][xx,yy,zz,ihyd])+" "
	    s += "\n"
	else:
	  raise IOError("Did not understand lay/lev: "+layLev)

    
    # write stuff to file
    f = open(profileFile, 'w')
    f.write(s)
    f.close()
    return

  def createProfile(self,**kwargs):
    '''
    Function to create Pamtra Profiles.
    
    Variables ending on _lev mean level values, variables without are layer values (vectors one entry shorter!)!
    
    Everything is needed in SI units, relhum is in Pa/PA not %
    
    The following variables are mandatroy:
    hgt_lev, (temp_lev or temp), (press_lev or press) and (relhum_lev OR relhum)
    
    The following variables are optional and guessed if not provided:  "timestamp","lat","lon","lfrac","wind10u","wind10v","hgt_lev","hydro_q","hydro_n","hydro_reff","obs_height"
    
    hydro_q, hydro_reff and hydro_n can also provided as hydro_q+no001, hydro_q+no002 etc etc

    #'''
      
      
    #we don't want any masked arrays here:
        
    for key in kwargs.keys():
      if type(kwargs[key]) == np.ma.core.MaskedArray:
        kwargs[key] = kwargs[key].filled(np.nan)  
      
      
    #make sure that an descriptor file exists already!
    if not self.df.data.shape[0] > 0:
      warnings.warn("No descriptor file defined. Assuming dry profile without hydrometeors")
    if "hydro_q" in kwargs.keys():
      assert self.df.data.shape[0] > 0
      assert self.df.nhydro == kwargs["hydro_q"].shape[-1]
    elif "hydro_n" in kwargs.keys():
      assert self.df.data.shape[0] > 0
      assert  self.df.nhydro == kwargs["hydro_n"].shape[-1]
    elif "hydro_reff" in kwargs.keys():
      assert self.df.data.shape[0] > 0
      assert self.df.nhydro == kwargs["hydro_reff"].shape[-1]
 
    allVars = self.default_p_vars  
      
    for key in kwargs.keys():
      if (key not in allVars and key.split("+")[0] not in allVars) or (key == "model_i") or (key == "model_j"):
        raise TypeError("Could not parse "+key)
    
    
    if not (("hgt_lev" in kwargs.keys() or "hgt" in kwargs.keys()) and 
        ("temp_lev" in kwargs.keys() or "temp" in kwargs.keys()) and 
        ("press_lev" in kwargs.keys() or "press" in kwargs.keys()) and 
        ("relhum_lev" in kwargs.keys() or  "relhum" in kwargs.keys())):#"q" in kwargs.keys()
      raise TypeError("I need hgt(_lev) and temp(_lev) and press(_lev) and relhum(_lev)!")
    
    if "hgt" not in kwargs.keys():
      kwargs["hgt"] = (kwargs["hgt_lev"][...,1:] + kwargs["hgt_lev"][...,:-1])/2.
    
    
    #assert self.df.nhydro > 0
  
    noDims = len(np.shape(kwargs["hgt"]))
    if noDims == 1:
      self.p["ngridx"] = 1
      self.p["ngridy"] = 1
    elif noDims == 2:
      self.p["ngridx"] = np.shape(kwargs["hgt"])[0]
      self.p["ngridy"] = 1
    elif noDims == 3:
      self.p["ngridx"] = np.shape(kwargs["hgt"])[0]
      self.p["ngridy"] = np.shape(kwargs["hgt"])[1]
    else:
      print "Too many dimensions!"
      raise IOError
    
    self.p["max_nlyrs"] = np.shape(kwargs["hgt"])[-1] 
    self.p["nlyrs"] = np.array(np.sum(kwargs["hgt"]!=missingNumber,axis=-1))
    hgtVar = "hgt"

    try:
      self.p["noutlevels"] = np.shape(kwargs["obs_height"])[-1]
    except (KeyError,IndexError):
      self.p["noutlevels"] = 2

    #if np.any(self.p["nlyrs"] != self.p["max_nlyrs"]):
      #self._radiosonde = True
    #else:
      #self._radiosonde = False
    
   
    
    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dout = (self.p["ngridx"],self.p["ngridy"],self.nmlSet["noutlevels"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin+1)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
    
    self.p["nlyrs"] = self.p["nlyrs"].reshape(self._shape2D)

    for key in ["hgt_lev","temp_lev","press_lev","relhum_lev"]:
      if key in kwargs.keys():
        self.p[key]= kwargs[key].reshape(self._shape3Dplus)

    for key in ["hgt","temp","press","relhum"]:
      if key in kwargs.keys():
        self.p[key]= kwargs[key].reshape(self._shape3D)

    if "hgt_lev" not in kwargs.keys():
      self.p["hgt_lev"] = np.empty(self._shape3Dplus)
      self.p["hgt_lev"][...,0] = self.p["hgt"][...,0] - (self.p["hgt"][...,1] - self.p["hgt"][...,0])*0.5
      self.p["hgt_lev"][...,1:-1] = 0.5 * (self.p["hgt"][...,:-1] + self.p["hgt"][...,1:])
      self.p["hgt_lev"][...,-1] = self.p["hgt"][...,-1] + (self.p["hgt"][...,-1] - self.p["hgt"][...,-2])*0.5


    self.p["model_i"] = np.array(np.where(np.logical_not(np.isnan(self.p[hgtVar][:,:,0])))[0]).reshape(self._shape2D) +1
    self.p["model_j"] = np.array(np.where(np.logical_not(np.isnan(self.p[hgtVar][:,:,0])))[1]).reshape(self._shape2D) +1

    if "timestamp" not in kwargs.keys():
      self.p["unixtime"] = np.ones(self._shape2D)* int(time.time())
      warnings.warn("timestamp set to now", Warning)
    else:
      if type(kwargs["timestamp"]) in (int,np.int32,np.int64,float,np.float32,np.float64) :
        self.p["unixtime"] = np.ones(self._shape2D,dtype=int)*kwargs["timestamp"]
      elif (type(kwargs["timestamp"]) == np.ndarray) or (type(kwargs["timestamp"]) == np.ma.masked_array):
        if kwargs["timestamp"].dtype in (int,np.int32,np.int64,float,np.float32,np.float64):
          self.p["unixtime"] = kwargs["timestamp"].reshape(self._shape2D)
        else:
          raise TypeError("timestamp entries have to be int or float objects")
      elif (type(kwargs["timestamp"]) == datetime):
        self.p["unixtime"] = np.ones(self._shape2D,dtype=int)*calendar.timegm(kwargs["timestamp"].timetuple())
      else:
        raise TypeError("timestamp has to be int, float or datetime object")

    for environment, preset in [["lat",50.938056],["lon",6.956944],["lfrac",1],["wind10u",0],["wind10v",0],["groundtemp",np.nan]]:
      if environment not in kwargs.keys():
        self.p[environment] = np.ones(self._shape2D)*preset
        warnings.warn("%s set to %s"%(environment,preset,), Warning)
      else:
        if type(kwargs[environment]) in (int,np.int32,np.int64,float,np.float32,np.float64):
          self.p[environment] = np.ones(self._shape2D) * kwargs[environment]
        else:
          self.p[environment] = kwargs[environment].reshape(self._shape2D)

    for environment, preset in [["obs_height",[833000.,0.]]]:
      if environment not in kwargs.keys():
        self.p[environment] = np.ones(self._shape3Dout)*preset
        warnings.warn("%s set to %s"%(environment,preset,), Warning)
      else:
        if type(kwargs[environment]) in (int,np.int32,np.int64,float,np.float32,np.float64):
          self.p[environment] = np.ones(self._shape3Dout) * kwargs[environment]
        else:
          self.p[environment] = kwargs[environment].reshape(self._shape3Dout)
    
    for qValue in ["hydro_q","hydro_reff","hydro_n"]:
      if qValue not in kwargs.keys()  and qValue+"+no000" not in kwargs.keys():
        self.p[qValue] = np.zeros(self._shape4D)
        warnings.warn(qValue + " set to 0", Warning)
      elif qValue+"+no000" in kwargs.keys():
        self.p[qValue] = np.zeros(self._shape4D)
        for hh in xrange(self.df.nhydro):
          self.p[qValue][...,hh] = kwargs[qValue+"+no"+str(hh).zfill(3)].reshape(self._shape3D)
      else:
        self.p[qValue] = kwargs[qValue].reshape(self._shape4D)

    for qValue in ["airturb"]:
      if qValue not in kwargs.keys():
        self.p[qValue] = np.zeros(self._shape3D)
        warnings.warn(qValue + " set to 0", Warning)
      else:
        self.p[qValue] = kwargs[qValue].reshape(self._shape3D)

    for qValue in ["wind_w"]:
      if qValue not in kwargs.keys():
        self.p[qValue] = np.zeros(self._shape3D) + np.nan
        warnings.warn(qValue + " set to nan", Warning)
      else:
        self.p[qValue] = kwargs[qValue].reshape(self._shape3D)        
              
    if "radar_prop" in kwargs.keys():
      self.p["radar_prop"] = kwargs["radar_prop"].reshape(self._shape2D+tuple([2]))
    else:
      self.p["radar_prop"] = np.zeros(self._shape2D+tuple([2]))*np.nan

    ##clean up: remove all nans
    #for key in self.p.keys():
      ##apply only to arrays, lists or tupels
      #if np.sum(np.shape(self.p[key])) > 0:
        #self.p[key][np.isnan(self.p[key])] = missingNumber
    
    return


  def _checkData(self):
    maxLimits = {"relhum_lev":200,"hydro_q":0.05,"temp_lev":320,"press_lev":110000,"relhum":200,"temp":320,"press":110000}
    minLimits = {"relhum_lev":0,"hydro_q": 0 ,"temp_lev":170,"press_lev":1,"relhum":0,"temp":170,"press":1}
    for key in self.p.keys():
      if type(self.p[key]) != np.ndarray:
        continue
      p_data = self.p[key][~np.isnan(self.p[key])]
      if key in maxLimits.keys() and len(p_data)>0:
        if np.max(p_data) > maxLimits[key]:
          raise ValueError("unrealistic value for "+ key + ": " +str(np.max(p_data)) + ", maximum is " + str(maxLimits[key]))
      if key in minLimits.keys() and len(p_data)>0:
        if len(p_data[p_data != missingNumber]) > 0:
          if np.min(p_data[p_data != missingNumber]) < minLimits[key]:
            raise ValueError("unrealistic value for "+ key + ": " +str(np.min(p_data)) + ", minimum is " + str(minLimits[key]))

  def filterProfiles(self,condition):
    '''
    discard profiles, which do not fullfill "condition" (2D boolean array)
    
    Note: If the initial field is a 2D grid, the field is flattend after application
    
    Applicate between CreateProfile and runPamtra
    
    '''
    if condition.shape != self._shape2D and condition.shape != (self._shape2D[0]*self._shape2D[1],):
      raise ValueError("shape mismatch, shape of condition must be 2D field!")
    
    condition = condition.reshape(self._shape2D)
    
    #make sure we did not removed everything!  
    assert np.sum(condition)>=1
    
    #create a new shape!
    self.p["ngridx"] = np.sum(condition)
    self.p["ngridy"] = 1
    
    try: self.df.fs_nbin =  self.df.dataFullSpec["d_ds"].shape[-1]
    except: self.df.fs_nbin = 0
    
    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dout = (self.p["ngridx"],self.p["ngridy"],self.nmlSet["noutlevels"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin+1)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
        
    for key in ["unixtime","nlyrs","lat","lon","lfrac","model_i","model_j","wind10u","wind10v","groundtemp"]:
      if key in self.p.keys(): self.p[key] = self.p[key][condition].reshape(self._shape2D)

    for key in ["obs_height"]:
      if key in self.p.keys(): self.p[key] = self.p[key][condition].reshape(self._shape3Dout)

    for key in ["hydro_q","hydro_n","hydro_reff"]:
      if key in self.p.keys(): self.p[key] = self.p[key][condition].reshape(self._shape4D)
      
    for key in ["hgt_lev","temp_lev","press_lev","relhum_lev"]:
      if key in self.p.keys(): self.p[key] = self.p[key][condition].reshape(self._shape3Dplus)
    
    for key in ["airturb",'temp', 'press', 'relhum','hgt','wind_w']:
      if key in self.p.keys(): self.p[key] = self.p[key][condition].reshape(self._shape3D)   
      
    if "radar_prop" in self.p.keys(): self.p["radar_prop"] = self.p["radar_prop"][condition].reshape(self._shape2D+tuple([2]))

    #make sure we did not forget something
    for key in self.p.keys():
      if type(self.p[key]) == np.ndarray:
        assert self.p[key].shape[0] == self.p["lat"].shape[0]
        assert self.p[key].shape[1] == self.p["lat"].shape[1]

    
 
    for key in self.df.data4D.keys():
      self.df.data4D[key] = self.df.data4D[key][condition].reshape(self._shape4D)
      
    for key in self.df.dataFullSpec.keys():
      if key == "d_bound_ds": shape5D = self._shape5Dplus
      else: shape5D = self._shape5D
      self.df.dataFullSpec[key] = self.df.dataFullSpec[key][condition].reshape(shape5D)
      

      
    return


  def tileProfiles(self,rep2D):
    '''
    repeat the profiles of Pamtra rep2D times. E.g. rep2D=(1,100) changes a Pamtra shape from (10,2) to (10,200).
    '''
    
    assert len(rep2D) == 2
    
    
    #create a new shape!
    self.p["ngridx"] = self.p["ngridx"] * rep2D[0]
    self.p["ngridy"] = self.p["ngridy"] * rep2D[1]
    
    try: self.df.fs_nbin =  self.df.dataFullSpec["d_ds"].shape[-1]
    except: self.df.fs_nbin = 0
    
    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin+1)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
        
    rep3D =  rep2D + (1,)
    rep4D =  rep2D + (1,1,)
    rep5D =  rep2D + (1,1,1,)
        
    for key in ["unixtime","nlyrs","lat","lon","lfrac","model_i","model_j","wind10u","wind10v","obs_height","groundtemp"]:
      if key in self.p.keys(): self.p[key] = np.tile(self.p[key], rep2D)

    for key in ["hydro_q","hydro_n","hydro_reff"]:
      if key in self.p.keys(): self.p[key] = np.tile(self.p[key], rep4D)
      
    for key in ["hgt_lev","temp_lev","press_lev","relhum_lev","airturb",'temp', 'press', 'relhum','hgt','wind_w',"radar_prop"]:
      if key in self.p.keys(): self.p[key] = np.tile(self.p[key], rep3D)
             
 
    for key in self.df.data4D.keys():
      self.df.data4D[key] = np.tile(self.df.data4D[key],rep4D)
      
    for key in self.df.dataFullSpec.keys():
      self.df.dataFullSpec[key] = np.tile(self.df.dataFullSpec[key],rep5D)
      

      
    return




  def addProfile(self,profile, axis=0):
    """
    Add additional dictionary "profiles" with profile information to axis "axis". Number of height bins must be equal.
    """
    for key in self.p.keys():
      if len(np.shape(self.p[key])) >= 2:
        self.p[key] = np.concatenate((self.p[key], profile[key]),axis=axis)
          
    self.p["ngridx"] = np.shape(self.p["hgt_lev"])[0]
    self.p["ngridy"] = np.shape(self.p["hgt_lev"])[1]

    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin+1)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
    
    
  def rescaleHeights(self,new_hgt_lev, new_hgt=[]):

    # sort height vectors
    old_hgt_lev = self.p["hgt_lev"]
    old_hgt = self.p["hgt"]
    if len(new_hgt) == 0: new_hgt = (new_hgt_lev[...,1:] + new_hgt_lev[...,:-1])/2.
       
#    self.p["max_nlyrs"] = len(new_hgt_lev) -1
    self.p["max_nlyrs"] = np.shape(new_hgt_lev)[2] -1
    
    #new shape
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin+1)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
    assert len(self.df.dataFullSpec.keys()) == 0

    def extrap(x, xp, yp):
      """np.interp function with linear extrapolation"""
      y = np.interp(x, xp, yp)
      y[x < xp[0]] = yp[0] + (x[x<xp[0]]-xp[0]) * (yp[0]-yp[1]) / (xp[0]-xp[1])
      y[x > xp[-1]]= yp[-1] + (x[x>xp[-1]]-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2])
      return y
  
    for key in ["hydro_q","hydro_n","hydro_reff",]:
      if key in self.p.keys():
        #make new array
        newP = np.ones(self._shape4D) * missingNumber
        for x in xrange(self._shape2D[0]):
          for y in xrange(self._shape2D[1]):
            for h in xrange(self.df.nhydro):
            #interpolate!
              newP[x,y,:,h] = extrap(new_hgt[x,y,:],old_hgt[x,y,:],self.p[key][x,y,:,h])
        #save new array
        self.p[key] = newP
        #and mark all entries below -1 as missing Number!
        self.p[key][self.p[key]<-1] = missingNumber

    for key in self.df.data4D:
      #make new array
      newP = np.ones(self._shape4D) * missingNumber
      for x in xrange(self._shape2D[0]):
        for y in xrange(self._shape2D[1]):
          for h in xrange(self.df.nhydro):
          #interpolate!
            newP[x,y,:,h] = extrap(new_hgt[x,y,:],old_hgt[x,y,:],self.df.data4D[key][x,y,:,h])
      #save new array
      self.df.data4D[key] = newP

      
    for key in ["hgt_lev","temp_lev","relhum_lev"]:
      if key in self.p.keys():
        newP = np.ones(self._shape3Dplus) * missingNumber
        for x in xrange(self._shape2D[0]):
          for y in xrange(self._shape2D[1]):
#            newP[x,y] = np.interp(new_hgt_lev,old_hgt_lev[x,y],self.p[key][x,y])
#            newP[x,y] = np.interp(new_hgt_lev[x,y],old_hgt_lev[x,y],self.p[key][x,y])
            newP[x,y] = extrap(new_hgt_lev[x,y],old_hgt_lev[x,y],self.p[key][x,y])
        self.p[key] = newP
        if key != "hgt_lev": self.p[key][self.p[key]<-1] = missingNumber
      
    for key in ["airturb","wind_w","hgt","temp","relhum"]:
      if key in self.p.keys():
        newP = np.ones(self._shape3D) * missingNumber
        for x in xrange(self._shape2D[0]):
          for y in xrange(self._shape2D[1]):
#            newP[x,y] = np.interp(new_hgt,old_hgt[x,y],self.p[key][x,y])
#            newP[x,y] = np.interp(new_hgt[x,y],old_hgt[x,y],self.p[key][x,y])
            newP[x,y] = extrap(new_hgt[x,y],old_hgt[x,y],self.p[key][x,y])
        self.p[key] = newP
        if key != "hgt": self.p[key][self.p[key]<-1] = missingNumber
      
      
    for key in ["press_lev"]:
      if key in self.p.keys():
        newP = np.ones(self._shape3Dplus) * missingNumber
        for x in xrange(self._shape2D[0]):
          for y in xrange(self._shape2D[1]):
#            newP[x,y] = np.exp(np.interp(new_hgt_lev,old_hgt_lev[x,y],np.log(self.p[key][x,y])))
#            newP[x,y] = np.exp(np.interp(new_hgt_lev[x,y],old_hgt_lev[x,y],np.log(self.p[key][x,y])))
            newP[x,y] = np.exp(extrap(new_hgt_lev[x,y],old_hgt_lev[x,y],np.log(self.p[key][x,y])))
        self.p[key] = newP
        self.p[key][self.p[key]<-1] = missingNumber
        
    for key in ["press"]:
      if key in self.p.keys():
        newP = np.ones(self._shape3D) * missingNumber
        for x in xrange(self._shape2D[0]):
          for y in xrange(self._shape2D[1]):
#            newP[x,y] = np.exp(np.interp(new_hgt,old_hgt[x,y],np.log(self.p[key][x,y])))
#            newP[x,y] = np.exp(np.interp(new_hgt[x,y],old_hgt[x,y],np.log(self.p[key][x,y])))
            newP[x,y] = np.exp(extrap(new_hgt[x,y],old_hgt[x,y],np.log(self.p[key][x,y])))
        self.p[key] = newP
        self.p[key][self.p[key]<-1] = missingNumber

        
        
    self.p["nlyrs"] = np.sum(self.p["hgt_lev"] != missingNumber,axis=-1) -1

    return
    
  def addSpectralBroadening(self,EDR,wind_uv,beamwidth_deg,integration_time,frequency,kolmogorov = 0.5):
    """
    Fills array self.p["airturb"] according to input data. Array are broadcasted to self._shape3D
    Only broadening by horizontal wind and turbuilence is considered as of now.
    
    Input:
      EDR(array) : Eddy dissipation rate (SI units)
      wind_uv(array) : horizontal wind field (SI units)
      beamwidth_deg(float) : full-width half-radiated one-way (deg)
      integration_time(float) : radar integration time
      frequency(float) : frequency used to determien lower integration bound
      kolmogorov(float,optional) : Kolmogorov constant, default 0.5

    """
    
    if "hgt" not in self.p.keys():
      self.p["hgt"] = (self.p["hgt_lev"][...,:-1] + self.p["hgt_lev"][...,1:])/2.
    
    
    heightVec = self.p["hgt"]
    lamb = 3e8/(frequency * 1e9)
    beamWidth=beamwidth_deg/2./180.*np.pi #half width half radiated http://www.wmo.int/pages/prog/gcos/documents/gruanmanuals/Z_instruments/mmcr_handbook.pdf
    L_s = (wind_uv * integration_time) + 2*heightVec*np.sin(beamWidth)#RIGHT, formular of shupe or oconnor is for full width beam width
    L_lambda = lamb / 2.
    sig_B = np.sqrt(wind_uv**2*beamWidth**2/2.76)
    sig_T = np.sqrt(3*kolmogorov/2. * (EDR/(2.*np.pi))**(2./3.) * (L_s**(2./3.) - L_lambda**(2./3.)))
    self.p["airturb"][:] = np.sqrt(sig_B**2 + sig_T**2)

    return
  
  def addIntegratedValues(self):
    
    self._shape3Dhyd = (self.p["ngridx"],self.p["ngridy"],self.df.nhydro)
    self.p['hydro_wp'] = np.zeros(self._shape3Dhyd)
    self._calcMoistRho() # provies as well dz, sum_hydro_q, and q within dict() self._helperP
    self.p['iwv'] =  np.sum(self._helperP['vapor']*self._helperP["rho_moist"]*self._helperP["dz"],axis=-1)
    #nothing to do without hydrometeors:
    if np.all(self.p['hydro_q']==0):
      self.p['hydro_wp'] = np.zeros(self._shape3Dhyd)
    else:
      for i in range(self.df.nhydro):
	self.p['hydro_wp'][...,i] = np.sum(self.p['hydro_q'][...,i]*self._helperP["rho_moist"]*self._helperP["dz"],axis=-1)

    return

  def _calcMoistRho(self):
    self._helperP = dict()
    self._helperP['dz'] = self.p['hgt_lev'][...,1::]-self.p['hgt_lev'][...,0:-1]
    self._helperP['vapor'] = rh2q(self.p['relhum']/100.,self.p['temp'],self.p['press'])
    self._helperP['sum_hydro_q'] = np.sum(self.p['hydro_q'],axis=-1)
    self._helperP['rho_moist'] = moist_rho_rh(self.p['press'],self.p['temp'],self.p['relhum']/100.,self._helperP['sum_hydro_q'])

    return
    
  def addCloudShape(self):
    """
    adds cloud base and cloud top to the data to an existing pamtra profile
    clouds are detected with rh >= 0.95 and T >= 253.15 K
    if multiple cloud layers are found, only the uppermost and the lowest one are saved.
    """
    self.p["cloudTop"] = np.ones(self._shape2D)*missingNumber
    self.p["cloudBase"] = np.ones(self._shape2D)*missingNumber
    
    for x in xrange(self._shape2D[0]):
      for y in xrange(self._shape2D[1]):
        i_top, i_base, i_cloud =  detect_liq_cloud(self.p["hgt_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["temp_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["relhum_lev"][x,y,0:self.p["nlyrs"][x,y]])
        if len(i_top)> 0:
          self.p["cloudTop"][x,y] = self.p["hgt_lev"][x,y,i_top[-1]+1]
          self.p["cloudBase"][x,y] = self.p["hgt_lev"][x,y,i_base[0]]
    return

  def addPseudoAdiabaticLWC(self):
    """
    adds liquid water content to the data to an existing pamtra profile
    using a simple cloud model with pseudo adiabatic lapse rate
    clouds are detected with rh >= 0.95 and T >= 253.15 K
    existing cwc_q values are overriden!
    """

    #calculate CWC
    
    raise NotImplementedError("not yet avaiable in pamtra v 1.0")    
    
    self.p["cwc_q"][:] = 0

    for x in xrange(self._shape2D[0]):
      for y in xrange(self._shape2D[1]):
        i_top, i_base, i_cloud =  detect_liq_cloud(self.p["hgt_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["temp_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["relhum_lev"][x,y,0:self.p["nlyrs"][x,y]])
        for k in np.arange(len(i_top)):
          i_cloud_k = np.arange(i_base[k],i_top[k]+1) 
          self.p["cwc_q"][x,y,i_cloud_k[:-1]] = mod_ad(self.p["temp_lev"][x,y,i_cloud_k], self.p["press_lev"][x,y,i_cloud_k], self.p["hgt_lev"][x,y,i_cloud_k], 1)

    self.p["cwp"][:] = 0
    
    self._calcMoistRho() #provides also temp,press,dz and q!
    self.p["cwp"][:] =  np.sum(self.p["cwc_q"]*self._helperP["rho_moist"]*self._helperP["dz"],axis=-1)

    return

  def addPIA(self,direction,applyAtmo=True, applyHydros=True):
    """
    adds two-way path integrated attenuation to result dictionary
    
    Parameters
    ----------    
    direction :{'bottom-up', 'top-down'}
      Assumed direction
    applyAtmo : bool, optional
      use attenuation of the atmosphere (default: true)
    applyHydros : bool, optional
      use attenuation of the hydrometeors (default: true)
      
    Notes
    -----
    Adds "PIA" to `self.r`
    """
        
    shapePIA = np.shape(self.r["Att_hydro"])
    if applyAtmo:
      Att_atmo = np.zeros(shapePIA)
      Att_atmo.T[:] = self.r["Att_atmo"].T
      Att_atmo[Att_atmo==missingNumber] = 0
    else: 
      Att_atmo = np.zeros(shapePIA)
    if applyHydros:
      Att_hydro = self.r["Att_hydro"]
      Att_hydro[Att_hydro==missingNumber] = 0    
    else:
      Att_hydro = np.zeros(shapePIA)
      
    self.r["PIA"] = np.zeros(shapePIA) - missingNumber
    if direction == "bottom-up":
      for hh in range(self.p["max_nlyrs"]):
        self.r["PIA"][:,:,hh,:] = Att_atmo[:,:,hh,:] +Att_hydro[:,:,hh,:] + 2*np.sum(Att_atmo[:,:,:hh,:] +Att_hydro[:,:,:hh,:],axis=2)#only half for the current layer
    elif direction == "top-down":
      for hh in range(self.p["max_nlyrs"]-1,-1,-1):
        self.r["PIA"][:,:,hh,:] = Att_atmo[:,:,hh,:] +Att_hydro[:,:,hh,:] + 2*np.sum(Att_atmo[:,:,hh+1:,:] +Att_hydro[:,:,hh+1:,:],axis=2)#only half for the current layer
    
    return

   
  def createFullProfile(self,timestamp,lat,lon,lfrac,wind10u,wind10v,
      obs_height,
      hgt_lev,press_lev,temp_lev,relhum_lev,
      hydro_q,hydro_n,hydro_reff,radar_prop):
    
    '''
    create comple PAmtra Profile
    
    No Extras, no missing values are guessed. the Data is only reshaped
    '''
    
    
    noDims = len(np.shape(hgt_lev))
    
    if noDims == 1:
      self.p["ngridx"] = 1
      self.p["ngridy"] = 1
    elif noDims == 2:
      self.p["ngridx"] = np.shape(hgt_lev)[0]
      self.p["ngridy"] = 1
    elif noDims == 3:
      self.p["ngridx"] = np.shape(hgt_lev)[0]
      self.p["ngridy"] = np.shape(hgt_lev)[1]
    else:
      print "Too many dimensions!"
      raise IOError
    
    self.p["nlyrs"] = np.sum(hgt_lev!=missingNumber,axis=-1) -1
    self.p["max_nlyrs"] = np.shape(hgt_lev)[-1] -1
    self.p["noutlevels"] = np.shape(obs_height)[-1]
    
    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dout = (self.p["ngridx"],self.p["ngridy"],self.p["noutlevels"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin+1)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
    self.p["unixtime"] = timestamp.reshape(self._shape2D)
        
    self.p["lat"] = lat.reshape(self._shape2D)
    self.p["lon"] = lon.reshape(self._shape2D)
    self.p["lfrac"] = lfrac.reshape(self._shape2D)
    self.p["model_i"] = np.array(np.where(lat.reshape(self._shape2D))[0]).reshape(self._shape2D) +1
    self.p["model_j"] = np.array(np.where(lon.reshape(self._shape2D))[1]).reshape(self._shape2D) +1
    self.p["wind10u"] = wind10u.reshape(self._shape2D)
    self.p["wind10v"] = wind10v.reshape(self._shape2D)
    self.p["obs_height"] = obs_height.reshape(self._shape3Dout)

    self.p["hydro_q"] = hydro_q.reshape(self._shape4D)
    self.p["hydro_n"] = hydro_n.reshape(self._shape4D)
    self.p["hydro_reff"] = hydro_reff.reshape(self._shape4D)
    
    self.p["hgt_lev"] = hgt_lev.reshape(self._shape3Dplus)
    self.p["temp_lev"] = temp_lev.reshape(self._shape3Dplus)
    self.p["press_lev"] = press_lev.reshape(self._shape3Dplus)
    self.p["relhum_lev"] = relhum_lev.reshape(self._shape3Dplus)    
    
    self.p["radar_prop"] = radar_prop.reshape(self._shape2D)    
    return   
    

  def runPamtra(self,freqs,checkData=True):
    '''
    run Pamtra from python
    '''
    tttt = time.time()

    
    if type(freqs) in (int,np.int32,np.int64,float,np.float32,np.float64): freqs = [freqs]
    
    self.set["freqs"] = freqs
    self.set["nfreqs"] = len(freqs)
    self.set["radar_pol"] = self.nmlSet["radar_polarisation"].split(",")
    self.set["radar_npol"] = len(self.set["radar_pol"])
    self.set["att_pol"] = ["N"]
    self.set["att_npol"] = len(self.set["att_pol"])

    assert self.set["nfreqs"] > 0
    assert self.set["radar_npol"] > 0
    assert self.set["att_npol"] > 0
    assert self.p["noutlevels"] > 0
    assert np.prod(self._shape2D)>0
    
    if checkData: self._checkData()

    if self.nmlSet["add_obs_height_to_layer"]: self._addObservationLevels()
 
    
    fortResults, self.fortError, fortObject = PamtraFortranWrapper(self.set,self.nmlSet,self.df.data,self.df.data4D,self.df.dataFullSpec,self.p)
    self.r = fortResults
    self.fortObject = fortObject
    
     
    self.r["nmlSettings"] = self.nmlSet

    if self.set["pyVerbose"] > 0: print "pyPamtra runtime:", time.time() - tttt
    return
 
  def runParallelPamtra(self,freqs,pp_local_workers="auto",pp_deltaF=1,pp_deltaX=0,pp_deltaY = 0, activeFreqs="auto", passiveFreqs="auto",checkData=True,timeout=None):
    '''
    run Pamtra from python
    '''
    
    if type(freqs) in (int,np.int32,np.int64,float,np.float32,np.float64): freqs = [freqs]

    self.set["freqs"] = freqs
    self.set["nfreqs"] = len(freqs)
    self.set["radar_pol"] = self.nmlSet["radar_polarisation"].split(",")
    self.set["radar_npol"] = len(self.set["radar_pol"])
    self.set["att_pol"] = ["N"]
    self.set["att_npol"] = len(self.set["att_pol"])

    assert self.set["nfreqs"] > 0
    assert self.set["radar_npol"] > 0
    assert self.set["att_npol"] > 0
    
    if pp_deltaF==0: pp_deltaF = self.set["nfreqs"]
    if pp_deltaX==0: pp_deltaX = self.p["ngridx"]
    if pp_deltaY==0: pp_deltaY = self.p["ngridy"]

    if hasattr(self, "fortObject"): del self.fortObject
    self.fortError = 0
    
    if pp_local_workers == "auto": pp_local_workers = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=pp_local_workers,maxtasksperchild=100)
    tttt = time.time()

    assert self.set["nfreqs"] > 0
    assert np.prod(self._shape2D)>0
    
    if checkData: self._checkData()

    if self.nmlSet["add_obs_height_to_layer"]: self._addObservationLevels()

    jobs = list()
    self.pp_resultData = list()
    pp_i = 0
    self._prepareResults()
    try: 
      for pp_startF in np.arange(0,self.set["nfreqs"],pp_deltaF):
        pp_endF = pp_startF + pp_deltaF
        if pp_endF > self.set["nfreqs"]: pp_endF = self.set["nfreqs"]
        for pp_startX in np.arange(0,self.p["ngridx"],pp_deltaX):
          pp_endX = pp_startX + pp_deltaX
          if pp_endX > self.p["ngridx"]: pp_endX = self.p["ngridx"]
          pp_ngridx = pp_endX - pp_startX
          for pp_startY in np.arange(0,self.p["ngridy"],pp_deltaY):
            pp_endY = pp_startY + pp_deltaY
            if pp_endY > self.p["ngridy"]: pp_endY = self.p["ngridy"]
            pp_ngridy = pp_endY - pp_startY
      
            if self.set["pyVerbose"] > 0: print "submitting job ", pp_i, pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY
      
            indices = [pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY]       
            profilePart, dfPart,dfPart4D,dfPartFS, settings = self._sliceProfile(*indices)
            jobs.append(pool.apply_async(parallelPamtraFortranWrapper,(
            indices,
            settings,self.nmlSet,dfPart,dfPart4D,dfPartFS,profilePart),{"returnModule":False}))#),callback=self.pp_resultData.append)
           
            
            pp_i += 1
            if self.set["pyVerbose"] > 0: print "submitted job: ", pp_i

            
      #pool.close()
      #pool.join()
    except KeyboardInterrupt:
      pool.terminate()
      pool.join()
      print "TERMINATED: KeyboardInterrupt"

    if self.set["pyVerbose"] > 0: print "waiting for all jobs to finish"
    for jj,job in enumerate(jobs):
      #import pdb;pdb.set_trace()
      try: self._joinResults(job.get(timeout=timeout))
      except multiprocessing.TimeoutError: 
        print "KILLED pool due to timeout of job", jj+1
      if self.set["pyVerbose"] > 0: print "got job", jj+1

    self.r["nmlSettings"] = self.nmlSet

    pool.terminate()
    if self.set["pyVerbose"] > 0: print "pyPamtra runtime:", time.time() - tttt
    del pool, jobs
    return    
 
  def runPicklePamtra(self,freqs,picklePath="pyPamJobs",pp_deltaF=1,pp_deltaX=0,pp_deltaY = 0, activeFreqs="auto", passiveFreqs="auto",checkData=True,timeout=None,maxWait =3600):
    '''
    run Pamtra from python
    '''
    import hashlib
    
    os.system("mkdir -p %s"%picklePath)
    
    if type(freqs) in (int,np.int32,np.int64,float,np.float32,np.float64): freqs = [freqs]

    self.set["freqs"] = freqs
    self.set["nfreqs"] = len(freqs)
    self.set["radar_pol"] = self.nmlSet["radar_polarisation"].split(",")
    self.set["radar_npol"] = len(self.set["radar_pol"])
    self.set["att_pol"] = ["N"]
    self.set["att_npol"] = len(self.set["att_pol"])

    assert self.set["nfreqs"] > 0
    assert self.set["radar_npol"] > 0
    assert self.set["att_npol"] > 0
    
    if pp_deltaF==0: pp_deltaF = self.set["nfreqs"]
    if pp_deltaX==0: pp_deltaX = self.p["ngridx"]
    if pp_deltaY==0: pp_deltaY = self.p["ngridy"]

    if hasattr(self, "fortObject"): del self.fortObject
    self.fortError = 0
    
    tttt = time.time()

    assert self.set["nfreqs"] > 0
    assert np.prod(self._shape2D)>0
    
    if checkData: self._checkData()

    if self.nmlSet["add_obs_height_to_layer"]: self._addObservationLevels()
    
    
    jobs = list()
    self.pp_resultData = list()
    pp_i = 0
    self._prepareResults()
    try: 
      for pp_startF in np.arange(0,self.set["nfreqs"],pp_deltaF):
        pp_endF = pp_startF + pp_deltaF
        if pp_endF > self.set["nfreqs"]: pp_endF = self.set["nfreqs"]
        for pp_startX in np.arange(0,self.p["ngridx"],pp_deltaX):
          pp_endX = pp_startX + pp_deltaX
          if pp_endX > self.p["ngridx"]: pp_endX = self.p["ngridx"]
          pp_ngridx = pp_endX - pp_startX
          for pp_startY in np.arange(0,self.p["ngridy"],pp_deltaY):
            pp_endY = pp_startY + pp_deltaY
            if pp_endY > self.p["ngridy"]: pp_endY = self.p["ngridy"]
            pp_ngridy = pp_endY - pp_startY
      
            if self.set["pyVerbose"] > 0: print "submitting job ", pp_i, pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY
      
            indices = [pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY]       
            profilePart, dfPart,dfPart4D,dfPartFS, settings = self._sliceProfile(*indices)
            #import pdb;pdb.set_trace()
            inputPickle = pickle.dumps((indices,settings,self.nmlSet,dfPart,dfPart4D,dfPartFS,profilePart))
            md5 = hashlib.md5(inputPickle).hexdigest()
            fname = "%s/%d_%04d_%s"%(picklePath, time.time(), pp_i, md5)
            jobs.append(fname)
            with open(fname+".job.tmp", 'w') as f:
              f.write(inputPickle)
            os.rename(fname+".job.tmp",fname+".job")
            pp_i += 1
            if self.set["pyVerbose"] > 0: print "wrote job: ", pp_i

      startTime = time.time()
      
      for mm, md5 in enumerate(jobs):
        fname = "%s.result"%(md5)
        while True:
          if ((time.time() - startTime) > maxWait):
            print("\rWaiting too long for job %d: %s"%(mm,fname))
            break
          try:
            with open(fname, 'r') as f:
              resultPickle = pickle.load(f)
            os.remove(fname)
            if resultPickle[0] is not None: 
              self._joinResults(resultPickle)
              sys.stdout.write("\rgot job %d: %s from %s"%(mm,fname,resultPickle[-1]) +" "*3 )  
              sys.stdout.flush()
            else:
              sys.stdout.write("\rjob broken %d: %s from %s"%(mm,fname,resultPickle[-1]))
              sys.stdout.flush()
            break
          except EOFError:
            sys.stdout.write("\rjob broken %d: %s from %s"%(mm,fname,resultPickle[-1]))
            sys.stdout.flush()
            break
          except (OSError, IOError):
            time.sleep(1)
            sys.stdout.write("\rwaiting for job %d: %s"%(mm,fname) +" "*3  )
            sys.stdout.flush()
    except KeyboardInterrupt:
      print "clean up"
      for mm, md5 in enumerate(jobs):
        for fname in glob.glob("%s.*"%(md5)):      
          try: 
            os.remove(fname)
            print fname, "removed."
          except OSError: continue
      raise KeyboardInterrupt
    #pool.terminate()
    if self.set["pyVerbose"] > 0: print "pyPamtra runtime:", time.time() - tttt
    #del pool, jobs
    print(" ")
    return    
 
  def runPicklePamtraSFTP(self,freqs,host,user,localPicklePath="pyPamJobs",remotePicklePath="pyPamJobs",pp_deltaF=1,pp_deltaX=0,pp_deltaY = 0, activeFreqs="auto", passiveFreqs="auto",checkData=True,timeout=None,maxWait =3600):
    '''
    run Pamtra from python
    '''
    import hashlib
    
    os.system("mkdir -p %s"%localPicklePath)
    
    cluster = sftp2Cluster(host,user)
    
    if type(freqs) in (int,np.int32,np.int64,float,np.float32,np.float64): freqs = [freqs]

    self.set["freqs"] = freqs
    self.set["nfreqs"] = len(freqs)
    self.set["radar_pol"] = self.nmlSet["radar_polarisation"].split(",")
    self.set["radar_npol"] = len(self.set["radar_pol"])
    self.set["att_pol"] = ["N"]
    self.set["att_npol"] = len(self.set["att_pol"])

    assert self.set["nfreqs"] > 0
    assert self.set["radar_npol"] > 0
    assert self.set["att_npol"] > 0
    
    if pp_deltaF==0: pp_deltaF = self.set["nfreqs"]
    if pp_deltaX==0: pp_deltaX = self.p["ngridx"]
    if pp_deltaY==0: pp_deltaY = self.p["ngridy"]

    if hasattr(self, "fortObject"): del self.fortObject
    self.fortError = 0
    
    tttt = time.time()

    assert self.set["nfreqs"] > 0
    assert np.prod(self._shape2D)>0
    
    if checkData: self._checkData()

    if self.nmlSet["add_obs_height_to_layer"]: self._addObservationLevels()

    jobs = list()
    self.pp_resultData = list()
    pp_i = 0
    self._prepareResults()
    try: 
      for pp_startF in np.arange(0,self.set["nfreqs"],pp_deltaF):
        pp_endF = pp_startF + pp_deltaF
        if pp_endF > self.set["nfreqs"]: pp_endF = self.set["nfreqs"]
        for pp_startX in np.arange(0,self.p["ngridx"],pp_deltaX):
          pp_endX = pp_startX + pp_deltaX
          if pp_endX > self.p["ngridx"]: pp_endX = self.p["ngridx"]
          pp_ngridx = pp_endX - pp_startX
          for pp_startY in np.arange(0,self.p["ngridy"],pp_deltaY):
            pp_endY = pp_startY + pp_deltaY
            if pp_endY > self.p["ngridy"]: pp_endY = self.p["ngridy"]
            pp_ngridy = pp_endY - pp_startY
      
            if self.set["pyVerbose"] > 0: print "submitting job ", pp_i, pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY
      
            indices = [pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY]       
            profilePart, dfPart,dfPart4D,dfPartFS, settings = self._sliceProfile(*indices)
            #import pdb;pdb.set_trace()
            inputPickle = pickle.dumps((indices,settings,self.nmlSet,dfPart,dfPart4D,dfPartFS,profilePart))
            md5 = hashlib.md5(inputPickle).hexdigest()
            fname = "%d_%04d_%s"%(time.time(), pp_i, md5)
            jobs.append(fname)
            cluster.put(remotePicklePath+"/"+fname+".job.tmp",inputPickle)
            cluster.mv(remotePicklePath+"/"+fname+".job.tmp",remotePicklePath+"/"+fname+".job")
            pp_i += 1
            if self.set["pyVerbose"] > 0: print "wrote job: ", pp_i

            
      del cluster
          
      startTime = time.time()
      
      for mm, md5 in enumerate(jobs):
        fname = "%s/%s.result"%(localPicklePath,md5)
        while True:
          if ((time.time() - startTime) > maxWait):
            print("\rWaiting too long for job %d: %s"%(mm,fname))
            break
          try:
            with open(fname, 'r') as f:
              resultPickle = pickle.load(f)
            os.remove(fname)
            if resultPickle[0] is not None: 
              self._joinResults(resultPickle)
              sys.stdout.write("\rgot job %d: %s from %s"%(mm,fname,resultPickle[-1]) +" "*3 )  
              sys.stdout.flush()
            else:
              sys.stdout.write("\rjob broken %d: %s from %s"%(mm,fname,resultPickle[-1]))
              sys.stdout.flush()
            break
          except EOFError:
            sys.stdout.write("\rjob broken %d: %s from %s"%(mm,fname,resultPickle[-1]))
            sys.stdout.flush()
            break
          except (OSError, IOError):
            time.sleep(1)
            sys.stdout.write("\rwaiting for job %d: %s"%(mm,fname) +" "*3  )
            sys.stdout.flush()
    except KeyboardInterrupt:
      print "clean up"
      cluster = sftp2Cluster(host,user)
      for mm, md5 in enumerate(jobs):
          cluster.rm(remotePicklePath+"/"+md5+".job")
          cluster.rm(remotePicklePath+"/"+md5+".job.tmp")
          os.remove(localPicklePath+"/"+md5+".result")
      del cluster
      raise KeyboardInterrupt
      
    #pool.terminate()
    if self.set["pyVerbose"] > 0: print "pyPamtra runtime:", time.time() - tttt
    #del pool, jobs
    print(" ")
    
    
    return

  def _addObservationLevels(self):
    """Adds observation levels to the height grid of each profile.
    Observation heights are 3-dimensional arrays (nx,ny,nout) and can be provided by either the 
    ascii input files for each profile or as profile variable self.p["obs_height"].
    
    The function is called from runPamtra and its relatives.
    """
    reduceObsTop = 0
    reduceObsBot = 0
    for i in range(np.shape(self.p["obs_height"])[2]): 
      if np.all(self.p["obs_height"][:,:,i] >= self.p["hgt_lev"][:,:,-1]):
	reduceObsTop =+1
      if np.all(self.p["obs_height"][:,:,i] <= self.p["hgt_lev"][:,:,0]):
	reduceObsBot =+1
#    import pdb; pdb.set_trace()
    new_hgt_lev = np.ones((self._shape3Dplus[0],self._shape3Dplus[1],self._shape3Dplus[2]+np.shape(self.p["obs_height"])[2]-reduceObsTop-reduceObsBot)) * np.nan
    for i in range(np.shape(self.p["obs_height"])[0]):
      for j in range(np.shape(self.p["obs_height"])[1]):
	tmp_hgt_lev = self.p["hgt_lev"][i,j,:]
	# Loop over obs heights to insert
	for k in range(reduceObsTop,np.shape(self.p["obs_height"])[2]-reduceObsBot):
#	  if (self.p["obs_height"][i,j,k] < tmp_hgt_lev[-1]) and (self.p["obs_height"][i,j,k] > tmp_hgt_lev[0]):
	    if np.any(self.p["obs_height"][i,j,k] == tmp_hgt_lev[1:-1]):
	      tmp_hgt_lev = np.sort(np.append(tmp_hgt_lev[:],self.p["obs_height"][i,j,k]+1.))
	    elif self.p["obs_height"][i,j,k] >= tmp_hgt_lev[-1]:
	      tmp_hgt_lev = np.sort(np.append(tmp_hgt_lev[:],tmp_hgt_lev[-1]-1.))
	    elif self.p["obs_height"][i,j,k] <= tmp_hgt_lev[0]:
	      tmp_hgt_lev = np.sort(np.append(tmp_hgt_lev[:],tmp_hgt_lev[0]+1.))
	    else:
	      tmp_hgt_lev = np.sort(np.append(tmp_hgt_lev[:],self.p["obs_height"][i,j,k]))
	new_hgt_lev[i,j,:] = tmp_hgt_lev
    ##      assert np.all(self.p["obs_height"][:,:,i] > self.p["hgt_lev"][:,:,0])
	  #if not np.all(self.p["obs_height"][:,:,i] > self.p["hgt_lev"][:,:,0]):
	    #print "Setting some observation heights to surface level!"
	    #self.p["obs_height"][(self.p["obs_height"][...,i] > self.p["hgt_lev"][...,0]),i] = self.p["hgt_lev"][(self.p["obs_height"][...,i] > self.p["hgt_lev"][...,0]),0]
	  #if np.any(self.p["obs_height"][:,:,i] < self.p["hgt_lev"][:,:,-1]) and np.all(self.p["obs_height"][:,:,i] > self.p["hgt_lev"][:,:,0]):
	    #new_hgt_lev = np.sort(np.concatenate((self.p["hgt_lev"],self.p["obs_height"][:,:,i].reshape(self.p["ngridx"],self.p["ngridy"],1)),axis=2),axis=2)
	    #if not np.all(np.diff(new_hgt_lev) > 0.):
	      #new_hgt_lev[np.diff(new_hgt_lev) == 0.] = new_hgt_lev[np.diff(new_hgt_lev) == 0.]-1.
	    #self.rescaleHeights(new_hgt_lev)

    self.rescaleHeights(new_hgt_lev)
    
    return
 
  def _sliceProfile(self,pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY):
    
   
    profilePart = dict()
    for key in self.p.keys():
      #import pdb;pdb.set_trace()
      if type(self.p[key]) is not np.ndarray:
        profilePart[key] = self.p[key]
      else:
        profilePart[key] = self.p[key][pp_startX:pp_endX,pp_startY:pp_endY]
        
    profilePart["ngridx"] = pp_endX - pp_startX
    profilePart["ngridy"] = pp_endY - pp_startY
       
    dfData = self.df.data
    dfData4D = dict()
    for key in self.df.data4D.keys():
      dfData4D[key] = self.df.data4D[key][pp_startX:pp_endX,pp_startY:pp_endY]
    dfDataFS = dict()
    for key in self.df.dataFullSpec.keys():
      dfDataFS[key] = self.df.dataFullSpec[key][pp_startX:pp_endX,pp_startY:pp_endY]
       
    settings = deepcopy(self.set)
    settings["nfreqs"] = pp_endF - pp_startF
    settings["freqs"] = self.set["freqs"][pp_startF:pp_endF]
       
    return profilePart, dfData, dfData4D, dfDataFS, settings
    
  def _prepareResults(self):
    
    try: maxNBin = np.max(self.df.data["nbin"]) 
    except: 
      try: 
        maxNBin = np.max(self.df.data4D["nbin"])
      except:
        maxNBin = self.df.dataFullSpec["d_ds"].shape[-1]
    radar_spectrum_length = self.nmlSet["radar_nfft"]
    
    
    self.r = dict()
    self.r["Ze"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],self.set["radar_npol"],self.nmlSet["radar_npeaks"],))*missingNumber
    self.r["Att_hydro"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],self.set["att_npol"],))*missingNumber
    self.r["Att_atmo"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
    self.r["radar_hgt"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],))*missingNumber
    if self.nmlSet["radar_mode"]=="spectrum":
      self.r["radar_spectra"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],self.set["radar_npol"],radar_spectrum_length,))*missingNumber
    else:
      self.r["radar_spectra"] = np.array([missingNumber])
    self.r["radar_snr"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],self.set["radar_npol"],self.nmlSet["radar_npeaks"],))*missingNumber
    self.r["radar_moments"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],self.set["radar_npol"],self.nmlSet["radar_npeaks"],4,))*missingNumber
    self.r["radar_slopes"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],self.set["radar_npol"],self.nmlSet["radar_npeaks"],2,))*missingNumber
    self.r["radar_edges"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],self.set["radar_npol"],self.nmlSet["radar_npeaks"],2,))*missingNumber
    self.r["radar_quality"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],self.set["radar_npol"],self.nmlSet["radar_npeaks"],),dtype=int)*missingNumber
    self.r["tb"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["noutlevels"],self._nangles*2.,self.set["nfreqs"],self._nstokes))*missingNumber
    if self.nmlSet["save_psd"]:
      self.r["psd_area"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,maxNBin))*missingNumber
      self.r["psd_n"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,maxNBin))*missingNumber
      self.r["psd_d"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,maxNBin))*missingNumber
      self.r["psd_mass"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,maxNBin))*missingNumber
    else: #save memory
      self.r["psd_area"] = np.array([missingNumber])
      self.r["psd_n"] = np.array([missingNumber])
      self.r["psd_d"] = np.array([missingNumber])
      self.r["psd_mass"] =np.array([missingNumber])
    if self.nmlSet["save_ssp"]:
      self.r["kextatmo"] = np.ones((self.p["max_nlyrs"]))*missingNumber
      self.r["scatter_matrix"] = np.ones((self.p["max_nlyrs"],self._nstokes,self._nangles,self._nstokes,self._nangles,4))*missingNumber
      self.r["extinct_matrix"] = np.ones((self.p["max_nlyrs"],self._nstokes,self._nstokes,self._nangles,2))*missingNumber
      self.r["emission_vector"] = np.ones((self.p["max_nlyrs"],self._nstokes,self._nangles,2))*missingNumber
    else: #save memory
      self.r["kextatmo"] = np.array([missingNumber])
      self.r["scatter_matrix"] = np.array([missingNumber])
      self.r["extinct_matrix"] = np.array([missingNumber])
      self.r["emission_vector"] = np.array([missingNumber])

    self.r["radar_pol"] = self.set["radar_pol"]
    self.r["att_pol"] = self.set["att_pol"]
    return 
    
  def _joinResults(self,resultList):
    '''
    Collect the data of parallel pyPamtra    
    '''  
    
    [pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY], results, fortError, host = resultList
    if self.set["pyVerbose"] > 2: print "Callback started %s:"%host

    self.fortError += fortError   
       
    self.r["pamtraVersion"] = results["pamtraVersion"]
    self.r["pamtraHash"] = results["pamtraVersion"]
    self.r["Ze"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF] = results["Ze"]
    self.r["Att_hydro"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF] = results["Att_hydro"] 
    self.r["Att_atmo"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF] = results["Att_atmo"]
    self.r["radar_hgt"][pp_startX:pp_endX,pp_startY:pp_endY,:]= results["radar_hgt"]
    self.r["tb"][pp_startX:pp_endX,pp_startY:pp_endY,:,:,pp_startF:pp_endF,:]= results["tb"]
    if self.nmlSet["radar_mode"]=="spectrum":
      self.r["radar_spectra"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF]= results["radar_spectra"]
    self.r["radar_snr"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF]= results["radar_snr"]
    self.r["radar_moments"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF]= results["radar_moments"]
    self.r["radar_slopes"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF]= results["radar_slopes"]
    self.r["radar_edges"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF]= results["radar_edges"]
    self.r["radar_quality"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF]= results["radar_quality"]
    self.r["radar_vel"]= results["radar_vel"]
    self.r["angles_deg"]= results["angles_deg"]
    if self.nmlSet["save_psd"]:
      for key in ["psd_d","psd_n","psd_mass","psd_area"]:
        self.r[key][pp_startX:pp_endX,pp_startY:pp_endY] = results[key]
    if self.nmlSet["save_ssp"]:
      for key in ["kextatmo","scatter_matrix","extinct_matrix","emission_vector"]:
        self.r[key] = results[key]
    
    return
    
  def writeResultsToNumpy(self,fname,seperateFiles=False):
    
    if not seperateFiles:
      '''
      write the complete state of the session (profile,results,settings to a file
      '''
      f = open(fname, "w")
      pickle.dump([self.r,self.p,self.nmlSet,self.set,self.df.data,self.df.data4D,self.df.dataFullSpec], f)
      f.close()
    else:
      '''
      write the complete state of the session (profile,results,settings to several files
      '''      
      os.makedirs(fname)
      for dic in ["r","p", "nmlSet","set","df.data","df.data4D","df.dataFullSpec"]:
	for key in self.__dict__[dic].keys():
	  if self.set["pyVerbose"]>1: print "saving: "+fname+"/"+dic+"%"+key+"%"+".npy"  
	  data = self.__dict__[dic][key]
	  if  type(data) == np.ma.core.MaskedArray:
	    data = data.filled(-9999)
	  
	  if type(data) in [str,OrderedDict,int,float,dict,list]:
	    np.save(fname+"/"+dic+"%"+key+"%"+".npy",data)
	  elif data.dtype == np.float64:
	    np.save(fname+"/"+dic+"%"+key+"%"+".npy",data.astype("f4"))
	  elif data.dtype == np.int64:
	    np.save(fname+"/"+dic+"%"+key+"%"+".npy",data.astype("i4"))
	  else:
	    np.save(fname+"/"+dic+"%"+key+"%"+".npy",data)
    return

    
    
  def loadResultsFromNumpy(self,fname):
    '''
    load the complete state of the session (profile,results,settings from (a) file(s)
    '''
    if os.path.isdir(fname):
      try: 
	for key in ["r","p","nmlSet","set","df.data","df.data4D","df.dataFullSpec"]:
	  self.__dict__[key] = dict()
        self.__dict__["nmlSet"] = OrderedDict()
	for fnames in os.listdir(fname):
	  key,subkey,dummy = fnames.split("%")
	  self.__dict__[key][subkey] = np.load(fname+"/"+fnames)
      except:
	print formatExceptionInfo()
	raise IOError ("Could not read data from dir")
    else:  
      try: 
	f = open(fname, "r")
	[self.r,self.p,self.nmlSet,self.set,self.df.data,self.df.data4D,self.df.dataFullSpec] = pickle.load(f)
	f.close()
      except:
	print formatExceptionInfo()	
	raise IOError ("Could not read data from file")
      
      self.df.nhydro = len(self.df.data)
      try: self.df.fs_nbin =  self.df.dataFullSpec["d_bound_ds"].shape[-1]
      except: self.df.fs_nbin = 0
      self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
      self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
      self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
      self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
      self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin+1)    
      self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
      
      return
  

  
  
  def writeResultsToNetCDF(self,fname,profileVars="all",wpNames=[],ncForm="NETCDF3_CLASSIC"):
    '''
    write the results to a netcdf file
    
    Input:
    
    fname: str filename with path
    profileVars list of variables of the profile to be saved. "all" saves all implmented ones
    ncForm: str netcdf file format, possible values are NETCDF3_CLASSIC, NETCDF3_64BIT, NETCDF4_CLASSIC, and NETCDF4 for the python-netcdf4 package (netcdf3 gives netcdf4 files for newer ubuntu versions!!") NETCDF3 takes the "old" Scientific.IO.NetCDF module, which is a bit more convinient
    '''
    
    #most syntax is identical, but there is one nasty difference regarding the fillValue...
    if ncForm in ["NETCDF3_CLASSIC", "NETCDF3_64BIT", "NETCDF4_CLASSIC", "NETCDF4"]:
      import netCDF4 as nc
      pyNc = True
    elif ncForm in ["NETCDF3"]:
      try:
        import Scientific.IO.NetCDF as nc
        pyNc = False
      except:
        #fallback for cheops with the same syntax as netcdf4!
        import netCDF3 as nc
        pyNc = True
    else:
      raise ValueError("Unknown nc form "+ncForm)
            
    if ncForm in ["NETCDF4_CLASSIC", "NETCDF4"]:  
      warnings.warn("NETCDF4 is buggy, use NETCDF3_CLASSIC if possible", Warning)
      
    try: 
      self.r
      self.r["pamtraVersion"]
    except:
      raise IOError ("run runPamtra first!")

    if pyNc: cdfFile = nc.Dataset(fname,"w",format= ncForm)
    else: cdfFile = nc.NetCDFFile(fname,"w")
    
    #write meta data
    cdfFile.history = "Created with pyPamtra (Version: "+self.r["pamtraVersion"]+", Git Hash: "+self.r["pamtraHash"]+") by "+self.nmlSet["creator"]+" (University of Cologne, IGMK) at " + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    
    cdfFile.properties = str(self.r["nmlSettings"])
    #make frequsnions
    cdfFile.createDimension('grid_x',int(self.p["ngridx"]))
    cdfFile.createDimension('grid_y',int(self.p["ngridy"]))
    cdfFile.createDimension('frequency',int(self.set["nfreqs"]))
    
    
    if (self.r["nmlSettings"]["passive"]):
      cdfFile.createDimension('angles',len(self.r["angles_deg"]))
      cdfFile.createDimension('outlevels',int(self.p["noutlevels"]))
      cdfFile.createDimension('passive_polarisation',int(self._nstokes))
      
    if (self.r["nmlSettings"]["radar_mode"] in ["spectrum","moments"]):
      cdfFile.createDimension('nfft',int(self.r["nmlSettings"]["radar_nfft"])) 
      
    if (self.r["nmlSettings"]["active"]):
      cdfFile.createDimension('heightbins',int(self.p["max_nlyrs"]))
      cdfFile.createDimension('radar_polarisation',int(self.set["radar_npol"]))
      cdfFile.createDimension('attenuation_polarisation',int(self.set["att_npol"]))
      cdfFile.createDimension('radar_peak_number',int(self.nmlSet["radar_npeaks"]))
      
    dim2d = ("grid_x","grid_y",)
    dim3d = ("grid_x","grid_y","heightbins",)
    dim3dout = ("grid_x","grid_y","outlevels",)
    dim4d = ("grid_x","grid_y","heightbins","frequency")
    dim5d_att = ("grid_x","grid_y","heightbins","frequency","attenuation_polarisation")
    dim6d_rad = ("grid_x","grid_y","heightbins","frequency","radar_polarisation","radar_peak_number")
    dim6d_rad_spec = ("grid_x","grid_y","heightbins","frequency","radar_polarisation","nfft")
    dim6d_pas = ("grid_x","grid_y","outlevels","angles","frequency","passive_polarisation")


    attUnit = "dBz"
    zeUnit = "dBz"
 
    #create and write dim variables
    
    fillVDict = dict()
    #little cheat to avoid hundreds of if, else...
    if pyNc: fillVDict["fill_value"] = missingNumber
    
    nc_frequency = cdfFile.createVariable('frequency','f',('frequency',),**fillVDict)
    nc_frequency.units = 'GHz'
    nc_frequency[:] = self.set["freqs"]
    if not pyNc: nc_frequency._fillValue =missingNumber
    
    nc_gridx = cdfFile.createVariable('grid_x','f',('grid_x',),**fillVDict)#= missingNumber)
    nc_gridx.units = '-'
    nc_gridx[:] = np.arange(self.p["ngridx"],dtype="f")
    if not pyNc: nc_gridx._fillValue =missingNumber
    
    nc_gridy = cdfFile.createVariable('grid_y','f',('grid_y',),**fillVDict)
    nc_gridy.units = '-'
    nc_gridy[:] = np.arange(self.p["ngridy"],dtype="f")
    if not pyNc: nc_gridy._fillValue =missingNumber
    
    
    
    if (self.r["nmlSettings"]["active"]):

      nc_heightbins = cdfFile.createVariable('heightbins', 'i',("heightbins",),**fillVDict)
      nc_heightbins.units = "-"
      nc_heightbins[:] = np.arange(0,self.p["max_nlyrs"],dtype="i")
      if not pyNc: nc_heightbins._fillValue =int(missingNumber)

      nc_height = cdfFile.createVariable('height', 'f',dim3d,**fillVDict)
      nc_height.units = "m"
      nc_height[:] = np.array(self.r["radar_hgt"],dtype="f")
      if not pyNc: nc_height._fillValue =missingNumber
      
      #nc_act_pol = cdfFile.createVariable('radar_polarisation', str,("radar_polarisation",))
      #nc_act_pol.units = "-"
      ##dataTmp = np.empty(self.set["radar_npol"],dtype="O")
      ##import pdb;pdb.set_trace()
      ##for dd in xrange(self.set["radar_npol"]):
        ##dataTmp[dd] = self.set["radar_pol"][dd]
      #nc_act_pol[:] = np.array(self.set["radar_pol"],dtype="O")

      #nc_att_pol = cdfFile.createVariable('attenuation_polarisation', 'char',("attenuation_polarisation"),**fillVDict)
      #nc_att_pol.units = "-"
      #nc_att_pol[:] = np.array(self.set["att_pol"],dtype="S1")
      
    if (self.r["nmlSettings"]["passive"]):
      nc_angle = cdfFile.createVariable('angles','f',('angles',),**fillVDict)
      nc_angle.units = 'deg'
      nc_angle[:] = np.array(self.r["angles_deg"],dtype="f")
      if not pyNc: nc_angle._fillValue =missingNumber
      
      nc_stokes = cdfFile.createVariable('passive_polarisation', 'S1',("passive_polarisation",),**fillVDict)
      nc_stokes.units = "-"
      nc_stokes[:] = "HV"
      
      nc_out = cdfFile.createVariable('outlevels', 'f',dim3dout,**fillVDict)
      nc_out.units = "m over sea level"
      nc_out[:] = np.array(self.p["obs_height"],dtype="f")
      if not pyNc: nc_out._fillValue =missingNumber
      
    #create and write variables
    
    nc_model_i = cdfFile.createVariable('model_i', 'i',dim2d,**fillVDict)
    nc_model_i.units = "-"
    nc_model_i[:] = np.array(self.p["model_i"],dtype="i")
    if not pyNc: nc_model_i._fillValue =int(missingNumber)
    
    nc_model_j = cdfFile.createVariable('model_j', 'i',dim2d,**fillVDict)
    nc_model_j.units = "-"
    nc_model_j[:] = np.array(self.p["model_j"],dtype="i")
    if not pyNc: nc_model_j._fillValue =int(missingNumber)
      
    nc_nlyrs = cdfFile.createVariable('nlyr', 'i',dim2d,**fillVDict)
    nc_nlyrs.units = "-"
    nc_nlyrs[:] = np.array(self.p["nlyrs"],dtype="i")
    if not pyNc: nc_nlyrs._fillValue =int(missingNumber)
    
    nc_time = cdfFile.createVariable('datatime', 'i',dim2d,**fillVDict)
    nc_time.units = "seconds since 1970-01-01 00:00:00"
    nc_time[:] = np.array(self.p["unixtime"],dtype="i")
    if not pyNc: nc_time._fillValue =int(missingNumber)
    
    nc_longitude = cdfFile.createVariable('longitude', 'f',dim2d,**fillVDict)
    nc_longitude.units = "deg.dec"
    nc_longitude[:] = np.array(self.p["lon"],dtype="f")
    if not pyNc: nc_longitude._fillValue =missingNumber
    
    nc_latitude = cdfFile.createVariable('latitude', 'f',dim2d,**fillVDict)
    nc_latitude.units = "deg.dec"
    nc_latitude[:] = np.array(self.p["lat"],dtype="f")
    if not pyNc: nc_latitude._fillValue =missingNumber
      
    nc_lfrac = cdfFile.createVariable('lfrac', 'f',dim2d,**fillVDict)
    nc_lfrac.units = "-"
    nc_lfrac[:] = np.array(self.p["lfrac"],dtype="f")
    if not pyNc: nc_lfrac._fillValue =missingNumber
    
    
    if (self.r["nmlSettings"]["active"]):
         
      nc_Ze = cdfFile.createVariable('Ze', 'f',dim6d_rad,**fillVDict)
      nc_Ze.units = zeUnit
      nc_Ze[:] = np.array(self.r["Ze"],dtype='f')
      if not pyNc: nc_Ze._fillValue =missingNumber

      nc_Attenuation_Hydrometeors = cdfFile.createVariable('Attenuation_Hydrometeors', 'f',dim5d_att,**fillVDict)
      nc_Attenuation_Hydrometeors.units = attUnit
      nc_Attenuation_Hydrometeors[:] = np.array(self.r["Att_hydro"],dtype='f')
      if not pyNc: nc_Attenuation_Hydrometeors._fillValue =missingNumber
        
      nc_Attenuation_Atmosphere = cdfFile.createVariable('Attenuation_Atmosphere', 'f',dim4d,**fillVDict)
      nc_Attenuation_Atmosphere.units = attUnit
      nc_Attenuation_Atmosphere[:] = np.array(self.r["Att_atmo"],dtype='f')
      if not pyNc: nc_Attenuation_Atmosphere._fillValue =missingNumber
      
      if "PIA" in self.r.keys():
        nc_PIA = cdfFile.createVariable('Path_Integrated_Attenuation', 'f',dim4d,**fillVDict)
        nc_PIA.units = "dB"
        nc_PIA[:] = np.array(self.r["PIA"],dtype='f')
        if not pyNc: nc_PIA._fillValue =missingNumber

      
      if ((self.r["nmlSettings"]["radar_mode"] == "spectrum") or (self.r["nmlSettings"]["radar_mode"] == "moments")):
	nc_snr=cdfFile.createVariable('Radar_SNR', 'f',dim6d_rad,**fillVDict)
	nc_snr.units="dB"
	nc_snr[:] = np.array(self.r["radar_snr"],dtype='f')
	if not pyNc: nc_snr._fillValue =missingNumber

	nc_fvel=cdfFile.createVariable('Radar_MeanDopplerVel', 'f',dim6d_rad,**fillVDict)
	nc_fvel.units="m/s"
	nc_fvel[:] = np.array(self.r["radar_moments"][...,0],dtype='f')
	if not pyNc: nc_fvel._fillValue =missingNumber
	
	nc_specw=cdfFile.createVariable('Radar_SpectrumWidth', 'f',dim6d_rad,**fillVDict)
	nc_specw.units="m/s"
	nc_specw[:] = np.array(self.r["radar_moments"][...,1],dtype='f')
	if not pyNc: nc_specw._fillValue =missingNumber
	
	nc_skew=cdfFile.createVariable('Radar_Skewness', 'f',dim6d_rad,**fillVDict)
	nc_skew.units="-"
	nc_skew[:] = np.array(self.r["radar_moments"][...,2],dtype='f')
	if not pyNc: nc_skew._fillValue =missingNumber
	
	nc_kurt=cdfFile.createVariable('Radar_Kurtosis', 'f',dim6d_rad,**fillVDict)
	nc_kurt.units="-"
	nc_kurt[:] = np.array(self.r["radar_moments"][...,3],dtype='f')
	if not pyNc: nc_kurt._fillValue =missingNumber
	
	nc_lslop=cdfFile.createVariable('Radar_LeftSlope', 'f',dim6d_rad,**fillVDict)
	nc_lslop.units="dB/(m/s)"
	nc_lslop[:] = np.array(self.r["radar_slopes"][...,0],dtype='f')
	if not pyNc: nc_lslop._fillValue =missingNumber
	
	nc_rslop=cdfFile.createVariable('Radar_RightSlope', 'f',dim6d_rad,**fillVDict)
	nc_rslop.units="dB/(m/s)"
	nc_rslop[:] = np.array(self.r["radar_slopes"][...,1],dtype='f')
	if not pyNc: nc_rslop._fillValue =missingNumber

        nc_lslop=cdfFile.createVariable('Radar_LeftEdge', 'f',dim6d_rad,**fillVDict)
        nc_lslop.units="m/s"
        nc_lslop[:] = np.array(self.r["radar_edges"][...,0],dtype='f')
        if not pyNc: nc_lslop._fillValue =missingNumber
        
        nc_rslop=cdfFile.createVariable('Radar_RightEdge', 'f',dim6d_rad,**fillVDict)
        nc_rslop.units="m/s"
        nc_rslop[:] = np.array(self.r["radar_edges"][...,1],dtype='f')
        if not pyNc: nc_rslop._fillValue =missingNumber	
	
	nc_qual=cdfFile.createVariable('Radar_Quality', 'i',dim6d_rad,**fillVDict)
	nc_qual.units="bytes"
	nc_qual.description="1st byte: aliasing; 2nd byte: 2nd peak present; 7th: no peak found"
	nc_qual[:] = np.array(self.r["radar_quality"],dtype='i')
	if not pyNc: nc_qual._fillValue =missingNumber
	
        #nc_qual=cdfFile.createVariable('Radar_Polarisation', 'i',dim6d_rad,**fillVDict)
        #nc_qual.units="bytes"
        #nc_qual.description="1st byte: aliasing; 2nd byte: 2nd peak present; 7th: no peak found"
        #nc_qual[:] = np.array(self.r["radar_quality"],dtype='i')
        #if not pyNc: nc_qual._fillValue =missingNumber	
	
	
	if ((self.r["nmlSettings"]["radar_mode"] == "spectrum")): 

	  nc_vel=cdfFile.createVariable('Radar_Velocity', 'f',("nfft",),**fillVDict)
	  nc_vel.units="m/s"
	  nc_vel[:] = np.array(self.r["radar_vel"],dtype='f')
	  if not pyNc: nc_vel._fillValue =missingNumber
	  
	  nc_spec=cdfFile.createVariable('Radar_Spectrum', 'f',dim6d_rad_spec,**fillVDict)
	  nc_spec.units="dBz"
	  nc_spec[:] = np.array(self.r["radar_spectra"],dtype='f')
	  if not pyNc: nc_spec._fillValue =missingNumber
	     
    if (self.r["nmlSettings"]["passive"]):
      nc_tb = cdfFile.createVariable('tb', 'f',dim6d_pas,**fillVDict)
      nc_tb.units = "K"
      nc_tb[:] = np.array(self.r["tb"],dtype='f')
      if not pyNc: nc_tb._fillValue =missingNumber
        
    #profile data
    
    if ("iwv" in profileVars or profileVars =="all") and ("iwv" in self.p.keys()):
      nc_iwv = cdfFile.createVariable('iwv', 'f',dim2d,**fillVDict)
      nc_iwv.units = "kg/m^2"
      nc_iwv[:] = np.array(self.p["iwv"],dtype="f")
      if not pyNc: nc_iwv._fillValue =missingNumber

    if ("hydro_wp" in profileVars or profileVars =="all") and ("hydro_wp" in self.p.keys()):
      for i,wp in enumerate(wpNames):
	nc_wp= cdfFile.createVariable(wp, 'f',dim2d,**fillVDict)
	nc_wp.units = "kg/m^2"
	nc_wp[:] = np.array(self.p["hydro_wp"][...,i],dtype="f")
	if not pyNc: nc_wp._fillValue =missingNumber

    #if ("cwp" in profileVars or profileVars =="all") and ("cwp" in self.p.keys()):
      #nc_cwp= cdfFile.createVariable('cwp', 'f',dim2d,**fillVDict)
      #nc_cwp.units = "kg/m^2"
      #nc_cwp[:] = np.array(self.p["cwp"],dtype="f")
      #if not pyNc: nc_cwp._fillValue =missingNumber
  
    #if ("iwp" in profileVars or profileVars =="all") and ("iwp" in self.p.keys()):
      #nc_iwp = cdfFile.createVariable('iwp', 'f',dim2d,**fillVDict)
      #nc_iwp.units = "kg/m^2"
      #nc_iwp[:] = np.array(self.p["iwp"],dtype="f")
      #if not pyNc: nc_iwp._fillValue =missingNumber

    #if ("rwp" in profileVars or profileVars =="all") and ("rwp" in self.p.keys()):
      #nc_rwp = cdfFile.createVariable('rwp', 'f',dim2d,**fillVDict)
      #nc_rwp.units = "kg/m^2"
      #nc_rwp[:] = np.array(self.p["rwp"],dtype="f")
      #if not pyNc: nc_rwp._fillValue =missingNumber

    #if ("swp" in profileVars or profileVars =="all") and ("swp" in self.p.keys()):
      #nc_swp = cdfFile.createVariable('swp', 'f',dim2d,**fillVDict)
      #nc_swp.units = "kg/m^2"
      #nc_swp[:] = np.array(self.p["swp"],dtype="f")
      #if not pyNc: nc_swp._fillValue =missingNumber

    #if ("gwp" in profileVars or profileVars =="all") and ("gwp" in self.p.keys()):
      #nc_gwp = cdfFile.createVariable('gwp', 'f',dim2d,**fillVDict)
      #nc_gwp.units = "kg/m^2"
      #nc_gwp[:] = np.array(self.p["gwp"],dtype="f")
      #if not pyNc: nc_gwp._fillValue =missingNumber

    #if ("hwp" in profileVars or profileVars =="all") and ("hwp" in self.p.keys()):
      #nc_hwp = cdfFile.createVariable('hwp', 'f',dim2d,**fillVDict)
      #nc_hwp.units = "kg/m^2"
      #nc_hwp[:] = np.array(self.p["hwp"],dtype="f")
      #if not pyNc: nc_hwp._fillValue =missingNumber

    if ("cloudBase" in profileVars or profileVars =="all") and ("cloudBase" in self.p.keys()):
      nc_cb = cdfFile.createVariable('cloudBase', 'f',dim2d,**fillVDict)
      nc_cb.units = "m"
      nc_cb[:] = np.array(self.p["cloudBase"],dtype="f")
      if not pyNc: nc_cb._fillValue =missingNumber

    if ("cloudTop" in profileVars or profileVars =="all") and ("cloudTop" in self.p.keys()):
      nc_ct = cdfFile.createVariable('cloudTop', 'f',dim2d,**fillVDict)
      nc_ct.units = "m"
      nc_ct[:] = np.array(self.p["cloudTop"],dtype="f")
      if not pyNc: nc_ct._fillValue =missingNumber

    cdfFile.close()
    if self.set["pyVerbose"] > 0: print fname,"written"
    
  def averageResultTbs(self,translatorDict):
    """
    average several frequencies of the passive observations to account for channel width. Replaces self.r["tb"]
    
    Input:
    translatorDict: dict with list of old and new frequencies. Note that the "new" frequency is only for naming of the channel. e.g.
    
    translatorDict = {
      90.0: [90.0],
      120.15: [117.35, 120.15],
    }
    
    
    As of know this works only in passive mode. #
    The old tbs and freqs are stored in self.r["not_averaged_tb"] and self.set["not_averaged_freqs"].
    """
    
    assert not self.nmlSet["active"]
    
    
    self.r["not_averaged_tb"]  = deepcopy(self.r["tb"])
    self.set["not_averaged_freqs"] = deepcopy(self.set["freqs"])
    
    self.set["freqs"] = sorted(translatorDict.keys())
    self.set["nfreqs"] = len(self.set["freqs"])
    self.r["tb"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["noutlevels"],self._nangles*2.,self.set["nfreqs"],self._nstokes))*missingNumber

    for ff, (freqNew, freqList) in enumerate(sorted(translatorDict.items())):
      assert np.where(np.array(self.set["freqs"]) == freqNew)[0][0] == ff
      indices = []
      for freq in freqList:
        indexFound = np.where(freq == np.array(self.set["not_averaged_freqs"]))[0]
        assert len(indexFound) == 1
        indices.append(indexFound[0])
      if len(indices) > 1:
        self.r["tb"][:,:,:,:,ff,:] = np.mean(self.r["not_averaged_tb"][:,:,:,:,indices,:],axis=4)
      else:
        self.r["tb"][:,:,:,:,ff,:] = self.r["not_averaged_tb"][:,:,:,:,indices[0],:]


  def dumpForRT4(self,runFile):
    """
    Write a files to run the passive model with the pure RT4 code
    by Frank Evans.
    """
    def extrap(x, xp, yp):
      """np.interp function with linear extrapolation"""
      y = np.interp(x, xp, yp)
      y[x < xp[0]] = yp[0] + (x[x<xp[0]]-xp[0]) * (yp[0]-yp[1]) / (xp[0]-xp[1])
      y[x > xp[-1]]= yp[-1] + (x[x>xp[-1]]-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2])
      return y

    # atmosphere is height[km] temp[K] gaseaous extinction scatfile for the low below
    s = ''
    for xx in range(self._shape2D[0]):
      for yy in range(self._shape2D[1]):
#	self.p['temp_lev'][xx,yy,:] = extrap(self.p['hgt_lev'][xx,yy,:],self.p['hgt'][xx,yy,:],self.p['temp'][xx,yy,:])
	ge = np.zeros((self._shape3D[2]))
	ge[1:self._shape3D[2]] = self.r["kextatmo"]*1.e3
	mu = (np.cos(np.deg2rad(self.r['angles_deg'][16::])))
	sm = self.r['scatter_matrix']*1.e3
	em = self.r['extinct_matrix']*1.e3
	ev = self.r['emis_vector']*1.e3
	for zz in range(self._shape3D[2])[::-1]:
	  if zz < self._shape3D[2]:
	    if (np.all(self.p['hydro_q'][xx,yy,zz-1] == 0.)):
	      scat_file = '            '
	    else:
	      scat_file = 'scat_'+'%03d'%(zz-1)+'.txt'
	      ns = '%.12e'%0.0
	      sf = open(scat_file,'w')
	      ss = ''
	      ss += '  16   0 \'LOBATTO        \'\n\n'
	      for i in range(16):
		for j in range(16):
		  ss += '%.9f'%mu[i]+'    '+'%.9f'%mu[j]+'    0'+'\n'
		  for k in range(2):
		    ss += '%.12e'%sm[zz-1,k,j,0,i,0]+' '+'%.12e'%sm[zz-1,k,j,1,i,0]+' '+ns+' '+ns+'\n'
		  ss += ns+' '+ns+' '+ns+' '+ns+'\n'
		  ss += ns+' '+ns+' '+ns+' '+ns+'\n'
		for j in range(16):
		  ss += '%.9f'%mu[i]+'    '+'%.9f'%(-1.*mu[j])+'    0'+'\n'
		  for k in range(2):
		    ss += '%.12e'%sm[zz-1,k,j,0,i,2]+' '+'%.12e'%sm[zz-1,k,j,1,i,2]+' '+ns+' '+ns+'\n'
		  ss += ns+' '+ns+' '+ns+' '+ns+'\n'
		  ss += ns+' '+ns+' '+ns+' '+ns+'\n'
	      for i in range(16):
		for j in range(16):
		  ss += '%.9f'%(-1.*mu[i])+'    '+'%.9f'%mu[j]+'    0'+'\n'
		  for k in range(2):
		    ss += '%.12e'%sm[zz-1,k,j,0,i,1]+' '+'%.12e'%sm[zz-1,k,j,1,i,1]+' '+ns+' '+ns+'\n'
		  ss += ns+' '+ns+' '+ns+' '+ns+'\n'
		  ss += ns+' '+ns+' '+ns+' '+ns+'\n'
		for j in range(16):
		  ss += '%.9f'%(-1.*mu[i])+'    '+'%.9f'%(-1.*mu[j])+'    0'+'\n'
		  for k in range(2):
		    ss += '%.12e'%sm[zz-1,k,j,0,i,3]+' '+'%.12e'%sm[zz-1,k,j,1,i,3]+' '+ns+' '+ns+'\n'
		  ss += ns+' '+ns+' '+ns+' '+ns+'\n'
		  ss += ns+' '+ns+' '+ns+' '+ns+'\n'
	      ss += '\n'
	      for i in range(16):
		ss += '%.9f'%mu[i]+'\n'
		for k in range(2):
		  ss += '%.12e'%em[zz-1,k,0,i,0]+' '+'%.12e'%em[zz-1,k,1,i,0]+' '+ns+' '+ns+'\n'
		ss += ns+' '+ns+' '+ns+' '+ns+'\n'
		ss += ns+' '+ns+' '+ns+' '+ns+'\n'
	      for i in range(16):
		ss += '%.9f'%(-1.*mu[i])+'\n'
		for k in range(2):
		  ss += '%.12e'%em[zz-1,k,0,i,1]+' '+'%.12e'%em[zz-1,k,1,i,1]+' '+ns+' '+ns+'\n'
		ss += ns+' '+ns+' '+ns+' '+ns+'\n'
		ss += ns+' '+ns+' '+ns+' '+ns+'\n'
	      ss += '\n'
	      for i in range(16):
		ss += '%.9f'%mu[i]+'  '
		ss += '%.12e'%ev[zz-1,0,i,0]+' '+'%.12e'%ev[zz-1,1,i,0]+' '+ns+' '+ns+'\n'
	      for i in range(16):
		ss += '%.9f'%(-1.*mu[i])+'  '
		ss += '%.12e'%ev[zz-1,0,i,1]+' '+'%.12e'%ev[zz-1,1,i,1]+' '+ns+' '+ns+'\n'
	      sf.write(ss)
	      sf.close()
	  else:
	    scat_file = '            '
	  s += '%3.2f'%(self.p["hgt_lev"][xx,yy,zz]/1.e3)+" "+'%3.2f'%self.p["temp_lev"][xx,yy,zz]+" "+'%3.6f'%ge[zz]+'  \''+scat_file+'\'\n'
    
    
    #firstTime = datetime.datetime.utcfromtimestamp(self.p["unixtime"][0,0])
    #year=str(firstTime.year)
    #mon=str(firstTime.month).zfill(2)
    #day=str(firstTime.day).zfill(2)
    #hhmm=datetime.datetime.strftime(firstTime,"%H%M")
    
    #if "iwv" not in self.p.keys():
      #self.p['iwv'] = np.ones((self._shape2D[0],self._shape2D[1]))*-9999.
      #self.p['hydro_wp'] = np.ones((self._shape2D[0],self._shape2D[1],self.df.nhydro))*-9999.
      #self.p['hydro_tn'] = np.ones((self._shape2D[0],self._shape2D[1],self.df.nhydro))*-9999.
      ##self.addIntegratedValues()

    #s = str(self._shape2D[0])+" "+str(self._shape2D[1])+" "+str(self._shape3D[2])+" "+str(self._shape3Dout[2])+"\n"
    
    #for xx in range(self._shape2D[0]):
      #for yy in range(self._shape2D[1]):
	#s += year+" "+mon+" "+day+" "+hhmm+" "+str(self._shape3D[2])+" "+str(xx+1)+" "+str(yy+1)+"\n"
	#s += ' '.join(['%9e'%height for height in self.p['obs_height'][xx,yy,:]])+"\n"
	#s += '%3.2f'%self.p["lat"][xx,yy]+" "+'%3.2f'%self.p["lon"][xx,yy]+" "+str(self.p["lfrac"][xx,yy])+" "+str(self.p["wind10u"][xx,yy])+" "+str(self.p["wind10v"][xx,yy])+" "+str(self.p['groundtemp'][xx,yy])+" "+str(self.p['hgt_lev'][xx,yy,0])+"\n"
	#s += str(self.p["iwv"][xx,yy])
	#for ihyd in range(self.df.nhydro):
	  #if self.df.data['moment_in'][ihyd] == 1:
	    #s +=" "+'%9e'%self.p["hydro_wp"][xx,yy,ihyd]
	  #if self.df.data['moment_in'][ihyd] == 13:
	    #s +=" "+'%9e'%self.p["hydro_wp"][xx,yy,ihyd]+" "+'%9e'%self.p["hydro_tn"][xx,yy,ihyd]
	#s += "\n"
	#if levLay == 'lev':
	  #s += '%6.1f'%self.p["hgt_lev"][xx,yy,0]+" "+'%6.1f'%self.p["press_lev"][xx,yy,0]+" "+'%3.2f'%self.p["temp_lev"][xx,yy,0]+" "+'%1.4f'%(self.p["relhum_lev"][xx,yy,0])+"\n"
	  #for zz in range(1,self._shape3D[2]+1):
	    #s += '%6.1f'%self.p["hgt_lev"][xx,yy,zz]+" "+'%6.1f'%self.p["press_lev"][xx,yy,zz]+" "+'%3.2f'%self.p["temp_lev"][xx,yy,zz]+" "+'%1.4f'%(self.p["relhum_lev"][xx,yy,zz])+" "
	    #for ihyd in range(self.df.nhydro):
	      #if self.df.data['moment_in'][ihyd] == 1:
		#s +=str('%9e'%self.p["hydro_q"][xx,yy,zz,ihyd])+" "
	      #if self.df.data['moment_in'][ihyd] == 13:
		#s +=str('%9e'%self.p["hydro_n"][xx,yy,zz,ihyd])+" "+str('%9e'%self.p["hydro_q"][xx,yy,zz,ihyd])+" "
	    #s += "\n"
	#elif levLay == 'lay':
	  #for zz in range(0,self._shape3D[2]):
	    #s += '%6.1f'%self.p["hgt"][xx,yy,zz]+" "+'%6.1f'%self.p["press"][xx,yy,zz]+" "+'%3.2f'%self.p["temp"][xx,yy,zz]+" "+'%1.4f'%(self.p["relhum"][xx,yy,zz])+" "
	    #for ihyd in range(self.df.nhydro):
	      #if self.df.data['moment_in'][ihyd] == 1:
		#s +=str('%9e'%self.p["hydro_q"][xx,yy,zz,ihyd])+" "
	      #if self.df.data['moment_in'][ihyd] == 13:
		#s +=str('%9e'%self.p["hydro_n"][xx,yy,zz,ihyd])+" "+str('%9e'%self.p["hydro_q"][xx,yy,zz,ihyd])+" "
	    #s += "\n"
	#else:
	  #raise IOError("Did not understand lay/lev: "+layLev)

    
    # write stuff to file
    f = open(runFile, 'w')
    f.write(s)
    f.close()
    return
    
