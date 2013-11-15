# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import datetime
import csv
import pickle

import time,calendar,datetime
import warnings
import sys
import traceback
import os
import random
import string
import itertools
from copy import deepcopy
from matplotlib import mlab
import multiprocessing

import namelist #parser for fortran namelist files

import meteoSI
try: 
  import pyPamtraLibWrapper 
except: 
  warnings.warn("pyPamtraLib not available", Warning)



missingNumber=-9999.
#logging.basicConfig(filename='example.log',level=logging.DEBUG)



class pamDescriptorFile(object):
  #class with the descriptor file content in data. In case you want to use 4D data, use data4D instead, the coreesponding column in data is automatically removed.
  
         
  def __init__(self,parent):
    self.names =np.array(["hydro_name", "as_ratio", "liq_ice", "rho_ms", "a_ms", "b_ms", "alpha_as", "beta_as", "moment_in", "nbin", "dist_name", "p_1", "p_2", "p_3", "p_4", "d_1", "d_2", "scat_name", "vel_size_mod"])
    self.types = ["S15",float,int,float,float,float,float,float,int,int,"S15",float,float,float,float,float,float, "S15", "S15"]  
    self.data = np.recarray((0,),dtype=zip(self.names, self.types))
    self.data4D = dict()
    self.nhydro = 0
    self.parent = parent
    self.dataFullSpec = dict()
    self.fs_nbin = 0
    return
    
   
  def readFile(self,fname):
    f = open(fname,"r")
    for row in csv.reader(f,delimiter=" ",skipinitialspace=True):
      #skipp comments
      if row[0][0] == "!":
        continue
      #make sure line is complete
      assert len(row) == 19  

      for ii, item in enumerate(row):
      #python does not like double types
        if self.types[ii] == float:
          row[ii] = item.replace("d", "e")
      #fortran strips quotes
        if self.types[ii] == "S15":
          row[ii] = item.replace("'", "").replace('"', '')
  
      self.addHydrometeor(row)
      print ', '.join(row), len(row)
      
    f.close()  
    #make rec array

  def writeFile(self,fname):
    assert len(data[0]) == len(self.names)
    return mlab.rec2csv(data, fname, delimiter=' ',withheader=False)
      
  def addHydrometeor(self,hydroTuple):
    assert hydroTuple[0] not in self.data["hydro_name"]
    assert len(self.dataFullSpec.keys()) == 0
    
    self.data = np.append(self.data,np.array(tuple(hydroTuple),dtype=zip(self.names,self.types)))
    self.nhydro += 1
    self.parent._shape4D = (self.parent.p["ngridx"],self.parent.p["ngridy"],self.parent.p["max_nlyrs"],self.nhydro)
    for key in ["hydro_q","hydro_reff","hydro_n"]:
      if key in self.parent.p.keys():    
        self.parent.p[key] = np.concatenate([self.parent.p[key],np.ones(self.parent._shape3D + tuple([1]))*missingNumber],axis=-1)
    return
    
  def removeHydrometeor(self,hydroName):
    assert len(self.dataFullSpec.keys()) == 0
    if hydroName == "all":
      self.__init__()
      return
    removed = False#
    hydroIndex = np.where(self.parent.df.data["hydro_name"]==hydroName)[0][0]
    self.dataNew = np.recarray((0,),dtype=zip(self.names, self.types))
    for ii in range(self.data.shape[0]):
      if self.data[ii][0] == hydroName:
        removed = True
        continue
      else:
        self.dataNew = np.append(self.dataNew,self.data[ii])
    self.data = self.dataNew
    
    if not removed:
      raise ValueError("Did not find "+hydroName)
    else:
      self.nhydro -= 1
      self.parent._shape4D = (self.parent.p["ngridx"],self.parent.p["ngridy"],self.parent.p["max_nlyrs"],self.nhydro)
      for key in ["hydro_q","hydro_reff","hydro_n", "nbin"]:
        if key in self.parent.p.keys():
          self.parent.p[key] = np.delete(self.parent.p[key],hydroIndex,axis=-1)
    return

    
  def add4D(self, key, arr):
    print "changing", key, "to 4D"
    self.data = mlab.rec_drop_fields(self.data,[key])
    self.data4D[key] = arr.reshape(self.parent._shape4D)
    
  def remove4D(self, key, val):
    self.data = mlab.rec_append_fields(self.data,key,val)
    del data4D[key]
    
  def addFullSpectra(self):
    #assert len(self.data4D.keys()) == 0
    assert len(self.data) >0
    assert np.min(self.data["nbin"]) == np.max(self.data["nbin"])

    self.fs_nbin =  np.max(self.data["nbin"])
    
    self.dataFullSpec["delta_d_ds"] = np.ones(self.parent._shape4D+(self.fs_nbin-1,))
    self.dataFullSpec["density2scat"] = np.ones(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["diameter2scat"] = np.ones(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["d_bound_ds"] = np.ones(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["f_ds"] = np.ones(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["mass_ds"] = np.ones(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["area_ds"] = np.ones(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["as_ratio"] = np.ones(self.parent._shape4D+(self.fs_nbin,))
  
    #for name in self.names:
      #if name not in ["hydro_name","liq_ice","scat_name","vel_size_mod"]:
        #self.data = mlab.rec_drop_fields(self.data,[name])
    #return
    
class pyPamtra(object):

  '''
  Class for pamtra calculations. Initalisations fill dictonary 'set' with default values and 
  also fills 'dimensions' and 'units'
  
  '''
  def __init__(self):
    #set setting default values
  
  
    self._prepareNmlUnitsDimensions()
  
    self._nstokes = 2
    self._noutlevels = 2
    self._nangles = 32
    
    self.df = pamDescriptorFile(self)
    
    self.p = dict()
    self.r = dict()
    
    self.p["ngridx"] = 0
    self.p["ngridy"] = 0
    self.p["max_nlyrs"] = 0
    
    return
  
  def _prepareNmlUnitsDimensions(self):
  
    self.default_p_vars = ["timestamp","lat","lon","lfrac","wind10u","wind10v","hgt_lev","press_lev","temp_lev","relhum_lev","q","hydro_q","hydro_n","hydro_reff","wind10u","wind10v","obs_height", "ngridy","ngridx","max_nlyrs","nlyrs","model_i","model_j","unixtime"]
  
    self.nmlSet = dict() #settings which are required for the nml file. keeping the order is important for fortran
    self.nmlSet["hydro_threshold"]=  1.e-10   # [kg/kg] 
    #set namelist defaults#
    # sec inoutput_mode
    self.nmlSet["write_nc"]= True
    self.nmlSet["dump_to_file"]= False
    self.nmlSet["tmp_path"]= '/tmp/'
    self.nmlSet["data_path"]= 'data/'
    self.nmlSet["crm_case"]= ''
    self.nmlSet["crm_data"]= ''
    self.nmlSet["crm_data2"]= ''
    self.nmlSet["crm_constants"]= ''
    self.nmlSet["jacobian_mode"]= False #profile 1,1 is reference, for all other colums only layers with different values are calculated
    self.nmlSet["save_psd"]= False #also saves the PSDs used for radiative transfer
    # sec output
    self.nmlSet["obs_height"]= 833000.
    self.nmlSet["units"]= 'T'
    self.nmlSet["outpol"]= 'VH'
    self.nmlSet["file_desc"]= ''
    self.nmlSet["creator"]= 'Pamtrauser'
    # sec run_mode
    self.nmlSet["active"]= True
    self.nmlSet["passive"]= True
    self.nmlSet["radar_mode"]= "simple" #"splitted"|"moments"|"spectrum"
    # sec surface params
    self.nmlSet["ground_type"]= 'S'
    self.nmlSet["salinity"]= 33.0
    self.nmlSet["emissivity"]= 0.6
    # sec gas_abs_mod
    self.nmlSet["lgas_extinction"]= True
    self.nmlSet["gas_mod"]= 'R98'
    # sec hyd_opts
    self.nmlSet["lhyd_extinction"]= True
    self.nmlSet["lphase_flag"]=  True
    self.nmlSet["hydro_fullspec"] = False
    self.nmlSet["hydro_limit_density_area"] = True
    # radar_simulator
    #number of FFT points in the Doppler spectrum [typically 256 or 512]
    self.nmlSet["radar_nfft"]= 256
    #number of average spectra for noise variance reduction, typical range [1 150]
    self.nmlSet["radar_no_ave"]= 150
    #MinimumNyquistVelocity in m/sec
    self.nmlSet["radar_max_v"]= 7.885
    #MaximumNyquistVelocity in m/sec
    self.nmlSet["radar_min_v"]= -7.885
    #radar noise offset in same unit as Ze 10*log10(mm⁶/m³). noise is calculated with noise"]=  radar_pnoise0 + 20*log10(range)
    self.nmlSet["radar_pnoise0"]= -84.031043312334901 # value for BArrow MMCR for 2008, 4, 15,
    self.nmlSet["radar_airmotion"]=  False
    self.nmlSet["radar_airmotion_model"]=  "step" #"constant","linear","step"
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
    self.nmlSet["radar_npeaks"] = 3
    self.nmlSet["radar_noise_distance_factor"]=  2.0
    self.nmlSet["radar_receiver_uncertainty_std"]=  0.e0 #dB

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
    self.dimensions["obs_height"] = ["ngridx","ngridy"]

    self.dimensions["iwv"] = ["ngridx","ngridy"]
    self.dimensions["cwp"] = ["ngridx","ngridy"]
    self.dimensions["iwp"] = ["ngridx","ngridy"]
    self.dimensions["rwp"] = ["ngridx","ngridy"]
    self.dimensions["swp"] = ["ngridx","ngridy"]
    self.dimensions["gwp"] = ["ngridx","ngridy"]
    self.dimensions["hwp"] = ["ngridx","ngridy"]
    
    self.dimensions["hgt_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
    self.dimensions["temp_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
    self.dimensions["p_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
    self.dimensions["relhum_lev"] = ["ngridx","ngridy","max_nlyrs+1"]
    
    self.dimensions["hydro_q"] = ["ngridx","ngridy","max_nlyrs","nhydro"]
    self.dimensions["hydro_n"] = ["ngridx","ngridy","max_nlyrs","nhydro"]
    self.dimensions["hydro_reff"] = ["ngridx","ngridy","max_nlyrs","nhydro"]
    
    self.dimensions["radar_hgt"] = ["ngridx","ngridy","max_nlyrs"]
    self.dimensions["Ze"] = ["gridx","gridy","lyr","frequency"]
    self.dimensions["Att_hydros"] = ["gridx","gridy","lyr","frequency"]
    self.dimensions["Att_atmo"] = ["gridx","gridy","lyr","frequency"]
    self.dimensions["tb"] = ["gridx","gridy","outlevels","angles","frequency","stokes"]

    
    self.units = dict()
    
    self.units["unixtime"] = "seconds since 1970-01-01 00:00"
    
    self.units["ngridx"] = "-"
    self.units["ngridy"] = "-"
    self.units["nlyrs"] = "-"
    
    self.units["model_i"] = "-"
    self.units["model_j"] = "-"
    self.units["lat"] = "deg.dec"
    self.units["lon"] = "deg.dec"
    self.units["lfrac"] = "-"
    self.units["wind10u"] = "m/s"
    self.units["wind10v"] = "m/s"
    self.units["obs_height"] = "m"

    self.units["iwv"] = "kg/m^2"
    self.units["cwp"] = "kg/m^2"
    self.units["iwp"] = "kg/m^2"
    self.units["rwp"] = "kg/m^2"
    self.units["swp"] = "kg/m^2"
    self.units["gwp"] = "kg/m^2"
    self.units["hwp"] = "kg/m^2"
    
    self.units["hgt_lev"] = "m"
    self.units["temp_lev"] = "K"
    self.units["p_lev"] = "Pa"
    self.units["relhum_lev"] = "1"
    
    self.units["hydro_q"] = "kg/kg"
    self.units["hydro_n"] = "-"
    self.units["hydro_reff"] = "m"
    
    self.units["radar_hgt"] = "m"
    self.units["Ze"] = "dBz OR mm^6 m^-3"
    self.units["Att_hydros"] = "dB OR linear"
    self.units["Att_atmo"] = "dB OR linear"
    self.units["tb"] = "K"  
  
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

    nmlFile = namelist.Namelist(inputFile)
    if nmlFile == {}:
      raise IOError("file not found: "+inputFile)

    for key in nmlFile.keys():
      for subkey in nmlFile[key]["par"][0].keys():
        if subkey.lower() in self._nmlDefaultSet.keys():
          if nmlFile[key]["par"][0][subkey][0] == ".true.": value = True
          elif nmlFile[key]["par"][0][subkey][0] == ".false.": value = False
          else: value = nmlFile[key]["par"][0][subkey][0]
          self.nmlSet[subkey.lower()] = value
          if self.set["pyVerbose"] > 1: print subkey.lower(), ":", value
        elif self.set["pyVerbose"] > 0:
          print "Setting '"+ subkey.lower() +"' from '"+ inputFile +"' skipped."
    if self.set["pyVerbose"] > 1: print "reading nml file done: ", inputFile
    return
    
  def readPamtraProfile(self,inputFile,n_moments=1):
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
    
    if self.p["max_nlyrs"] > 200:
      warnings.warn("Too many layers for pamtra (max:200): " + str(self.p["max_nlyrs"]),Warning)
    
    
    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],4+n_moments)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],4+n_moments,0)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],4+n_moments,1)

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
    for key in ["cwc_q","iwc_q","rwc_q","swc_q","gwc_q","hwc_q","cwc_n","iwc_n","rwc_n","swc_n","gwc_n","hwc_n"] :
      del self.p[key]
      
    return
    
  def writePamtraProfile(self,filename):
    
    #the ASCII format has no support for changing dates, thus tkae the first one for all!
    raise NotImplementedError("not yet implemented again in v1.0")
    
    
    #firstTime = datetime.datetime.utcfromtimestamp(self.p["unixtime"][0,0])
    #year=str(firstTime.year)
    #mon=str(firstTime.month).zfill(2)
    #day=str(firstTime.day).zfill(2)
    #hhmm=datetime.datetime.strftime(firstTime,"%H%M")
    
    #if "iwv" not in self.p.keys():
      #self.addIntegratedValues()

    #s = ""
    #s += year+" "+mon+" "+day+" "+hhmm+" "+str(self._shape2D[0])+" "+str(self._shape2D[1])+" "+str(self._shape3D[2])+" "+str(self.p["deltax"])+" "+str(self.p["deltay"])+"\n"
    
    #for xx in range(self._shape2D[0]):
      #for yy in range(self._shape2D[1]):
	#s += str(xx+1)+" "+str(yy+1)+"\n"
	#s += '%3.2f'%self.p["lat"][xx,yy]+" "+'%3.2f'%self.p["lon"][xx,yy]+" "+str(self.p["lfrac"][xx,yy])+" "+str(self.p["wind10u"][xx,yy])+" "+str(self.p["wind10v"][xx,yy])+"\n"
	#s += str(self.p["iwv"][xx,yy])+" "+str(self.p["cwp"][xx,yy])+" "+str(self.p["iwp"][xx,yy])+" "+str(self.p["rwp"][xx,yy])+" "+str(self.p["swp"][xx,yy])+" "+str(self.p["gwp"][xx,yy])
	#if (self.nmlSet["n_moments"]) == 2: s += " "+str(self.p["hwp"][xx,yy])
	#s += "\n"
	#s += '%6.1f'%self.p["hgt_lev"][xx,yy,0]+" "+'%6.1f'%self.p["press_lev"][xx,yy,0]+" "+'%3.2f'%self.p["temp_lev"][xx,yy,0]+" "+'%1.4f'%(self.p["relhum_lev"][xx,yy,0]*100)+"\n"
	#for zz in range(1,self._shape3D[2]+1):
	  #s += '%6.1f'%self.p["hgt_lev"][xx,yy,zz]+" "+'%6.1f'%self.p["press_lev"][xx,yy,zz]+" "+'%3.2f'%self.p["temp_lev"][xx,yy,zz]+" "+'%1.4f'%(self.p["relhum_lev"][xx,yy,zz]*100)+\
	  #" "+'%9e'%self.p["cwc_q"][xx,yy,zz-1]+" "+'%9e'%self.p["iwc_q"][xx,yy,zz-1]+" "+'%9e'%self.p["rwc_q"][xx,yy,zz-1]+" "+'%9e'%self.p["swc_q"][xx,yy,zz-1]+" "+'%9e'%self.p["gwc_q"][xx,yy,zz-1]
	  #if (self.nmlSet["n_moments"]) == 2: s += " "+'%9e'%self.p["hwc_q"][xx,yy,zz-1]
	  #s += "\n"

    
    ## write stuff to file
    #f = open(filename, 'w')
    #f.write(s)
    #f.close()
    return

  def createProfile(self,**kwargs):
    '''
    Function to create Pamtra Profiles.
    
    Variables ending on _lev mean level values, variables without are layer values (height vector one entry shorter!)!
    
    Everything is needed in SI units, relhum is in Pa/PA not %
    
    The following variables are mandatroy:
    hgt_lev, temp_lev, press_lev and (relhum_lev OR q)
    
    The following variables are optional and guessed if not provided:  "timestamp","lat","lon","lfrac","wind10u","wind10v","hgt_lev","cwc_q","iwc_q","rwc_q","swc_q","gwc_q"
    
    The integrated values are skipped if not provided, they have no effect on pamtra!:
    "iwv","cwp","iwp","rwp","swp","gwp"
    
    The following variables are only needed together with the 2 moments scheme!
    "hwc_q","hwp", "cwc_n","iwc_n","rwc_n","swc_n","gwc_n","hwc_n"
    #'''
      
    allVars = self.default_p_vars  
      
    for key in kwargs.keys():
      if key not in allVars:
        raise TypeError("Could not parse "+key)
    
    if not ("hgt_lev" in kwargs.keys() and "temp_lev" in kwargs.keys() and "press_lev" in kwargs.keys() and ("relhum_lev" in kwargs.keys() or "q" in kwargs.keys())):
      raise TypeError("I need hgt_lev and temp_lev and press_lev and (relhum_lev or q)!")
    
    #assert self.df.nhydro > 0
    
    noDims = len(np.shape(kwargs["hgt_lev"]))
    
    if noDims == 1:
      self.p["ngridx"] = 1
      self.p["ngridy"] = 1
    elif noDims == 2:
      self.p["ngridx"] = np.shape(kwargs["hgt_lev"])[0]
      self.p["ngridy"] = 1
    elif noDims == 3:
      self.p["ngridx"] = np.shape(kwargs["hgt_lev"])[0]
      self.p["ngridy"] = np.shape(kwargs["hgt_lev"])[1]
    else:
      print "Too many dimensions!"
      raise IOError
    
    self.p["max_nlyrs"] = np.shape(kwargs["hgt_lev"])[-1] -1
    
    if self.p["max_nlyrs"] > 200:
      warnings.warn("Too many layers for pamtra (max:200): " + str(self.p["max_nlyrs"]),Warning)
    
    #if np.any(self.p["nlyrs"] != self.p["max_nlyrs"]):
      #self._radiosonde = True
    #else:
      #self._radiosonde = False
      
    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin-1)

    self.p["nlyrs"] = np.sum(kwargs["hgt_lev"]!=missingNumber,axis=-1) -1
    self.p["nlyrs"] = self.p["nlyrs"].reshape(self._shape2D)
    
    self.p["hgt_lev"] = kwargs["hgt_lev"].reshape(self._shape3Dplus)
    self.p["temp_lev"] = kwargs["temp_lev"].reshape(self._shape3Dplus)
    self.p["press_lev"] = kwargs["press_lev"].reshape(self._shape3Dplus)


    self.p["model_i"] = np.array(np.where(self.p["press_lev"][:,:,0])[0]).reshape(self._shape2D) +1
    self.p["model_j"] = np.array(np.where(self.p["press_lev"][:,:,0])[1]).reshape(self._shape2D) +1

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
      
    for environment, preset in [["lat",50.938056],["lon",6.956944],["lfrac",1],["wind10u",0],["wind10v",0]]:
      if environment not in kwargs.keys():
        self.p[environment] = np.ones(self._shape2D)*preset
        warnings.warn("%s set to %s"%(environment,preset,), Warning)
      else:
        if type(kwargs[environment]) in (int,np.int32,np.int64,float,np.float32,np.float64):
          self.p[environment] = np.ones(self._shape2D) * kwargs[environment]
        else:
          self.p[environment] = kwargs[environment].reshape(self._shape2D)
    
    for qValue in ["hydro_q","hydro_reff","hydro_n"]:
      if qValue not in kwargs.keys():
        self.p[qValue] = np.zeros(self._shape4D)
        warnings.warn(qValue + " set to 0", Warning)
      else:
        self.p[qValue] = kwargs[qValue].reshape(self._shape4D)

    for qValue in ["airturb"]:
      if qValue not in kwargs.keys():
        self.p[qValue] = np.zeros(self._shape3D)
        warnings.warn(qValue + " set to 0", Warning)
      else:
        self.p[qValue] = kwargs[qValue].reshape(self._shape3D)
        
    if "obs_height" in kwargs.keys():
      self.p["obs_height"] = kwargs["obs_height"].reshape(self._shape2D)

    if "relhum_lev" in kwargs.keys():
      self.p["relhum_lev"] = kwargs["relhum_lev"].reshape(self._shape3Dplus)
    else:
      self._calcRelhum_lev()

    #clean up: remove all nans
    for key in self.p.keys():
      #apply only to arrays, lists or tupels
      if np.sum(np.shape(self.p[key])) > 0:
        self.p[key][np.isnan(self.p[key])] = missingNumber
    
    return


  def _checkData(self):
    maxLimits = {"relhum_lev":200,"hydro_q":0.05,"temp_lev":320,"press_lev":110000}
    minLimits = {"relhum_lev":0,"hydro_q": 0 ,"temp_lev":170,"press_lev":1}
    for key in maxLimits.keys():
      if np.max(self.p[key]) > maxLimits[key]:
        raise ValueError("unrealistic value for "+ key + ": " +str(np.max(self.p[key])) + ", maximum is " + str(maxLimits[key]))
    for key in minLimits.keys():
      if len(self.p[key][self.p[key] != missingNumber]) > 0:
        if np.min(self.p[key][self.p[key] != missingNumber]) < minLimits[key]:
          raise ValueError("unrealistic value for "+ key + ": " +str(np.min(self.p[key])) + ", minimum is " + str(minLimits[key]))

  def filterProfiles(self,condition):
    '''
    discard profiles, which do not fullfill "condition" (2D boolean array)
    
    Note: If the initial field is a 2D grid, the field is flattend after application
    
    Applicate between CreateProfile and runPamtra
    
    '''
    if condition.shape != self._shape2D and condition.shape != (self._shape2D[0]*self._shape2D[1],):
      raise ValueError("shape mismatch, shape of condition must be 2D field!")
    
    condition = condition.reshape(self._shape2D)
    
    #create a new shape!
    self.p["ngridx"] = np.sum(condition)
    self.p["ngridy"] = 1
    
    try: self.df.fs_nbin =  self.df.dataFullSpec["d_bound_ds"].shape[-1]
    except: self.df.fs_nbin = 0
    
    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin-1)
        
    for key in ["unixtime","nlyrs","lat","lon","lfrac","model_i","model_j","wind10u","wind10v","obs_height"]:
      if key in self.p.keys(): self.p[key] = self.p[key][condition].reshape(self._shape2D)
      
    for key in ["hydro_q","hydro_n","hydro_reff"]:
      if key in self.p.keys(): self.p[key] = self.p[key][condition].reshape(self._shape4D)
      
    for key in ["hgt_lev","temp_lev","press_lev","relhum_lev"]:
      if key in self.p.keys(): self.p[key] = self.p[key][condition].reshape(self._shape3Dplus)
    
    for key in ["airturb",'temp', 'press', 'relhum','hgt']:
      if key in self.p.keys(): self.p[key] = self.p[key][condition].reshape(self._shape3D)    
 
    for key in self.df.data4D.keys():
      self.df.data4D[key] = self.df.data4D[key][condition].reshape(self._shape4D)
      
    for key in self.df.dataFullSpec.keys():
      if key == "delta_d_ds": shape5D = self._shape5Dplus
      else: shape5D = self._shape5D
      self.df.dataFullSpec[key] = self.df.dataFullSpec[key][condition].reshape(shape5D)
    return

  def rescaleHeights(self,new_hgt_lev):
  
    # sort height vectors
    old_hgt_lev = self.p["hgt_lev"]
    new_hgt = (new_hgt_lev[...,1:] + new_hgt_lev[...,:-1])/2.
    old_hgt = (old_hgt_lev[...,1:] + old_hgt_lev[...,:-1])/2.
    
    self.p["max_nlyrs"] = len(new_hgt_lev) -1
    
    #new shape
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin-1)
    assert len(df.data4D.keys() == 0) 
    assert len(df.dataFullSpec.keys() == 0) 

    for key in ["hydro_q","hydro_n","hydro_reff",]:
      #make new array
      newP = np.ones(self._shape3D) * missingNumber
      for x in xrange(self._shape2D[0]):
        for y in xrange(self._shape2D[1]):
          for h in xrange(self.df.nhydro):
          #interpolate!
            newP[x,y,:,h] = np.interp(new_hgt,old_hgt[x,y,:,h],self.p[key][x,y,:,h])
      #save new array
      self.p[key] = newP
      #and mark all entries below -1 as missing Number!
      self.p[key][self.p[key]<-1] = missingNumber
      
    for key in ["hgt_lev","temp_lev","relhum_lev"]:
      newP = np.ones(self._shape3Dplus) * missingNumber
      for x in xrange(self._shape2D[0]):
        for y in xrange(self._shape2D[1]):
          newP[x,y] = np.interp(new_hgt_lev,old_hgt_lev[x,y],self.p[key][x,y])
      self.p[key] = newP
      self.p[key][self.p[key]<-1] = missingNumber
      
    for key in ["airturb"]:
      newP = np.ones(self._shape3D) * missingNumber
      for x in xrange(self._shape2D[0]):
        for y in xrange(self._shape2D[1]):
          newP[x,y] = np.interp(new_hgt,old_hgt[x,y],self.p[key][x,y])
      self.p[key] = newP
      self.p[key][self.p[key]<-1] = missingNumber
      
      
    for key in ["press_lev"]:
      newP = np.ones(self._shape3Dplus) * missingNumber
      for x in xrange(self._shape2D[0]):
        for y in xrange(self._shape2D[1]):
          newP[x,y] = np.exp(np.interp(new_hgt_lev,old_hgt_lev[x,y],np.log(self.p[key][x,y])))
      self.p[key] = newP
      self.p[key][self.p[key]<-1] = missingNumber

    self.p["nlyrs"] = np.sum(self.p["hgt_lev"] != missingNumber,axis=-1) -1
    
    if self.p["max_nlyrs"] > 200:
      warnings.warn("Still too many layers for pamtra (max:200): " + str(self.p["max_nlyrs"]),Warning)


    return
    
  def addIntegratedValues(self):
    raise NotImplementedError("not yet avaiable in pamtra v 1.0")
    #for pDict,qValue,intValue in [[self._helperP,"q","iwv"],[self.p,"cwc_q","cwp"],[self.p,"iwc_q","iwp"],[self.p,"rwc_q","rwp"],[self.p,"swc_q","swp"],[self.p,"gwc_q","gwp"],[self.p,"hwc_q","hwp"]]:
      ##now we need q!
      #self._calcQ()
      ##nothing to do without hydrometeors:
      #if np.all(pDict[qValue]==0):
        #self.p[intValue] = np.zeros(self._shape2D)
      #else:
        #self._calcMoistRho() #provides also temp,press,dz and q!
        #self.p[intValue] =  np.sum(pDict[qValue]*self._helperP["rho_moist"]*self._helperP["dz"],axis=-1)
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
        i_top, i_base, i_cloud =  meteoSI.detect_liq_cloud(self.p["hgt_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["temp_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["relhum_lev"][x,y,0:self.p["nlyrs"][x,y]])
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
        i_top, i_base, i_cloud =  meteoSI.detect_liq_cloud(self.p["hgt_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["temp_lev"][x,y,0:self.p["nlyrs"][x,y]], self.p["relhum_lev"][x,y,0:self.p["nlyrs"][x,y]])
        for k in np.arange(len(i_top)):
          i_cloud_k = np.arange(i_base[k],i_top[k]+1) 
          self.p["cwc_q"][x,y,i_cloud_k[:-1]] = meteoSI.mod_ad(self.p["temp_lev"][x,y,i_cloud_k], self.p["press_lev"][x,y,i_cloud_k], self.p["hgt_lev"][x,y,i_cloud_k], 1)

    self.p["cwp"][:] = 0
    
    self._calcMoistRho() #provides also temp,press,dz and q!
    self.p["cwp"][:] =  np.sum(self.p["cwc_q"]*self._helperP["rho_moist"]*self._helperP["dz"],axis=-1)

    return

  #def addPIA(self)
  #todo: make pia for all hydrometeortypes!
    #assert self.r["nmlSettings"]["output"]["activeLogScale"]
    
    #Att_atmo = pam.r["Att_atmo"]
    #Att_atmo[Att_atmo==missingNumber] = 0
    #Att_hydro = pam.r["Att_hydro"]
    #Att_hydro[Att_hydro==missingNumber] = 0    
    #self.r["PIA"] = np.zeros(np.shape(Att_atmo))
    #for hh in range(self.p["max_nlyrs"]):
      #self.r["PIA"][:,:,hh,:] = Att_atmo[:,:,hh,:] +Att_hydro[:,:,hh,:] + 2*np.sum(Att_atmo[:,:,:hh,:] +Att_hydro[:,:,:hh,:],axis=2)#only half for the current layer
    #return

    
  def createFullProfile(self,timestamp,lat,lon,lfrac,wind10u,wind10v,
      obs_height,
      hgt_lev,press_lev,temp_lev,relhum_lev,
      hydro_q,hydro_n,hydro_reff):
    
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
    
    if self.p["max_nlyrs"] > 200:
      warnings.warn("Too many layers for pamtra (max:200): " + str(self.p["max_nlyrs"]),Warning)
    
    self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
    self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],)
    self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"]+1,)
    self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
    self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)
    self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin-1)
    self.p["unixtime"] = timestamp.reshape(self._shape2D)
        
    self.p["lat"] = lat.reshape(self._shape2D)
    self.p["lon"] = lon.reshape(self._shape2D)
    self.p["lfrac"] = lfrac.reshape(self._shape2D)
    self.p["model_i"] = np.array(np.where(lat.reshape(self._shape2D))[0]).reshape(self._shape2D) +1
    self.p["model_j"] = np.array(np.where(lon.reshape(self._shape2D))[1]).reshape(self._shape2D) +1
    self.p["wind10u"] = wind10u.reshape(self._shape2D)
    self.p["wind10v"] = wind10v.reshape(self._shape2D)
    self.p["obs_height"] = obs_height.reshape(self._shape2D)

    self.p["hydro_q"] = hydro_q.reshape(self._shape4D)
    self.p["hydro_n"] = hydro_n.reshape(self._shape4D)
    self.p["hydro_reff"] = hydro_reff.reshape(self._shape4D)
    
    self.p["hgt_lev"] = hgt_lev.reshape(self._shape3Dplus)
    self.p["temp_lev"] = temp_lev.reshape(self._shape3Dplus)
    self.p["press_lev"] = press_lev.reshape(self._shape3Dplus)
    self.p["relhum_lev"] = relhum_lev.reshape(self._shape3Dplus)    
    return   
    
  def runPamtra(self,freqs,checkData=True):
    '''
    run Pamtra from python
    '''
    tttt = time.time()

    
    if type(freqs) in (int,np.int32,np.int64,float,np.float32,np.float64): freqs = [freqs]
    
    self.set["freqs"] = freqs
    self.set["nfreqs"] = len(freqs)
    assert self.set["nfreqs"] > 0

    if checkData: self._checkData()

    fortResults, fortObject = pyPamtraLibWrapper.PamtraFortranWrapper(self.set,self.nmlSet,self.df,self.p)
    self.r = fortResults
    self.fortObject = fortObject
    
     
    self.r["nmlSettings"] = self.nmlSet

    if self.set["pyVerbose"] > 0: print "pyPamtra runtime:", time.time() - tttt
    return

    
  def runParallelPamtra(self,freqs,pp_local_workers="auto",pp_deltaF=1,pp_deltaX=0,pp_deltaY = 0, activeFreqs="auto", passiveFreqs="auto",checkData=True,timeout=86400):
    '''
    run Pamtra from python
    '''
    if pp_deltaF==0: pp_deltaF = self.set["nfreqs"]
    if pp_deltaX==0: pp_deltaX = self.p["ngridx"]
    if pp_deltaY==0: pp_deltaY = self.p["ngridy"]
    
    if hasattr(self, "fortObject"): del self.fortObject
    
    if pp_local_workers == "auto": pp_local_workers = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=pp_local_workers)
    tttt = time.time()

    
    if type(freqs) in (int,np.int32,np.int64,float,np.float32,np.float64): freqs = [freqs]
    
    self.set["freqs"] = freqs
    self.set["nfreqs"] = len(freqs)
    assert self.set["nfreqs"] > 0

    if checkData: self._checkData()

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
      
            print "submitting job ", pp_i, pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY
      
            indices = [pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY]       
            profilePart, dfPart, settings = self._sliceProfile(*indices)

            jobs.append(pool.apply_async(pyPamtraLibWrapper.parallelPamtraFortranWrapper,(
            indices,
            settings,self.nmlSet,dfPart,profilePart),{"returnModule":False}))#),callback=self.pp_resultData.append)
           
            
            pp_i += 1
            print "submitted job: ", pp_i

            
      #pool.close()
      #pool.join()
    except KeyboardInterrupt:
      pool.terminate()
      pool.join()
      print "TERMINATED: KeyboardInterrupt"

    print "waiting for all jobs to finish"
    for jj,job in enumerate(jobs):
      self._joinResults(job.get(timeout=timeout))
      print "got job", jj+1

    self.r["nmlSettings"] = self.nmlSet

    pool.terminate()
    if self.set["pyVerbose"] > 0: print "pyPamtra runtime:", time.time() - tttt
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
       
    dfData = deepcopy(self.df)   
    for key in dfData.data4D.keys():
      dfData.data4D[key] = self.df.data4D[key][pp_startX:pp_endX,pp_startY:pp_endY]
       
    for key in dfData.dataFullSpec.keys():
      dfData.dataFullSpec[key] = self.df.dataFullSpec[key][pp_startX:pp_endX,pp_startY:pp_endY]
       
    settings = deepcopy(self.set)
    settings["nfreqs"] = pp_endF - pp_startF
    settings["freqs"] = self.set["freqs"][pp_startF:pp_endF]
       
    return profilePart, dfData, settings
    
  def _prepareResults(self):
    
    try: maxNBin1 = np.max(self.df.data["nbin"]) + 1
    except: 
      try: 
        maxNBin1 = np.max(self.df.data4D["nbin"]) + 1
      except:
        maxNBin1 = self.df.dataFullSpec["d_bound_ds"].shape[-1]
    radar_spectrum_length = self.nmlSet["radar_nfft"]
    
    
    self.r = dict()
    self.r["Ze"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
    self.r["Att_hydro"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
    self.r["Att_atmo"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
    self.r["radar_hgt"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],))*missingNumber
    if self.nmlSet["radar_mode"]=="spectrum":
      self.r["radar_spectra"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],radar_spectrum_length,))*missingNumber
    else:
      self.r["radar_spectra"] = np.array([missingNumber])
    self.r["radar_snr"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],))*missingNumber
    self.r["radar_moments"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],4,))*missingNumber
    self.r["radar_slopes"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],2,))*missingNumber
    self.r["radar_edges"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],2,))*missingNumber
    self.r["radar_quality"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.set["nfreqs"],),dtype=int)*missingNumber
    self.r["tb"] = np.ones((self.p["ngridx"],self.p["ngridy"],self._noutlevels,self._nangles,self.set["nfreqs"],self._nstokes))*missingNumber
    if self.nmlSet["save_psd"]:
      self.r["psd_area"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,maxNBin1))*missingNumber
      self.r["psd_f"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,maxNBin1))*missingNumber
      self.r["psd_d_bound"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,maxNBin1))*missingNumber
      self.r["psd_mass"] = np.ones((self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,maxNBin1))*missingNumber
    else: #save memory
      self.r["psd_area"] = np.array([missingNumber])
      self.r["psd_f"] = np.array([missingNumber])
      self.r["psd_d_bound"] = np.array([missingNumber])
      self.r["psd_mass"] =np.array([missingNumber])

    
    return 
    
  def _joinResults(self,resultList):
    '''
    Collect the data of parallel pyPamtra    
    '''
    #if pp_ii == 1: import pdb;pdb.set_trace()
    #logging.debug(str(pp_ii)+": callback started ")
    #print "toll", results[0][0][1]
    #print "toller", np.shape(self.r["radar_snr"][pp_startX:pp_endX,pp_startY:pp_endY,:,pp_startF:pp_endF])    
    if self.set["pyVerbose"] > 2: print "Callback started:", 
    
    [pp_startF,pp_endF,pp_startX,pp_endX,pp_startY,pp_endY], results = resultList
       
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
      for key in ["psd_d_bound","psd_f","psd_mass","psd_area"]:
        self.r[key][pp_startX:pp_endX,pp_startY:pp_endY] = results[key]
    
    return
    
    
    
  def writeResultsToNumpy(self,fname,seperateFiles=False):
    
    if not seperateFiles:
      '''
      write the complete state of the session (profile,results,settings to a file
      '''
      f = open(fname, "w")
      pickle.dump([self.r,self.p,self.nmlSet,self.set,self.df.data,self.df.data4D,df.dataFullSpec], f)
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
	[self.r,self.p,self.nmlSet,self.set,self.df.data,self.df.data4D,df.dataFullSpec] = pickle.load(f)
	f.close()
      except:
	print formatExceptionInfo()	
	raise IOError ("Could not read data from file")
      
      self.df.nhydro = len(self.df.data)
      try: self.df.fs_nbin =  self.df.dataFullSpec["d_bound_ds"].shape[-1]
      except: self.df.fs_nbin = 0
      self._shape2D = (self.p["ngridx"],self.p["ngridy"],)
      self._shape3D = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"],)
      self._shape3Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["nlyrs"]+1,)
      self._shape4D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro)
      self._shape5D = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin)    
      self._shape5Dplus = (self.p["ngridx"],self.p["ngridy"],self.p["max_nlyrs"],self.df.nhydro,self.df.fs_nbin-1)
      
      return
  

  
  
  def writeResultsToNetCDF(self,fname,profileVars="all",ncForm="NETCDF3_CLASSIC"):
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
      
    try: 
      self.r
      self.r["pamtraVersion"]
    except:
      raise IOError ("run runPamtra first!")

    if pyNc: cdfFile = nc.Dataset(fname,"w",ncForm)
    else: cdfFile = nc.NetCDFFile(fname,"w")
    
    #write meta data
    cdfFile.history = "Created with pyPamtra (Version: "+self.r["pamtraVersion"]+", Git Hash: "+self.r["pamtraHash"]+") by "+self.nmlSet["creator"]+" (University of Cologne, IGMK) at " + datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    
    cdfFile.properties = str(self.r["nmlSettings"])
    #make frequsnions
    cdfFile.createDimension('grid_x',int(self.p["ngridx"]))
    cdfFile.createDimension('grid_y',int(self.p["ngridy"]))
    cdfFile.createDimension('frequency',int(self.set["nfreqs"]))
    
    if (self.r["nmlSettings"]["run_mode"]["passive"]):
      cdfFile.createDimension('angles',len(self.r["angles_deg"]))
      cdfFile.createDimension('outlevels',int(self._noutlevels))
      cdfFile.createDimension('stokes',int(self._nstokes))
      
    if (self.r["nmlSettings"]["run_mode"]["radar_mode"] in ["spectrum","moments"]):
      cdfFile.createDimension('nfft',int(self.r["nmlSettings"]["radar_simulator"]["radar_nfft"])) 
      
    if (self.r["nmlSettings"]["run_mode"]["active"]):
      cdfFile.createDimension('heightbins',int(self.p["max_nlyrs"]))
    
    dim2d = ("grid_x","grid_y",)
    dim3d = ("grid_x","grid_y","heightbins",)
    dim4d = ("grid_x","grid_y","heightbins","frequency")
    dim5d = ("grid_x","grid_y","heightbins","frequency","nfft")
    dim6d = ("grid_x","grid_y","outlevels","angles","frequency","stokes")
      
      
      

    attUnit = "dBz"
    zeUnit = "dBz"
 
    #create and write dim variables
    
    fillVDict = dict()
    #little cheat to avoid hundreds of if, else...
    if pyNc: fillVDict["fill_value"] = missingNumber
    
    nc_frequency = cdfFile.createVariable('frequency','f',('frequency',),**fillVDict)
    nc_frequency.units = 'GHz'
    nc_frequency[:] = self.set["freqs"]
    if not pyNc: nc_frequency._FillValue =missingNumber
    
    nc_gridx = cdfFile.createVariable('grid_x','f',('grid_x',),**fillVDict)#= missingNumber)
    nc_gridx.units = '-'
    nc_gridx[:] = np.arange(self.p["ngridx"],dtype="f")
    if not pyNc: nc_gridx._FillValue =missingNumber
    
    nc_gridy = cdfFile.createVariable('grid_y','f',('grid_y',),**fillVDict)
    nc_gridy.units = '-'
    nc_gridy[:] = np.arange(self.p["ngridy"],dtype="f")
    if not pyNc: nc_gridy._FillValue =missingNumber
    
    
    
    if (self.r["nmlSettings"]["run_mode"]["active"]):

      nc_heightbins = cdfFile.createVariable('heightbins', 'i',("heightbins",),**fillVDict)
      nc_heightbins.units = "-"
      nc_heightbins[:] = np.arange(0,self.p["max_nlyrs"],dtype="i")
      if not pyNc: nc_heightbins._FillValue =int(missingNumber)

      nc_height = cdfFile.createVariable('height', 'f',dim3d,**fillVDict)
      nc_height.units = "m"
      nc_height[:] = np.array(self.r["radar_hgt"],dtype="f")
      if not pyNc: nc_height._FillValue =missingNumber
      
    
    if (self.r["nmlSettings"]["run_mode"]["passive"]):
      nc_angle = cdfFile.createVariable('angles','f',('angles',),**fillVDict)
      nc_angle.units = 'deg'
      nc_angle[:] = np.array(self.r["angles"],dtype="f")
      if not pyNc: nc_angle._FillValue =missingNumber
      
      nc_stokes = cdfFile.createVariable('stokes', 'c',("stokes",),**fillVDict)
      nc_stokes.units = "-"
      nc_stokes[:] = "HV"
      
      nc_out = cdfFile.createVariable('outlevels', 'f',("outlevels",),**fillVDict)
      nc_out.units = "m over sea level (top of atmosphere value) OR m over ground (ground value)"
      nc_out[:] = np.array([self.r["nmlSettings"]["output"]["obs_height"],0],dtype="f")
      if not pyNc: nc_out._FillValue =missingNumber
      
    #create and write variables
    
    nc_model_i = cdfFile.createVariable('model_i', 'i',dim2d,**fillVDict)
    nc_model_i.units = "-"
    nc_model_i[:] = np.array(self.p["model_i"],dtype="i")
    if not pyNc: nc_model_i._FillValue =int(missingNumber)
    
    nc_model_j = cdfFile.createVariable('model_j', 'i',dim2d,**fillVDict)
    nc_model_j.units = "-"
    nc_model_j[:] = np.array(self.p["model_j"],dtype="i")
    if not pyNc: nc_model_j._FillValue =int(missingNumber)
      
    nc_nlyrs = cdfFile.createVariable('nlyr', 'i',dim2d,**fillVDict)
    nc_nlyrs.units = "-"
    nc_nlyrs[:] = np.array(self.p["nlyrs"],dtype="i")
    if not pyNc: nc_nlyrs._FillValue =int(missingNumber)
    
    nc_time = cdfFile.createVariable('datatime', 'i',dim2d,**fillVDict)
    nc_time.units = "seconds since 1970-01-01 00:00:00"
    nc_time[:] = np.array(self.p["unixtime"],dtype="i")
    if not pyNc: nc_time._FillValue =int(missingNumber)
    
    nc_longitude = cdfFile.createVariable('longitude', 'f',dim2d,**fillVDict)
    nc_longitude.units = "deg.dec"
    nc_longitude[:] = np.array(self.p["lon"],dtype="f")
    if not pyNc: nc_longitude._FillValue =missingNumber
    
    nc_latitude = cdfFile.createVariable('latitude', 'f',dim2d,**fillVDict)
    nc_latitude.units = "deg.dec"
    nc_latitude[:] = np.array(self.p["lat"],dtype="f")
    if not pyNc: nc_latitude._FillValue =missingNumber
      
    nc_lfrac = cdfFile.createVariable('lfrac', 'f',dim2d,**fillVDict)
    nc_lfrac.units = "-"
    nc_lfrac[:] = np.array(self.p["lfrac"],dtype="f")
    if not pyNc: nc_lfrac._FillValue =missingNumber
    
    
    if (self.r["nmlSettings"]["run_mode"]["active"]):
         
      nc_Ze = cdfFile.createVariable('Ze', 'f',dim4d,**fillVDict)
      nc_Ze.units = zeUnit
      nc_Ze[:] = np.array(self.r["Ze"],dtype='f')
      if not pyNc: nc_Ze._FillValue =missingNumber

      nc_Attenuation_Hydrometeors = cdfFile.createVariable('Attenuation_Hydrometeors', 'f',dim4d,**fillVDict)
      nc_Attenuation_Hydrometeors.units = attUnit
      nc_Attenuation_Hydrometeors[:] = np.array(self.r["Att_hydro"],dtype='f')
      if not pyNc: nc_Attenuation_Hydrometeors._FillValue =missingNumber
        
      nc_Attenuation_Atmosphere = cdfFile.createVariable('Attenuation_Atmosphere', 'f',dim4d,**fillVDict)
      nc_Attenuation_Atmosphere.units = attUnit
      nc_Attenuation_Atmosphere[:] = np.array(self.r["Att_atmo"],dtype='f')
      if not pyNc: nc_Attenuation_Atmosphere._FillValue =missingNumber
      
      if ((self.r["nmlSettings"]["run_mode"]["radar_mode"] == "spectrum") or (self.r["nmlSettings"]["run_mode"]["radar_mode"] == "moments")):
	nc_snr=cdfFile.createVariable('Radar_SNR', 'f',dim4d,**fillVDict)
	nc_snr.units="dB"
	nc_snr[:] = np.array(self.r["radar_snr"],dtype='f')
	if not pyNc: nc_snr._FillValue =missingNumber

	nc_fvel=cdfFile.createVariable('Radar_FallVelocity', 'f',dim4d,**fillVDict)
	nc_fvel.units="m/s"
	nc_fvel[:] = np.array(self.r["radar_moments"][...,0],dtype='f')
	if not pyNc: nc_fvel._FillValue =missingNumber
	
	nc_specw=cdfFile.createVariable('Radar_SpectralWidth', 'f',dim4d,**fillVDict)
	nc_specw.units="m/s"
	nc_specw[:] = np.array(self.r["radar_moments"][...,1],dtype='f')
	if not pyNc: nc_specw._FillValue =missingNumber
	
	nc_skew=cdfFile.createVariable('Radar_Skewness', 'f',dim4d,**fillVDict)
	nc_skew.units="-"
	nc_skew[:] = np.array(self.r["radar_moments"][...,2],dtype='f')
	if not pyNc: nc_skew._FillValue =missingNumber
	
	nc_kurt=cdfFile.createVariable('Radar_Kurtosis', 'f',dim4d,**fillVDict)
	nc_kurt.units="-"
	nc_kurt[:] = np.array(self.r["radar_moments"][...,3],dtype='f')
	if not pyNc: nc_kurt._FillValue =missingNumber
	
	nc_lslop=cdfFile.createVariable('Radar_LeftSlope', 'f',dim4d,**fillVDict)
	nc_lslop.units="dB/(m/s)"
	nc_lslop[:] = np.array(self.r["radar_slopes"][...,0],dtype='f')
	if not pyNc: nc_lslop._FillValue =missingNumber
	
	nc_rslop=cdfFile.createVariable('Radar_RightSlope', 'f',dim4d,**fillVDict)
	nc_rslop.units="dB/(m/s)"
	nc_rslop[:] = np.array(self.r["radar_slopes"][...,1],dtype='f')
	if not pyNc: nc_rslop._FillValue =missingNumber

        nc_lslop=cdfFile.createVariable('Radar_LeftEdge', 'f',dim4d,**fillVDict)
        nc_lslop.units="m/s"
        nc_lslop[:] = np.array(self.r["radar_edges"][...,0],dtype='f')
        if not pyNc: nc_lslop._FillValue =missingNumber
        
        nc_rslop=cdfFile.createVariable('Radar_RightEdge', 'f',dim4d,**fillVDict)
        nc_rslop.units="m/s"
        nc_rslop[:] = np.array(self.r["radar_edges"][...,1],dtype='f')
        if not pyNc: nc_rslop._FillValue =missingNumber	
	
	nc_qual=cdfFile.createVariable('Radar_Quality', 'i',dim4d,**fillVDict)
	nc_qual.units="bytes"
	nc_qual.description="1st byte: aliasing; 2nd byte: 2nd peak present; 7th: no peak found"
	nc_qual[:] = np.array(self.r["radar_quality"],dtype='i')
	if not pyNc: nc_qual._FillValue =missingNumber
	
	
	if ((self.r["nmlSettings"]["run_mode"]["radar_mode"] == "spectrum")): 

	  nc_vel=cdfFile.createVariable('Radar_Velocity', 'f',("nfft",),**fillVDict)
	  nc_vel.units="m/s"
	  nc_vel[:] = np.array(self.r["radar_vel"],dtype='f')
	  if not pyNc: nc_vel._FillValue =missingNumber
	  
	  nc_spec=cdfFile.createVariable('Radar_Spectrum', 'f',dim5d,**fillVDict)
	  nc_spec.units="dBz"
	  nc_spec[:] = np.array(self.r["radar_spectra"],dtype='f')
	  if not pyNc: nc_spec._FillValue =missingNumber
	     
    if (self.r["nmlSettings"]["run_mode"]["passive"]):
      nc_tb = cdfFile.createVariable('tb', 'f',dim6d,**fillVDict)
      nc_tb.units = "K"
      nc_tb[:] = np.array(self.r["tb"],dtype='f')
      if not pyNc: nc_tb._FillValue =missingNumber
        
    #profile data
    
    if ("iwv" in profileVars or profileVars =="all") and ("iwv" in self.p.keys()):
      nc_iwv = cdfFile.createVariable('iwv', 'f',dim2d,**fillVDict)
      nc_iwv.units = "kg/m^2"
      nc_iwv[:] = np.array(self.p["iwv"],dtype="f")
      if not pyNc: nc_iwv._FillValue =missingNumber

    if ("cwp" in profileVars or profileVars =="all") and ("cwp" in self.p.keys()):
      nc_cwp= cdfFile.createVariable('cwp', 'f',dim2d,**fillVDict)
      nc_cwp.units = "kg/m^2"
      nc_cwp[:] = np.array(self.p["cwp"],dtype="f")
      if not pyNc: nc_cwp._FillValue =missingNumber
  
    if ("iwp" in profileVars or profileVars =="all") and ("iwp" in self.p.keys()):
      nc_iwp = cdfFile.createVariable('iwp', 'f',dim2d,**fillVDict)
      nc_iwp.units = "kg/m^2"
      nc_iwp[:] = np.array(self.p["iwp"],dtype="f")
      if not pyNc: nc_iwp._FillValue =missingNumber

    if ("rwp" in profileVars or profileVars =="all") and ("rwp" in self.p.keys()):
      nc_rwp = cdfFile.createVariable('rwp', 'f',dim2d,**fillVDict)
      nc_rwp.units = "kg/m^2"
      nc_rwp[:] = np.array(self.p["rwp"],dtype="f")
      if not pyNc: nc_rwp._FillValue =missingNumber

    if ("swp" in profileVars or profileVars =="all") and ("swp" in self.p.keys()):
      nc_swp = cdfFile.createVariable('swp', 'f',dim2d,**fillVDict)
      nc_swp.units = "kg/m^2"
      nc_swp[:] = np.array(self.p["swp"],dtype="f")
      if not pyNc: nc_swp._FillValue =missingNumber

    if ("gwp" in profileVars or profileVars =="all") and ("gwp" in self.p.keys()):
      nc_gwp = cdfFile.createVariable('gwp', 'f',dim2d,**fillVDict)
      nc_gwp.units = "kg/m^2"
      nc_gwp[:] = np.array(self.p["gwp"],dtype="f")
      if not pyNc: nc_gwp._FillValue =missingNumber

    if ("hwp" in profileVars or profileVars =="all") and ("hwp" in self.p.keys()):
      nc_hwp = cdfFile.createVariable('hwp', 'f',dim2d,**fillVDict)
      nc_hwp.units = "kg/m^2"
      nc_hwp[:] = np.array(self.p["hwp"],dtype="f")
      if not pyNc: nc_hwp._FillValue =missingNumber

    if ("cloudBase" in profileVars or profileVars =="all") and ("cloudBase" in self.p.keys()):
      nc_cb = cdfFile.createVariable('cloudBase', 'f',dim2d,**fillVDict)
      nc_cb.units = "m"
      nc_cb[:] = np.array(self.p["cloudBase"],dtype="f")
      if not pyNc: nc_cb._FillValue =missingNumber

    if ("cloudTop" in profileVars or profileVars =="all") and ("cloudTop" in self.p.keys()):
      nc_ct = cdfFile.createVariable('cloudTop', 'f',dim2d,**fillVDict)
      nc_ct.units = "m"
      nc_ct[:] = np.array(self.p["cloudTop"],dtype="f")
      if not pyNc: nc_ct._FillValue =missingNumber

    cdfFile.close()
    if self.set["pyVerbose"] > 0: print fname,"written"

    
#some tools

def formatExceptionInfo(maxTBlevel=5):
  cla, exc, trbk = sys.exc_info()
  excName = cla.__name__
  try:
    excArgs = exc.__dict__["args"]
  except KeyError:
    excArgs = "<no args>"
  excTb = traceback.format_tb(trbk, maxTBlevel)
  return (excName, excArgs, excTb)

    
    