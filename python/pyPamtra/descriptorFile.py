# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import csv

missingNumber=-9999.
missingIntNumber=int(missingNumber)

class pamDescriptorFile(object):
  '''
  class with the descriptor file content in data. In case you want to use 4D data, use data4D instead, 
  the coreesponding column in data is automatically removed. Inializes with an empty array.
  '''
         
  def __init__(self,parent):
    self.names = np.array(["hydro_name", "as_ratio", "liq_ice", "rho_ms", "a_ms", "b_ms", "alpha_as", "beta_as", "moment_in", "nbin", "dist_name", "p_1", "p_2", "p_3", "p_4", "d_1", "d_2", "scat_name", "vel_size_mod","canting"])
    self.types = ["S15",float,int,float,float,float,float,float,int,int,"S15",float,float,float,float,float,float, "S30", "S30",float]  
    self.data = np.recarray((0,),dtype=zip(self.names, self.types))
    self.data4D = dict()
    self.nhydro = 0
    self.parent = parent
    self.dataFullSpec = dict()
    self.fs_nbin = 0
    return
    
   
  def readFile(self,fname):
    '''
    Read ASCII Pamtra descriptor file

    Parameters
    ----------
    fname : str
      filename
    '''
    f = open(fname,"r")
    for row in csv.reader(f,delimiter=" ",skipinitialspace=True):
      #skipp comments
      if row[0][0] == "!":
        continue
      #make sure line is complete
      assert len(row) == len(self.names)  

      for ii, item in enumerate(row):
      #python does not like double types
        if self.types[ii] == float:
          row[ii] = item.replace("d", "e")
      #fortran strips quotes
        if self.types[ii] == "S15":
          row[ii] = item.replace("'", "").replace('"', '')
  
      self.addHydrometeor(row)
      #print ', '.join(row), len(row)
      
    f.close()  
    #make rec array

  def writeFile(self,fname):
    '''
    Write ASCII Pamtra descriptor file

    Parameters
    ----------
    fname : str
      filename
    '''
    assert len(data[0]) == len(self.names)
    return mlab.rec2csv(data, fname, delimiter=' ',withheader=False)
      
  def addHydrometeor(self,hydroTuple):
    '''
    Add hydrometeor to df array. 

    Parameters
    ----------
    hydroTuple : tuple
      tuple describing the hydrometeor. Consists of
      "hydro_name", "as_ratio", "liq_ice", "rho_ms", "a_ms", "b_ms", "alpha_as", "beta_as", "moment_in", "nbin", "dist_name", "p_1", 
      "p_2", "p_3", "p_4", "d_1", "d_2", "scat_name", "vel_size_mod","canting"
    '''
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
    '''
    Remove hydrometeor from df array. 

    Parameters
    ----------
    hydroName : str
      Name of hydrometeor to be removed.
    '''

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
    '''
    Add 4D information for one hydrometeor property. I.e the 1D value with dimesnions (hydrometeors) in the df array is replaced by an array with (x,y,z,hydrometeor) dimensions. 

    Parameters
    ----------
    key : str
      Name of hydrometeor property to be altered.
    arr : array 
      New 4D values
    '''
    self.data = mlab.rec_drop_fields(self.data,[key])
    self.data4D[key] = arr.reshape(self.parent._shape4D)
    
  def remove4D(self, key, val):
    '''
    Remove 4D information for one hydrometeor property.

    Parameters
    ----------
    key : str
      Name of hydrometeor property to be altered.
    arr : array 
      New 1 D values.
    '''

    self.data = mlab.rec_append_fields(self.data,key,val)
    del data4D[key]
    
  def addFullSpectra(self):
    '''
    This function creates empty arrays which have 5 dimensions (x,y,z,hydrometeor,bin). 
    They are not returned but added to the df object. 
    Note that the hydrometeors have to be added to the df array before calling this 
    routine with dummy values. Information not avaialable here (e.g. scattering) are 
    still taken from the descriptor file. Note that all arrays of dataFullSpec mus be filled!
    Make sure that the array dimensions are not changed, otherwise Pamtra might crash.

    Returns
    -------
    dataFullSpec["rho_ds"] : 
      particle density at the middle of the bin in [kg/m3]
    dataFullSpec["d_ds"] : 
      particle diameter at the middle of the bin where the scattering properties shall be calculated
    dataFullSpec["d_bound_ds"] : 
      particle diameter at the edges of the bins. 5th dimension bin is one longer.
    dataFullSpec["n_ds"] : 
      particle number concentration in [#/m^3]. I.e. it is NOT normalized!
    dataFullSpec["mass_ds"] : 
      particle mass at the middle of the bin in [kg]
    dataFullSpec["area_ds"] : 
      particle cross section area at the middle of the bin in [m^2]
    dataFullSpec["as_ratio"] : 
      particle aspect ratio in the middle of the bin [no unit]
    dataFullSpec["canting"] : 
      particle canting the middle of the bin in [deg]

    Notes
    -----
    Do not forget to set nmlSet['hydro_fullspec'] to True.
    '''
    assert len(self.data) >0
    assert np.min(self.data["nbin"]) == np.max(self.data["nbin"])

    self.fs_nbin =  np.max(self.data["nbin"])
    self.parent._shape5Dplus = (self.parent.p["ngridx"],self.parent.p["ngridy"],self.parent.p["max_nlyrs"],self.nhydro,self.fs_nbin+1)
    self.parent._shape5D = (self.parent.p["ngridx"],self.parent.p["ngridy"],self.parent.p["max_nlyrs"],self.nhydro,self.fs_nbin)
    
    self.dataFullSpec["rho_ds"] = np.zeros(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["d_ds"] = np.zeros(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["d_bound_ds"] = np.zeros(self.parent._shape4D+(self.fs_nbin+1,))
    self.dataFullSpec["n_ds"] = np.zeros(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["mass_ds"] = np.zeros(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["area_ds"] = np.zeros(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["as_ratio"] = np.zeros(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["canting"] = np.zeros(self.parent._shape4D+(self.fs_nbin,))
  
