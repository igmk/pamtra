# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import csv

class pamDescriptorFile(object):
  #class with the descriptor file content in data. In case you want to use 4D data, use data4D instead, the coreesponding column in data is automatically removed.
  
         
  def __init__(self,parent):
    self.names = np.array(["hydro_name", "as_ratio", "liq_ice", "rho_ms", "a_ms", "b_ms", "alpha_as", "beta_as", "moment_in", "nbin", "dist_name", "p_1", "p_2", "p_3", "p_4", "d_1", "d_2", "scat_name", "vel_size_mod","canting"])
    self.types = ["S15",float,int,float,float,float,float,float,int,int,"S15",float,float,float,float,float,float, "S15", "S30",float]  
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
      assert len(row) == len(self.names)  

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
  
    #for name in self.names:
      #if name not in ["hydro_name","liq_ice","scat_name","vel_size_mod"]:
        #self.data = mlab.rec_drop_fields(self.data,[name])
    #return
