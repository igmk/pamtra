# -*- coding: utf-8 -*-

from __future__ import division, print_function
import numpy as np
import csv
from scipy.interpolate import interp1d

missingNumber=-9999.
missingIntNumber=int(missingNumber)

class pamDescriptorFile(object):
  '''
  class with the descriptor file content in data. In case you want to use 4D data, use data4D instead, 
  the coreesponding column in data is automatically removed. Inializes with an empty array.
  '''
         
  def __init__(self,parent):
    self.names = np.array(["hydro_name", "as_ratio", "liq_ice", "rho_ms", "a_ms", "b_ms", "alpha_as", "beta_as", "moment_in", "nbin", "dist_name", "p_1", "p_2", "p_3", "p_4", "d_1", "d_2", "scat_name", "vel_size_mod","canting"])
    self.types = ["U15",float,int,float,float,float,float,float,int,int,"U15",float,float,float,float,float,float, "U60", "U30",float]  
    self.data = np.recarray((0,),dtype=list(zip(self.names, self.types)))
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
        if self.types[ii] in ("S15", "U15"):
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
    Add hydrometeor to df array. See :any:`descriptorFile` for details.

    Parameters
    ----------
    hydroTuple : tuple
      tuple describing the hydrometeor. Consists of
      "hydro_name", "as_ratio", "liq_ice", "rho_ms", "a_ms", "b_ms", "alpha_as", "beta_as", "moment_in", "nbin", "dist_name", "p_1", 
      "p_2", "p_3", "p_4", "d_1", "d_2", "scat_name", "vel_size_mod","canting"
    '''
    assert hydroTuple[0] not in self.data["hydro_name"]
    assert len(list(self.dataFullSpec.keys())) == 0
    
    self.data = np.append(self.data,np.array(tuple(hydroTuple),dtype=list(zip(self.names,self.types))))
    self.nhydro += 1
    self.parent._shape4D = (self.parent.p["ngridx"],self.parent.p["ngridy"],self.parent.p["max_nlyrs"],self.nhydro)
    for key in ["hydro_q","hydro_reff","hydro_n"]:
      if key in list(self.parent.p.keys()):    
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

    assert len(list(self.dataFullSpec.keys())) == 0
    if hydroName == "all":
      self.__init__()
      return
    removed = False#
    hydroIndex = np.where(self.parent.df.data["hydro_name"]==hydroName)[0][0]
    self.dataNew = np.recarray((0,),dtype=list(zip(self.names, self.types)))
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
        if key in list(self.parent.p.keys()):
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
    still taken from the descriptor file. Note that all arrays of dataFullSpec mus be filled and
    you still need to define the hydrometeors if the df object! However, only name, liq_ice,
    nbins+1, scat_name, and vel_size_mod is used even though you have to define the other
    variables as well.
    
    Make sure that the array dimensions are not changed, otherwise Pamtra might crash.

    The following objects are added to the df object:
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
    dataFullSpec["fallvelocity"] : 
      particle fall velocity in [m/s]. Only used if any value > 0 m/s, 
      otherwise the fall velocity relation provided in the descriptor file is 
      used.

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
    self.dataFullSpec["fallvelocity"] = np.zeros(self.parent._shape4D+(self.fs_nbin,))
    self.dataFullSpec["rg_kappa_ds"] = np.full(self.parent._shape4D+(self.fs_nbin,), 0.19)
    self.dataFullSpec["rg_beta_ds"] = np.full(self.parent._shape4D+(self.fs_nbin,), 0.23)
    self.dataFullSpec["rg_gamma_ds"] = np.full(self.parent._shape4D+(self.fs_nbin,), 5./3.)
    self.dataFullSpec["rg_zeta_ds"] = np.full(self.parent._shape4D+(self.fs_nbin,), 1.)


"""
Helper functions to get SSRGA parameter from normalized rime mass M follwing Maherndl et al.[2022, in prep.]
"""



def ssrga_parameter(M):
    '''
    Calculate ssrga parameter kappa, beta, gamma, zeta1 and alpha_eff from the normalized rime mass M as defined in Seifert et al. (2019), 
    M = mass of rime / mass of size equivalent graupel assuming graupel density of 700kg/m³

    Input:
    M (unitless)

    Output
    kappa, beta, gamma, zeta1, alpha_eff
    '''

    p0 = 0.514
    p1 = np.array([0.16, -0.1, 4.06, -1.27, 0.127])
    p2 = np.array([0.187, 0.068, -7.45, 1.79, -0.091])
    p3 = np.array([0.575, 0.194, 5.42, 2.76, 0.067])

    alpha_eff = p1[0]*M**(2*p0)+p2[0]*M**p0+p3[0]
    kappa =   p1[1]*M**(2*p0)+p2[1]*M**p0+p3[1]
    beta =    p1[2]*M**(2*p0)+p2[2]*M**p0+p3[2]
    gamma =   p1[3]*M**(2*p0)+p2[3]*M**p0+p3[3]
    zeta1 =   p1[4]*M**(2*p0)+p2[4]*M**p0+p3[4]

    return kappa, beta, gamma, zeta1, alpha_eff

def riming_dependent_ssrga_name(M):
    '''
    Return string scattering name for ss-rayleigh-gans in form "ss-rayleigh-gans_kappa_beta_gamma_zeta"
    Calculate ssrga parameter kappa, beta, gamma, zeta1 and alpha_eff from the normalized rime mass M as defined in Seifert et al. (2019), 
    M = mass of rime / mass of size equivalent graupel assuming graupel density of 700kg/m³

    Input:
    M (unitless)

    Output
    scatering name string
    '''
    try:
        para = ssrga_parameter(M)
        stringname = 'ss-rayleigh-gans_' + str(np.round(ssrga_parameter(M), 3)[0]) + '_' + str(np.round(ssrga_parameter(M), 3)[1]) + '_' + str(np.round(ssrga_parameter(M), 3)[2]) + '_' + str(np.round(ssrga_parameter(M), 3)[3]) 
        
            
        if len(stringname) > 40:
            print('enter single M value as float, not array!')
    except:
        print('enter single M value as float!')
        
    return stringname

def riming_dependent_mass_size(M, monomer):
    '''
    Return mass size parameter a_m and b_m for given monomer type from the normalized rime mass M as defined in Seifert et al. (2019), 
    M = mass of rime / mass of size equivalent graupel assuming graupel density of 700kg/m³
    a_m and b_m are interpolated from Table 2 in Maherndl et al [in prep.], for M>0.8155 values for M=0.8155 are used
    Input:
    M (unitless)
    monomer type (column, dendrite, needle, plate, rosette, mean)

    Output
    a_m, b_m
    '''
    M_list = np.array([0.0, 0.0129,0.02045,0.03245,0.05145,0.08155,0.129,0.2045,0.3245,0.5145,0.8155])
    if hasattr(M, "__len__"): 
        M[M>M_list[-1]] = M_list[-1]

    else:
         if M > M_list[-1]:
            M = M_list[-1] 


    if monomer == 'column':
        a_m_list = np.array([0.0485, 0.0988, 0.210, 0.638, 1.91, 4.74, 12.5, 28.6, 56.2, 128., 166.])
        b_m_list = np.array([2.07, 2.17, 2.26, 2.40, 2.54, 2.64, 2.74, 2.82, 2.87, 2.93, 2.92])

        a_int = interp1d(M_list, a_m_list, kind='cubic')
        b_int = interp1d(M_list, b_m_list, kind='cubic')


        a_m = a_int(M)
        b_m = b_int(M)

    elif monomer == 'dendrite':
        a_m_list = np.array([0.0132, 0.388, 1.00, 2.77, 7.26, 16.4, 32.9, 59.4, 98.6, 173., 143.])
        b_m_list = np.array([2.09, 2.52, 2.63, 2.73, 2.82, 2.89, 2.93, 2.96, 2.97, 2.99, 2.90])

        a_int = interp1d(M_list, a_m_list, kind='cubic')
        b_int = interp1d(M_list, b_m_list, kind='cubic')

        a_m = a_int(M)
        b_m = b_int(M)

    elif monomer == 'needle':
        a_m_list = np.array([0.0254, 0.136, 0.336, 1.05, 3.01, 7.77, 19.1, 39.4, 78.0, 143., 184.])
        b_m_list = np.array([2.06, 2.28, 2.39, 2.53, 2.64, 2.74, 2.83, 2.88, 2.93, 2.96, 2.94])

        a_int = interp1d(M_list, a_m_list, kind='cubic')
        b_int = interp1d(M_list, b_m_list, kind='cubic')

        a_m = a_int(M)
        b_m = b_int(M)

    elif monomer == 'plate':
        a_m_list = np.array([0.0388, 0.219, 0.508, 1.44, 4.20, 9.97, 21.8, 42.8, 81.5, 160., 209.])
        b_m_list = np.array([2.14, 2.37, 2.47, 2.59, 2.71, 2.79, 2.85, 2.90, 2.94, 2.98, 2.95])

        a_int = interp1d(M_list, a_m_list, kind='cubic')
        b_int = interp1d(M_list, b_m_list, kind='cubic')

        a_m = a_int(M)
        b_m = b_int(M)

    elif monomer == 'rosette':
        a_m_list = np.array([0.0363, 0.277, 0.629, 1.80, 4.84, 11.6, 24.9, 46.4, 81.0, 182., 165.])
        b_m_list = np.array([2.13, 2.40, 2.50, 2.62, 2.73, 2.81, 2.87, 2.91, 2.94, 2.99, 2.92])

        a_int = interp1d(M_list, a_m_list, kind='cubic')
        b_int = interp1d(M_list, b_m_list, kind='cubic')

        a_m = a_int(M)
        b_m = b_int(M)

    elif monomer == 'mean':
        a_m_list = np.array([0.0324, 0.224, 0.537, 1.54, 4.27, 10.1, 22.2, 43.3, 79.0, 157., 173.])
        b_m_list = np.array([2.10, 2.35, 2.45, 2.57, 2.69, 2.77, 2.85, 2.89, 2.93, 2.97, 2.93])

        a_int = interp1d(M_list, a_m_list, kind='cubic')
        b_int = interp1d(M_list, b_m_list, kind='cubic')

        a_m = a_int(M)
        b_m = b_int(M)
    else:
        print('enter valid monomer type!')
        
    return a_m, b_m

def riming_dependent_area_size(M, monomer):
    '''
    Return area size parameter a_A and b_A for given monomer type from the normalized rime mass M as defined in Seifert et al. (2019), 
    M = mass of rime / mass of size equivalent graupel assuming graupel density of 700kg/m³
    a_A and b_A are interpolated from Table 3 in Maherndl et al [in prep.], for M>0.8155 values for M=0.8155 are used
    Input:
    M (unitless)
    monomer type (column, dendrite, needle, plate, rosette, mean)

    Output
    a_A, b_A
    '''
    M_list = np.array([0.0, 0.0129,0.02045,0.03245,0.05145,0.08155,0.129,0.2045,0.3245,0.5145,0.8155])
    if hasattr(M, "__len__"): 
        M[M>M_list[-1]] = M_list[-1]

    else:
         if M > M_list[-1]:
            M = M_list[-1]       

    if monomer == 'column':
        a_A_list = np.array([0.0508, 0.0392, 0.0506, 0.0831, 0.127, 0.161, 0.227, 0.298, 0.359, 0.427, 0.485])
        b_A_list = np.array([1.77, 1,73, 1.77, 1.83, 1.88, 1.90, 1.94, 1.97, 1.98, 1.99, 1.99])

        a_int = interp1d(M_list, a_A_list, kind='cubic')
        b_int = interp1d(M_list, b_A_list, kind='cubic')


        a_A = a_int(M)
        b_A = b_int(M)

    elif monomer == 'dendrite':
        a_A_list = np.array([0.0762, 0.128, 0.144, 0.182, 0.206, 0.215, 0.250, 0.298, 0.341, 0.417, 0.453])
        b_A_list = np.array([1.87, 1.93, 1.93, 1.95, 1.95, 1.94, 1.95, 1.96, 1.97, 1.98, 1.98])

        a_int = interp1d(M_list, a_A_list, kind='cubic')
        b_int = interp1d(M_list, b_A_list, kind='cubic')

        a_A = a_int(M)
        b_A = b_int(M)

    elif monomer == 'needle':
        a_A_list = np.array([0.0428, 0.0577, 0.0681, 0.110, 0.154, 0.190, 0.250, 0.316, 0.375, 0.434, 0.496])
        b_A_list = np.array([1.78, 1.79, 1.83, 1.89, 1.92, 1.94, 1.96, 1.98, 1.99, 1.99, 1.99])

        a_int = interp1d(M_list, a_A_list, kind='cubic')
        b_int = interp1d(M_list, b_A_list, kind='cubic')

        a_A = a_int(M)
        b_A = b_int(M)

    elif monomer == 'plate':
        a_A_list = np.array([0.0715, 0.0597, 0.0689, 0.0935, 0.134, 0.159, 0.210, 0.278, 0.368, 0.423, 0.512])
        b_A_list = np.array([1.81, 1.78, 1.80, 1.83, 1.88, 1.89, 1.92, 1.95, 1.98, 1.98, 2.00])

        a_int = interp1d(M_list, a_A_list, kind='cubic')
        b_int = interp1d(M_list, b_A_list, kind='cubic')

        a_A = a_int(M)
        b_A = b_int(M)

    elif monomer == 'rosette':
        a_A_list = np.array([0.0643, 0.0748, 0.0902, 0.129, 0.168, 0.192, 0.235, 0.276, 0.338, 0.402, 0.455])
        b_A_list = np.array([1.81, 1.83, 1.85, 1.89, 1.91, 1.92, 1.94, 1.95, 1.96, 1.97, 1.98])

        a_int = interp1d(M_list, a_A_list, kind='cubic')
        b_int = interp1d(M_list, b_A_list, kind='cubic')

        a_A = a_int(M)
        b_A = b_int(M)

    elif monomer == 'mean':
        a_A_list = np.array([0.0611, 0.0699, 0.0843, 0.120, 0.158, 0.183, 0.234, 0.293, 0.356, 0.421, 0.480])
        b_A_list = np.array([1.81, 1.81, 1.83, 1.88, 1.91, 1.92, 1.94, 1.96, 1.97, 1.98, 1.99])

        a_int = interp1d(M_list, a_A_list, kind='cubic')
        b_int = interp1d(M_list, b_A_list, kind='cubic')

        a_A = a_int(M)
        b_A = b_int(M)
    else:
        print('enter valid monomer type!')
        
    return a_A, b_A    
    