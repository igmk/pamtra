import numpy as np
import sys
from matplotlib import pyplot as plt
from matplotlib import colors, colorbar
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
import seaborn as sns
import pickle
import imp

import datetime
from netCDF4 import Dataset
from numpy.linalg import inv,det
from pyPamtra.importer import ncToDict
from geoLib import *

from struct import unpack

# define instrument specifications
# FOV: field of view
# pol: polarization at nadir view

instruments = {'amsua':
          {'name':'AMSU-A',
          'FOV':3.33,
          'footprint':[48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.,48.],
          'freq_str':[],
          'pol':['V','V','V','V','H','H','V','H','H','H','H','H','H','H']},
        'amsub':
          {'name':'AMSU-B',
          'FOV':1.1,
          'footprint':[16.,16.,16.,16.,16.],
          'freq_str':[],
          'pol':['V','V','H','H','V']},
        'mhs':
          {'name':'MHS',
          'FOV':1.1,
          'footprint':[16.,16.,16.,16.,16.],
          'freq_str':['89.0','157.0','183.311$\pm$1.0','183.311$\pm$3.0','190.311'],
          'pol':['V','V','H','H','V']}}
spacecraft = {'noaa18':
        {'name':'NOAA-18',
         'height':833.}}

class swath(object):
  def __init__(self):
    pass
  
  def read_swath(self,fname,sensor='amsub',satellite='noaa18'):
    p1 = 1.1910658e-5
    p2 = 1.438833

    self.sensor = sensor
    self.name = instruments[sensor]['name']
    self.fov = instruments[sensor]['FOV']
    self.footprint = instruments[sensor]['footprint']
    self.pol = instruments[sensor]['pol']
    self.height = spacecraft[satellite]['height']
    fid = open(fname, 'rb')
    # read header 
    head = Head()
    head.read_1bhead(fid)
    self.nrec = head.last_scan_record
    nrec = self.nrec
    self.sat_id = head.sat_id
    self.sensor = head.sensor
    sensor = self.sensor
    self.nchan = head.nchan
    nchan = self.nchan
    self.npos = head.npos
    npos = self.npos
    self.tcount = head.tcount
    tcount = self.tcount
    self.start_utc_of_day = head.start_utc_of_day
    self.end_utc_of_day = head.end_utc_of_day
    self.tr_central_wave_number = head.tr_central_wave_number
#        self.hour = (self.start_utc_of_day + 8000 * arange(nrec)) / (1000.*3600.)
    self.scanline_number = np.zeros([nrec],dtype=np.int64)
    self.scanline_year = np.zeros([nrec],dtype=np.int64)
    self.scanline_day_of_year = np.zeros([nrec],dtype=np.int64)
    self.scanline_utc_time_of_day = np.zeros([nrec],dtype=np.int64)
    self.date_time = np.empty([nrec], dtype=object)
    self.mode = np.zeros([nrec],dtype=np.int64)
    self.altitude = np.zeros([nrec])
    self.sza = np.zeros([nrec, npos])
    self.vza = np.zeros([nrec, npos])
    self.lat = np.zeros([nrec, npos])
    self.lon = np.zeros([nrec, npos])
    self.counta = np.zeros([nrec, npos, nchan],dtype=np.int64)
    self.countc = np.zeros([nrec, tcount, nchan],dtype=np.int64)
    self.countw = np.zeros([nrec, tcount, nchan],dtype=np.int64)
    if (sensor == 'A'):
        self.countwd = np.zeros([nrec,17],dtype=np.int64)
        self.countwa = np.zeros([nrec,3],dtype=np.int64)
    self.rada = np.zeros([nrec, npos, nchan])
    self.radc = np.zeros([nrec, tcount, nchan])
    self.radw = np.zeros([nrec, tcount, nchan])
    self.ta = np.zeros([nrec, npos, nchan])
    self.tc = np.zeros([nrec, tcount, nchan])
    self.tw = np.zeros([nrec, tcount, nchan])
    if (sensor == 'A'):
        self.twd = np.zeros([nrec,17])
        self.twa = np.zeros([nrec,3])
    self.clcoef = np.zeros([nrec, 3, nchan])

    if (sensor == 'A'):
        dcoef = np.zeros([4,17])
        acoef = np.zeros([2,3])
        dcoef[0,:] = head.warm_digital_coeffs[0,:] * 1.e-4
        dcoef[1,:] = head.warm_digital_coeffs[1,:] * 1.e-9
        dcoef[2,:] = head.warm_digital_coeffs[2,:] * 1.e-16
        dcoef[3,:] = head.warm_digital_coeffs[3,:] * 1.e-20
        acoef = head.warm_analog_coeffs * 1.e-3

    def convert_to_h_m_s_m(milliseconds):
        """Return the tuple of hours, minutes, seconds, and milliseconds."""
    
        seconds, milliseconds = divmod(milliseconds,1000)
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
#        days, hours = divmod(hours, 24)
    
        return hours, minutes, seconds, milliseconds
    # read data
    for irec in range(nrec):
        scanline = Scanline()
        if (sensor == 'A'):
            scanline.read_amsua_1bdata(fid,irec)
        elif (sensor == 'B') | ( sensor == 'S'):
            scanline.read_amsubmhs_1bdata(fid,irec)
        self.scanline_number[irec] = scanline.scanline_number
        self.scanline_year[irec] = scanline.scanline_year
        self.scanline_day_of_year[irec] = scanline.scanline_day_of_year
        self.scanline_utc_time_of_day[irec] = scanline.scanline_utc_time_of_day
        hours, minutes, seconds, milliseconds = convert_to_h_m_s_m(self.scanline_utc_time_of_day[irec])
        dt = datetime.date(self.scanline_year[0],1,1)+datetime.timedelta(self.scanline_day_of_year[0])
        self.date_time[irec] = datetime.datetime(self.scanline_year[irec],dt.month,dt.day,hours,minutes,seconds,milliseconds)
        self.mode[irec] = scanline.mode >> 15
        self.altitude[irec] = scanline.altitude *1.e-1
        self.sza[irec,:] = scanline.sza * 1.e-2
        self.vza[irec,:] = scanline.vza * 1.e-2
        self.lat[irec,:] = scanline.lat * 1.e-4
        self.lon[irec,:] = scanline.lon * 1.e-4
        self.counta[irec,:,:] = scanline.counta
        self.countc[irec,:,:] = scanline.countc
        self.countw[irec,:,:] = scanline.countw
        if (sensor == 'A'):
            self.countwd[irec,:] = scanline.countwd
            self.countwa[irec,:] = scanline.countwa
            self.clcoef[irec,0,:] = scanline.pri_cal_coeffs[0,:] * 1.e-19
            self.clcoef[irec,1,:] = scanline.pri_cal_coeffs[1,:] * 1.e-13
            self.clcoef[irec,2,:] = scanline.pri_cal_coeffs[2,:] * 1.e-9
        elif (sensor == 'B') | ( sensor == 'S'):
            self.clcoef[irec,0,:] = scanline.pri_cal_coeffs[0,:] * 1.e-16
            self.clcoef[irec,1,:] = scanline.pri_cal_coeffs[1,:] * 1.e-10
            self.clcoef[irec,2,:] = scanline.pri_cal_coeffs[2,:] * 1.e-6
        for ichan in range(nchan):
            self.rada[irec,:,ichan] = (self.clcoef[irec,2,ichan] + 
                self.clcoef[irec,1,ichan] * self.counta[irec,:,ichan] + 
                self.clcoef[irec,0,ichan] * self.counta[irec,:,ichan]** 2)
            self.radc[irec,:,ichan] = (self.clcoef[irec,2,ichan] + 
                self.clcoef[irec,1,ichan] * self.countc[irec,:,ichan] + 
                self.clcoef[irec,0,ichan] * self.countc[irec,:,ichan]** 2)
            self.radw[irec,:,ichan] = (self.clcoef[irec,2,ichan] + 
                self.clcoef[irec,1,ichan] * self.countw[irec,:,ichan] + 
                self.clcoef[irec,0,ichan] * self.countw[irec,:,ichan]** 2)
        if (sensor == 'A'):
                self.twd[irec,:] = (dcoef[0,:] + 
                    dcoef[1,:] * self.countwd[irec,:] + 
                    dcoef[2,:] * self.countwd[irec,:]**2 + 
                    dcoef[3,:] * self.countwd[irec,:]**3)
                self.twa[irec,:] = (acoef[0,:] + 
                    acoef[1,:] * self.countwa[irec,:])

    if (sensor == 'A'): self.twa = self.twa * 2.e-2
    
    fid.close()
    
    for ichan in range(nchan):
        pl = head.tr_central_wave_number[ichan] * 1.e-6
        c1v3 = p1 * pl ** 3
        c2v  = p2 * pl
        tmp  = np.log(c1v3 / self.rada[:,:,ichan] + 1)
        self.ta[:,:,ichan] = c2v / tmp
        tmp  = np.log(c1v3 / self.radc[:,:,ichan] + 1)
        self.tc[:,:,ichan] = c2v / tmp
        tmp  = np.log(c1v3 / self.radw[:,:,ichan] + 1)
        self.tw[:,:,ichan] = c2v / tmp

    return

    def __del__(self):
        pass

class Head:
    def __init__(self):
        pass
    
    def read_1bhead(self, fid):
        fid.seek(534,0)      # if archive header record (512) + 22
#        fid.seek(22, 0)                                         # 22
        self.data_set_name = unpack('>42c',fid.read(42)) # get data set name
        self.sensor = self.data_set_name[6]
        # set sensor specific data
        if (self.sensor == 'A'):
            self.nchan = 15
            self.npos = 30
            self.tcount = 2
            self.warm_digital_coeffs = np.zeros([4,17],dtype=np.int64)
            self.warm_analog_coeffs = np.zeros([2,3],dtype=np.int64)
        elif (self.sensor == 'B') | (self.sensor == 'S'):
            self.nchan = 5
            self.npos = 90
            self.tcount = 4

        self.tr_central_wave_number = np.zeros(self.nchan,dtype=np.int64)

        fid.seek(8,1)
        self.sat_id = unpack('>H',fid.read(2))[0]                     # 74
        self.ins_id = unpack('>2B', fid.read(2))[0]                   # 76
        
        fid.seek(8, 1)                                                # 84
        self.start_year = unpack('>H', fid.read(2))[0]                # 86
        self.start_day_of_year = unpack('>H', fid.read(2))[0]         # 88
        self.start_utc_of_day = unpack('>I', fid.read(4))[0]          # 92

        fid.seek(4, 1)                                                # 96
        self.end_year = unpack('>H', fid.read(2))[0]                  # 98
        self.end_day_of_year = unpack('>H', fid.read(2))[0]           # 100
        self.end_utc_of_day = unpack('>I', fid.read(4))[0]            # 104
        
        if (self.sensor == 'A'):
            fid.seek(40, 1)                                           # 144
            self.last_scan_record = unpack('>H', fid.read(2))[0]      # 146

            fid.seek(542, 1)
            for i in range(self.nchan):                                      # 688
                self.tr_central_wave_number[i] = unpack('>I', fid.read(4))[0]# 868
                fid.seek(8,1)
                
            fid.seek(668, 1)
            temp = np.array(unpack('>40i', fid.read(4*10*4))) # 1536
            self.warm_digital_coeffs[:,0:10] = temp.reshape([10,4]).transpose() # 1696

            fid.seek(36, 1)                                          # 1732
            temp = np.array(unpack('>4i', fid.read(2*2*4)))
            self.warm_analog_coeffs[:,0:2] = temp.reshape([2,2]).transpose()  # 1748

            fid.seek(364, 1)                                         # 2112
            self.warm_digital_coeffs[:,16] = unpack('>i', fid.read(4)) # 2128
            temp = np.array(unpack('>24i', fid.read(4*6*4)))
            self.warm_digital_coeffs[:,10:16] = temp.reshape([6,4]).transpose() # 2224

            fid.seek(28, 1)                                          # 2252
            self.warm_analog_coeffs[:,2] = unpack('>2i', fid.read(8))      # 2260

        elif (self.sensor == 'B'):
            fid.seek(28, 1)
            self.last_scan_record = unpack('>H', fid.read(2))[0]   # 134
        
            fid.seek(190, 1)                                         # 324
            for i in range(self.nchan):
                self.tr_central_wave_number[i] = unpack('>I', fid.read(4))[0]
                fid.seek(8,1)   # 384
                
        elif (self.sensor == 'S'):
            fid.seek(28,1)
            self.last_scan_record = unpack('>H', fid.read(2))[0]   # 134

            fid.seek(282, 1)                                         # 416
            for i in range(self.nchan):
                self.tr_central_wave_number[i] = unpack('>I', fid.read(4))[0]
                fid.seek(8,1)   # 476
class Scanline:
    def __init__(self):
        pass
    
    def read_amsua_1bdata(self,fid,irec):
        L_HBLOCK = 2560
        L_SCANLINE = 2560
        nchan = 15
        npos = 30
        
        self.counta = np.zeros([npos, nchan],dtype=np.int64)                       # antenna
        self.countc = np.zeros([2, nchan],dtype=np.int64)                          # cold target
        self.countw = np.zeros([2, nchan],dtype=np.int64)                          # warm target
        self.countwd = np.zeros(17,dtype=np.int64)                              # warm load digital
                    # 1-4, a1-1; 5, a1-1 center;
                    # 6-9, a1-2; 10, a1-2 center;
                    # 11-16, a2; 17, a2 center;
        self.countwa = np.zeros(3,dtype=np.int64);                           # warm load analog 
                    # 1, a1-1; 2, a1-2; 3, a2;
        
        fid.seek(L_HBLOCK + (irec) * L_SCANLINE + 512, 0)       # 0
#        fid.seek(12, 1)                                          # 12
        self.scanline_number = unpack('>H', fid.read(2))[0]
        self.scanline_year = unpack('>H', fid.read(2))[0]
        self.scanline_day_of_year = unpack('>H', fid.read(2))[0]
        fid.seek(2,1)
        self.scanline_utc_time_of_day = unpack('>I', fid.read(4))[0] # Scan Line UTC Time of Day in milliseconds
        self.mode = unpack('>H', fid.read(2))[0]                    # 14
        fid.seek(66, 1)                                          # 60
        self.pri_cal_coeffs = np.array(unpack('>45i', fid.read(3*nchan*4))).reshape([nchan,3]).transpose() # 120
        
        fid.seek(210, 1)                                          # 210
        self.altitude = unpack('>H', fid.read(2))[0]
        temp = np.array(unpack('>90h', fid.read(npos*3*2)))                         # 752
        self.sza = temp[::3]
        self.vza = temp[1::3]
        
        temp = np.array(unpack('>60i', fid.read(npos*2*4)))                         # 1472
        self.lat = temp[::2]
        self.lon = temp[1::2]

        fid.seek(12, 1)                                           # 1480
        temp = np.array(unpack('>510H', fid.read(1020)))                           # 2560
        for ichan in range(2,nchan):
            self.counta[:,ichan] = temp[ichan+2::17]
        
        temp = np.array(unpack('>30H', fid.read(60)))                            # 2616
        self.countc[0,2:15] = temp[4:17]
        self.countc[1,2:15] = temp[17::]
        
        temp = np.array(unpack('>46H', fid.read(92)))                            # 2616
        self.countwd[0:10] = temp[35:45]
        
        temp = np.array(unpack('>30H', fid.read(60)))                            # 2664
        self.countw[0,2:15] = temp[4:17]
        self.countw[1,2:15] = temp[17::]

        fid.seek(16, 1)                                          # 2152
        temp = unpack('>28B', fid.read(28))                             # 2180
        self.countwa[0:2] = temp[4:6]

        fid.seek(12, 1)                                          # 2192
        temp = unpack('>120H', fid.read(240))                           # 2432
        self.counta[:,0] = temp[0+2::4]
        self.counta[:,1] = temp[1+2::4]

        temp = unpack('>6H', fid.read(12))                             # 2444
        self.countc[0,0:2] = temp[2:4]
        self.countc[1,0:2] = temp[4::]

        temp = unpack('>20H', fid.read(40))                            # 2484
        self.countwd[16] = temp[12]
        self.countwd[10:16] = temp[13:19]

        temp = unpack('>6H', fid.read(12))                             # 2496
        self.countw[0,0:2] = temp[2:4]
        self.countw[1,0:2] = temp[4:6]

        fid.seek(19, 1)                                          # 2515
        self.countwa[2] = unpack('>B', fid.read(1))[0]
        
    def read_amsubmhs_1bdata(self,fid,irec):
        L_HBLOCK = 3072
        L_SCANLINE = 3072
        nchan = 5
        npos = 90
        
        self.counta = np.zeros([npos, nchan],dtype=np.int64)                       # antenna
        self.countc = np.zeros([4, nchan],dtype=np.int64)                          # cold target
        self.countw = np.zeros([4, nchan],dtype=np.int64)                          # warm target
        self.lat = np.zeros(npos,dtype=np.int64)
        self.lon = np.zeros(npos,dtype=np.int64)
        
        fid.seek(L_HBLOCK + (irec) * L_SCANLINE + 512, 0)       # 0
#        fid.seek(12, 1)                                          # 12
        self.scanline_number = unpack('>H', fid.read(2))[0]
        self.scanline_year = unpack('>H', fid.read(2))[0]
        self.scanline_day_of_year = unpack('>H', fid.read(2))[0]
        fid.seek(2,1)
        self.scanline_utc_time_of_day = unpack('>I', fid.read(4))[0]
        self.mode = unpack('>H', fid.read(2))[0]                    # 14
        fid.seek(46, 1)                                          # 60
        self.pri_cal_coeffs = np.array(unpack('>15i', fid.read(3*nchan*4))).reshape([nchan,3]).transpose() # 120
        
        fid.seek(90, 1)                                          # 210
        self.altitude = unpack('>H', fid.read(2))[0]
        temp = np.array(unpack('>270h', fid.read(npos*3*2)))                         # 752
        self.sza = temp[::3]
        self.vza = temp[1::3]
        
        temp = np.array(unpack('>180i', fid.read(npos*2*4)))                         # 1472
        self.lat = temp[::2]
        self.lon = temp[1::2]

        fid.seek(8, 1)                                           # 1480
        temp = np.array(unpack('>540H', fid.read(1080)))                           # 2560
        for ichan in range(nchan):
            self.counta[:,ichan] = temp[ichan+1::6]
        
        
        fid.seek(8, 1)                                           # 2568
        temp = np.array(unpack('>24H', fid.read(48)))                            # 2616
        self.countc[0,:] = temp[1:6]
        self.countc[1,:] = temp[7:12]
        self.countc[2,:] = temp[13:18]
        self.countc[3,:] = temp[19::]
        
        temp = np.array(unpack('>24H', fid.read(48)))                            # 2664
        self.countw[0,:] = temp[1:6]
        self.countw[1,:] = temp[7:12]
        self.countw[2,:] = temp[13:18]
        self.countw[3,:] = temp[19::]


    def __del__(self):
        pass    
  


Re = 6371 # Radius of th Earth

def getDistanceByHaversine(loc1, loc2):
   "Haversine formula - give coordinates as (lat_decimal,lon_decimal) tuples"

   lat1, lon1 = loc1
   lat2, lon2 = loc2

   # convert to radians
   lon1 = np.deg2rad(lon1)
   lon2 = np.deg2rad(lon2)
   lat1 = np.deg2rad(lat1)
   lat2 = np.deg2rad(lat2)

   # haversine formula
   dlon = lon2 - lon1
   dlat = lat2 - lat1
   a = (np.sin(dlat/2))**2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon/2.0))**2
   c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0-a))
   km = Re * c
   return km

def distdim(dist,fromto):

    if fromto == 'km2deg':
        return km2deg(dist)
    if fromto == 'deg2km':
        return deg2km(dist)
    else:
        sys.exit('No from to string')

def deg2km(deg):
    return np.deg2rad(Re * deg)

def km2deg(km):    
    return np.rad2deg(km/Re)
#    return km * 180 / (Re * pi)

  


def gauss(mu,covar,x):
  [n, d] = np.shape(x)

  [j, k] = np.shape(covar)

  # Check that the covariance matrix is the correct dimension
  if (j != d) | (k !=d):
      sys.exit('Dimension of the covariance matrix and data should match')
      
  
  invcov = inv(covar)

  mu = mu.reshape(1, d)    # Ensure that mu is a row vector

  x = x - np.ones([n, 1])*mu

  fact = np.sum((np.dot(x,invcov)*x), 1)

  y = np.exp(-0.5*fact)

  y = y/np.sqrt((2*np.pi)**d*det(covar))
  
  return y

class satConv(object):
    def __init__(self):
        self.o = dict()
        self.s = dict()
        self.c = dict()
    
        return
  
    def getSatData(self,fileName,sensor,obsChan):
        """
        Reads the satellite data and fills self.o
        """
        sw = swath.swath()

        sw.read_swath(fileName,sensor=sensor)

        self.o['fileName'] = fileName
        self.o['instrument'] = sensor
        self.o['fov'] = sw.fov
        self.o['footprint'] = sw.footprint[obsChan]
        self.o['pol'] = sw.pol[obsChan]
        self.o['height'] = sw.height

        if sensor == 'mhs':
            self.o['lon'] = sw.lon[1500::,:]
            self.o['lat'] = sw.lat[1500::,:]
            self.o['vza'] = sw.vza[1500::,:]
            self.o['tb'] = sw.ta[1500::,:,obsChan]
        else:
            self.o['lon'] = sw.lon
            self.o['lat'] = sw.lat
            self.o['vza'] = sw.vza
            self.o['tb'] = sw.ta[:,:,obsChan]

        del sw
        
        return
  
    def getSimData(self,fileName,pamChan,pamChan2=None):
        """
        Reads the simulated data for the requested channel and fills self.s. If a second channel is given, the mean of both is calculated.

        """
        
        pamData = ncToDict(fileName)
        self.s['fileName'] = fileName
        self.s['settings'] = dict()
        self.s['lon'] = pamData['longitude']
        self.s['lat'] = pamData['latitude']
        if 'lfrac' in pamData.keys():
            self.s['sfc_type'] = pamData['lfrac']
        else:
            self.s['sfc_type'] = pamData['sfc_type']
        self.s['frequency'] = pamData['frequency']
        self.s['tb'] = pamData['tb'][:,:,:,:,pamChan,:]
        if pamChan2:
            self.s['tb2'] = pamData['tb'][:,:,:,:,pamChan,:]
            self.s['tb'] = (self.s['tb']+self.s['tb2'])*0.5
        for wp in ['iwv','cwp','iwp','rwp','swp','gwp','hwp']:
          if wp in pamData.keys(): self.s[wp] = pamData[wp]
        self.s['ang'] = (pamData['angles']-180.)*(-1.)
         
        return
  
    def makeConvolution(self):
        """
        Does the convolution for one case and frequency.

        #
        # algorithm description
        #-----------------------
        #
        # * create a box around the available model area. area is defined by 
            border = [[self.s['lat'][0,0],self.s['lon'][0,0]],[self.s['lat'][0,-1],self.s['lon'][0,-1]],[self.s['lat'][-1,-1],self.s['lon'][-1,-1]],[self.s['lat'][-1,0],self.s['lon'][-1,0]]]
        # * find all foot prints of the satellite swath within
        #   the model area
        # * reduce satellite data to model area
        # * calculate the area (pi*r**2) that is relevant for a single foot print
        # * find nearest simulation point
        # * grab all points within circle of radius max(r). the number of points can 
        #   be found by translating max(r) in number of indices according to model resolution
        # * find the angle of observation that is necessary (do all the relevant simulations 
        #   have the same angle?)
        # * for each relevant simulation point calculate distance to foot print and 
        #   calculate its weight according to gaussian bell shape function
        # * calculate the mixture of polarization
        # * sum up the product between tb and the weights and normalize it
        # * store brightness temperatures and foot print coordinates

        """

        # We need to construct a projection, to see whether satellite observations are within the simulation
        # area.
        # Mercator Projection
        m = Basemap(projection='merc', llcrnrlat=-80, urcrnrlat=80,
            llcrnrlon=-180, urcrnrlon=180, lat_ts=20, resolution='c')
        # Poly vertices for simulation area [[lat[0,0],lon[0,0]],[lat[0,-1],lon[0,-1]],[lat[-1,-1],lon[-1,-1]],[lat[-1,0],lon[-1,0]]]

        border = [[self.s['lat'][0,0],self.s['lon'][0,0]],[self.s['lat'][0,-1],self.s['lon'][0,-1]],[self.s['lat'][-1,-1],self.s['lon'][-1,-1]],[self.s['lat'][-1,0],self.s['lon'][-1,0]]]

        # Projected vertices
        p_projected = [m(x[1], x[0]) for x in border]

        # Create the Path
        simBoundary = Path(p_projected)

        # construct an array of the indices in x and y direction to assign a footprint
        a = np.arange(self.o['vza'].shape[1]/2)
        a2 = np.concatenate((a[::-1],a))
        a3 = a2.reshape(1,self.o['vza'].shape[1])
        oind = np.repeat(a3,self.o['vza'].shape[0],axis=0) # array (len(swath),number_of_views)

        # reduce satellite data by finding all satellite observations within the simulation area.
        # the observations not needed are masked (mask=True)!
        d1, d2 = self.o['lon'].shape

        x,y = m(self.o['lon'].ravel(),self.o['lat'].ravel())
        points = np.c_[x,y]
        tmp = simBoundary.contains_points(points).reshape((d1,d2)) # bool array with true for observation points in simulation area

        # pdb.set_trace()
        self.c['indObs'] = np.where(tmp==True) # indices of observations within simulation area (array([rows]),array([cols])))

        # reduce the observations 2-d matrix to a vector containing only the observations in model area 
        self.c['obsTB'] = self.o['tb'][tmp]
        self.c['lon'] = self.o['lon'][tmp]
        self.c['lat'] = self.o['lat'][tmp]
        self.c['vza'] = self.o['vza'][tmp]
        self.c['ind'] = oind[tmp]
        
        self.c['frequency'] = self.s['frequency']
        
        # draw a circle around each observation footprint
        theta = np.deg2rad(self.c['vza'])
        zenith_mu = np.cos(theta)

        d3 = len(theta) # length of vector with relevant observations
        footprint = np.zeros(self.o['vza'].shape[1])
        self.c['simTB'] = np.zeros(d3)
        for wp in ['iwv','cwp','iwp','rwp','swp','gwp','hwp']:
            if wp in self.s.keys(): self.c[wp] = np.zeros(d3)

        # calculate distance around each foot print that needs to be included
        if self.o['fov']:
            footprint = self.o['height'] * (np.tan(np.deg2rad(self.o['vza'][0,:])+np.deg2rad(self.o['fov']/2.))-np.tan(np.deg2rad(self.o['vza'][0,:])))
        else:
            footprint[:] = self.o['footprint']/2.
            #  footprint = obsHeight * (tan(theta+np.deg2rad(self.o['fov']/2.))-tan(theta))
        # footprint[:] = 3.*16./2.354
        # for the big gridspacing of hirham we need to make the circle around each observation bigger
        # footprint[:] = 7.
        # convert from km to degrees:
        degdist = distdim(footprint,'km2deg')
        
        # simulation angles in radians reversed
        angle_mu = np.cos(np.deg2rad(self.s['ang'][0:16]))[::-1]
    

        # pdb.set_trace()
        # loop over observations in model domain
        for i in range(d3):
            needed = (getDistanceByHaversine([self.c['lat'][i],self.c['lon'][i]],[self.s['lat'],self.s['lon']]) < footprint[self.c['ind'][i]])
            # print i,footprint[self.c['ind'][i]],np.sum(needed)#, alldistssorted[0:5]
            if needed.any():
                need = dict()
                need['lat'] = self.s['lat'][needed]
                need['lon'] = self.s['lon'][needed]

                need['dists'] = np.zeros(np.size(need['lat']))
                need['tb'] = np.zeros(np.size(need['lat']))
                for wp in ['iwv','cwp','iwp','rwp','swp','gwp','hwp']:
                    if wp in self.s.keys(): need[wp] = np.zeros(np.size(need['lat']))
                index_needed = needed.nonzero()
                # pdb.set_trace()
                if self.o['instrument'] in ['amsua','amsub','mhs']:
                    if self.o['pol'] == 'V':
                        for k in range(np.size(need['lat'])):
                            need['dists'][k] = getDistanceByHaversine([need['lat'][k],need['lon'][k]],[self.c['lat'][i],self.c['lon'][i]])
                            need['tb'][k] = (np.interp(zenith_mu[i],angle_mu,self.s['tb'][index_needed[0][k],index_needed[1][k],0,0:16,0][::-1])*0.5*(1+np.cos(2*theta[i]))
                                  +np.interp(zenith_mu[i],angle_mu,self.s['tb'][index_needed[0][k],index_needed[1][k],0,0:16,1][::-1])*0.5*(1-np.cos(2*theta[i])))
                    elif self.o['pol'] == 'H':
                        for k in range(np.size(need['lat'])):
                            need['dists'][k] = getDistanceByHaversine([need['lat'][k],need['lon'][k]],[self.c['lat'][i],self.c['lon'][i]])
                            need['tb'][k] = (np.interp(zenith_mu[i],angle_mu,self.s['tb'][index_needed[0][k],index_needed[1][k],0,0:16,1][::-1])*0.5*(1+np.cos(2*theta[i]))
                                  +np.interp(zenith_mu[i],angle_mu,self.s['tb'][index_needed[0][k],index_needed[1][k],0,0:16,0][::-1])*0.5*(1-np.cos(2*theta[i])))
                elif self.o['instrument'] in ['ssmi','ssmis']:
                    if self.o['pol'] == 'V':
                        for k in range(np.size(need['lat'])):
                            need['dists'][k] = getDistanceByHaversine([need['lat'][k],need['lon'][k]],[self.c['lat'][i],self.c['lon'][i]])
                            need['tb'][k] = np.interp(zenith_mu[i],angle_mu,self.s['tb'][index_needed[0][k],index_needed[1][k],0,0:16,0][::-1])
                    elif self.o['pol'] == 'H':
                        for k in range(np.size(need['lat'])):
                            need['dists'][k] = getDistanceByHaversine([need['lat'][k],need['lon'][k]],[self.c['lat'][i],self.c['lon'][i]])
                            need['tb'][k] = np.interp(zenith_mu[i],angle_mu,self.s['tb'][index_needed[0][k],index_needed[1][k],0,0:16,1][::-1])
                    elif self.o['pol'] == 'RC':
                        for k in range(np.size(need['lat'])):
                            need['dists'][k] = getDistanceByHaversine([need['lat'][k],need['lon'][k]],[self.c['lat'][i],self.c['lon'][i]])
                            need['tb'][k] = 0.5*(np.interp(zenith_mu[i],angle_mu,self.s['tb'][index_needed[0][k],index_needed[1][k],0,0:16,0][::-1])+np.interp(zenith_mu[i],angle_mu,self.s['tb'][index_needed[0][k],index_needed[1][k],0,0:16,1][::-1]))
                for wp in ['iwv','cwp','iwp','rwp','swp','gwp','hwp']:
                    if wp in self.s.keys(): 
                        for k in range(np.size(need['lat'])):
                            need[wp][k] = self.s[wp][index_needed[0][k],index_needed[1][k]]
                sigma = footprint[self.c['ind'][i]]/3.
                # sigma = 16./2.354
                bell = gauss(np.array([0]),np.array([sigma]).reshape(1,1),need['dists'].reshape(len(need['dists']),1))
                alph = 1/np.sum(bell)
                weights = alph * bell
                self.c['simTB'][i] = np.sum(weights*need['tb'].transpose())
                for wp in ['iwv','cwp','iwp','rwp','swp','gwp','hwp']:
                    if wp in self.s.keys(): self.c[wp][i] = np.sum(weights*need[wp].transpose())
      
            else:
                self.c['simTB'][i] = float('nan')
      
        return



    def _reduce(self,dat):

        d1,d2 = self.o['lon'].shape
        dat_tmp = np.zeros((d1,d2))
        dat_tmp[self.c['indObs']] = dat
        dat_ma = np.ma.masked_where(dat_tmp==0., dat_tmp)

        return dat_ma



sns.set(style="whitegrid", palette="pastel", color_codes=True)

def performConvolution(ssrgFile,mieFile,satFile,sensor,simChan,satChan):

    sc = satConv.satConv()

    sc.getSimData(ssrgFile,simChan)

    sc.getSatData(satFile,sensor,satChan)

    sc.makeConvolution()

    scm = satConv.satConv()

    scm.getSimData(mieFile,simChan)

    scm.getSatData(satFile,sensor,satChan)

    scm.makeConvolution()

    return sc,scm

def makePlot(fig, ax_pos_grid0, sc,scm, hyd, lat, lon, row, norm, normd, nhyd, cmap='jet'):

    freqs = ['50.3', '89.0', '157.0']
    lhyds = ['CWP+RWP','IWP', 'SWP']
    olabels = ['(a)','(b)','(c)']
    mlabels = ['(d)','(e)','(f)']
    slabels = ['(g)','(h)','(i)']
    domlabels = ['(j)','(k)','(l)']
    doslabels = ['(m)','(n)','(o)']
    dsmlabels = ['(p)','(q)','(r)']
    hlabels = ['(s)','(t)','(u)']
    ancoords = (0.02,0.85)

    ax = fig.add_subplot(ax_pos_grid0[row,0])

    map = Basemap(projection='merc',llcrnrlat=50.,urcrnrlat=60.,llcrnrlon=-42.5,urcrnrlon=-27.5,lat_ts=55.,resolution='c')

    X,Y = map(sc._reduce(sc.c['lon']),sc._reduce(sc.c['lat']))

    map.contourf(X,Y,sc._reduce(sc.c['obsTB']),25,norm=norm,cmap=cmap)
    if row == 0: plt.title('Observation')
    ax.text(-0.2, 0.5, freqs[row]+' GHz', fontsize=14, va='center', rotation=90,transform=ax.transAxes)
    ax.annotate(olabels[row],ancoords, xycoords='axes fraction',fontweight='bold')

    axm = fig.add_subplot(ax_pos_grid0[row,1])

    map = Basemap(projection='merc',llcrnrlat=50.,urcrnrlat=60.,llcrnrlon=-42.5,urcrnrlon=-27.5,lat_ts=55.,resolution='c')

    X,Y = map(sc._reduce(scm.c['lon']),sc._reduce(scm.c['lat']))

    map.contourf(X,Y,scm._reduce(scm.c['simTB']),25,norm=norm,cmap=cmap)
    if row == 0: plt.title('Mie')
    axm.annotate(mlabels[row],ancoords, xycoords='axes fraction',fontweight='bold')

    axs = fig.add_subplot(ax_pos_grid0[row,2])

    map = Basemap(projection='merc',llcrnrlat=50.,urcrnrlat=60.,llcrnrlon=-42.5,urcrnrlon=-27.5,lat_ts=55.,resolution='c')

    X,Y = map(sc._reduce(sc.c['lon']),sc._reduce(sc.c['lat']))

    map.contourf(X,Y,sc._reduce(sc.c['simTB']),25,norm=norm,cmap=cmap)
    if row == 0: plt.title('SSRGA')
    axs.annotate(slabels[row],ancoords, xycoords='axes fraction',fontweight='bold')

    if row == 2:
        axc_box = ax.get_position()
        axc1_box = axs.get_position()
        cax = fig.add_axes([axc_box.x0, axc_box.y0+0.12, axc1_box.x1-axc_box.x0, 0.02], rasterized=True)
        cb1 = colorbar.ColorbarBase(cax,norm=norm,orientation='horizontal',cmap='jet')

        cb1.set_label('$\mathrm{T_B}$  [K]')

    # the differences

    axdom = fig.add_subplot(ax_pos_grid0[row,3])

    map = Basemap(projection='merc',llcrnrlat=50.,urcrnrlat=60.,llcrnrlon=-42.5,urcrnrlon=-27.5,lat_ts=55.,resolution='c')

    X,Y = map(sc._reduce(sc.c['lon']),sc._reduce(sc.c['lat']))

    map.contourf(X,Y,sc._reduce(sc.c['obsTB']-scm.c['simTB']),25,norm=normd,cmap='bwr')

    if row == 0: plt.title('Obs.-Mie')
    axdom.annotate(domlabels[row],ancoords, xycoords='axes fraction',fontweight='bold')

    axdos = fig.add_subplot(ax_pos_grid0[row,4])

    map = Basemap(projection='merc',llcrnrlat=50.,urcrnrlat=60.,llcrnrlon=-42.5,urcrnrlon=-27.5,lat_ts=55.,resolution='c')

    X,Y = map(sc._reduce(sc.c['lon']),sc._reduce(sc.c['lat']))

    map.contourf(X,Y,sc._reduce(sc.c['obsTB']-sc.c['simTB']),25,norm=normd,cmap='bwr')

    if row == 0: plt.title('Obs.-SSRGA')
    axdos.annotate(doslabels[row],ancoords, xycoords='axes fraction',fontweight='bold')

    axdsm = fig.add_subplot(ax_pos_grid0[row,5])

    map = Basemap(projection='merc',llcrnrlat=50.,urcrnrlat=60.,llcrnrlon=-42.5,urcrnrlon=-27.5,lat_ts=55.,resolution='c')

    X,Y = map(sc._reduce(sc.c['lon']),sc._reduce(sc.c['lat']))

    map.contourf(X,Y,sc._reduce(sc.c['simTB']-scm.c['simTB']),25,norm=normd,cmap='bwr')

    if row == 0: plt.title('SSRGA-Mie')
    axdsm.annotate(dsmlabels[row],ancoords, xycoords='axes fraction',fontweight='bold')

    if row == 2:
        axd1_box = axdom.get_position()
        axd2_box = axdsm.get_position()
        caxd = fig.add_axes([axd1_box.x0, axd1_box.y0+0.12, axd2_box.x1 - axd1_box.x0, 0.02], rasterized=True)
        cb = colorbar.ColorbarBase(caxd,norm=normd,orientation='horizontal',cmap='bwr')
        cb.set_ticks([-80,-40, 0 ,40, 80])

        cb.set_label('$\Delta \mathrm{T_B}$ [K]')

    axhyd = fig.add_subplot(ax_pos_grid0[row,6])

    map = Basemap(projection='merc',llcrnrlat=50.,urcrnrlat=60.,llcrnrlon=-42.5,urcrnrlon=-27.5,lat_ts=55.,resolution='c')

    X,Y = map(lon,lat)

    map.contourf(X,Y,hyd,25,norm=nhyd,cmap='binary')

    axhyd.set_title(lhyds[row])
    axhyd.annotate(hlabels[row],ancoords, xycoords='axes fraction',fontweight='bold')

    if row == 2:
        axd1_box = axhyd.get_position()
        caxd = fig.add_axes([axd1_box.x0, axd1_box.y0+0.12, axd1_box.x1 - axd1_box.x0, 0.02], rasterized=True)
        cb = colorbar.ColorbarBase(caxd,norm=nhyd,orientation='horizontal',cmap='binary')

        cb.set_label('WP $[\mathrm{kg\;m^{-2}}]$')

    return

makeConv = False

fig = plt.figure(figsize=(12,9))

gs1 = gridspec.GridSpec(3,7)
gs1.update(bottom=0.05, top=0.95,left=0.1, right=0.90, hspace=-0.6, wspace=0.05)

satFile = ['NSS.AMAX.NN.D16270.S1743.E1938.B5850708.GC','NSS.MHSX.NN.D16270.S1743.E1938.B5850708.GC','NSS.MHSX.NN.D16270.S1743.E1938.B5850708.GC']
sensor = ['amsua','mhs','mhs']
fstrs = ['50','89','157']
simChan = [0,1,2]
satChan = [2,0,1]

if makeConv:
    # Create the convolution
    sc = dict()
    scm = dict()
    for i,fstr in enumerate(fstrs):
        sc[i], scm[i] = performConvolution('../data/ecmwf_ssrga.nc','../data/ecmwf_mie.nc','../data/'+satFile[i],sensor[i],simChan[i],satChan[i])

    with open('../data/conv_ssrg_v2.pkl', 'wb') as handle:
        pickle.dump(sc, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open('../data/conv_mie_v2.pkl', 'wb') as handle:
        pickle.dump(scm, handle, protocol=pickle.HIGHEST_PROTOCOL)

else:
    # Work on existing convolution
    with open('../data/conv_ssrg_v2.pkl', 'rb') as handle:
        sc = pickle.load(handle)

    with open('../data/conv_mie_v2.pkl', 'rb') as handle:
        scm = pickle.load(handle)

    with open('../data/hydros.pkl', 'rb') as handle:
        hyds = pickle.load(handle)
        lat = pickle.load(handle)
        lon = pickle.load(handle)

smin = np.zeros(len(fstrs))
smax = np.zeros(len(fstrs))
mmin = np.zeros(len(fstrs))
mmax = np.zeros(len(fstrs))
omin = np.zeros(len(fstrs))
omax = np.zeros(len(fstrs))
cmin = np.zeros(len(fstrs))
cmax = np.zeros(len(fstrs))
for i,fstr in enumerate(fstrs):
    smin[i] = np.nanmin(sc[i].c['simTB'])
    smax[i] = np.nanmax(sc[i].c['simTB'])
    mmin[i] = np.nanmin(scm[i].c['simTB'])
    mmax[i] = np.nanmax(scm[i].c['simTB'])
    omin[i] = np.nanmin(scm[i].c['obsTB'])
    omax[i] = np.nanmax(scm[i].c['obsTB'])
    cmin[i] = np.floor(min(smin[i],mmin[i],omin[i]))
    cmax[i] = np.ceil(max(smax[i],mmax[i],omax[i]))
hyd_max = np.nanmax([np.nanmax(hyds[0][:].ravel()),np.nanmax(hyds[1][:].ravel()),np.nanmax(hyds[2][:].ravel())])
norm = colors.Normalize(vmin=np.nanmin(cmin), vmax=np.nanmax(cmax))
normd = colors.Normalize(vmin=np.nanmin(smin-mmin), vmax=-1.*np.nanmin(smin-mmin))
nhyd = colors.Normalize(vmin=0.,vmax=hyd_max)

for i,fstr in enumerate(fstrs):
    makePlot(fig, gs1, sc[i], scm[i], hyds[i], lat, lon, i,norm,normd,nhyd)
