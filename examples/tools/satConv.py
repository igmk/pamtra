import swath
import datetime
from netCDF4 import Dataset
import numpy as np
from numpy.linalg import inv,det
from pyPamtra.importer import ncToDict
from matplotlib import pyplot, colors, colorbar
from mpl_toolkits.basemap import Basemap
from matplotlib.path import Path
from geoLib import *
import pdb
import imp

imp.reload(swath)


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

