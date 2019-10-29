# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager 
from copy import deepcopy
import warnings

try:
  import netCDF4 as nc
  pyNc = True
except:
  warnings.warn("not tested for other netcdf implementations!")
  try:
    import Scientific.IO.NetCDF as nc
    pyNc = False 
  except:
    import netCDF3 as nc
    pyNc = True

missingNumber =-9999.
def niceColors(length,cmap='hsv'):
  colors = list()
  cm = plt.get_cmap(cmap)
  for l in range(length):
    colors.append(cm(1.*l/length))
  return colors
  
def plotTB(data,polarisation="H",angle=180.,outlevel="space",frequency=0,xgrid="index",ygrid="index",levels=["min","max"],cmap=None,title="",fig =None):
  """
  plot brightness temperatures from pamtra
  
  input:
    data: netcdf file OR pyPamtra object
    polarisatuion: H or V
    outlevel: space or ground
    frequency: if real number frequency, if int, index of frequency
    xgrid: index or?
    ygrid index or? 
  """
  #first parse options for both netcdf and pyPamtra
  if polarisation=="H":
    stokesIndex =1
  elif polarisation=="V":
    stokesIndex =0
  else:
    raise ValueError("polarisation must be H or V")

  if outlevel=="space":
    levelIndex =0
  elif outlevel=="ground":
    levelIndex =1
  else:
    raise ValueError("outlevel must be space or ground")
    
  #now specificly netcdf OR pyPamtra
  if type(data) == str:
    #open netcdf etc...
    if pyNc: cdfFile = nc.Dataset(data,"r")
    else: cdfFile = nc.NetCDFFile(data,"r")
    
    try:
      angleIndex = np.where(angle==cdfFile.variables["angle"][:])[0][0]
    except IndexError: 
      raise ValueError("angle was not found")
    if type(frequency) in [int,np.int32,np.int64]:
      if frequency < len(cdfFile.variables["frequency"][:]):
        frequencyIndex=frequency
      else:
        raise ValueError("frequencyIndex was not found. Use integer number to specify frequency index, use real to specify frequency directly.")
    elif type(frequency) in [float,np.float32,np.float64]:
      try:
        frequencyIndex = np.where(cdfFile.variables["frequency"][:]==frequency)[0][0]
      except IndexError: 
        raise ValueError("frequency was not found. Use integer number to specify frequency index, use real to specify frequency directly.")
    else:
      raise ValueError("frequency type not understood! Use integer number to specify frequency index, use real to specify frequency directly.")
    
    TB = cdfFile.variables["tb"][:,:,levelIndex,angleIndex,frequencyIndex,stokesIndex]
    shapeTB = np.shape(TB)
    
    xgridArray = np.zeros(shapeTB)
    ygridArray = np.zeros(shapeTB)
    
    if xgrid == "index": #more options could be lat, lon, time....
      xgridArray.T[:] = range(shapeTB[0])
    else:
      raise ValueError("did not understand xgrid")
    
    if ygrid == "index":
      ygridArray[:] = range(shapeTB[1])
    else:
      raise ValueError("did not understand ygrid")
    
    
    cdfFile.close()
  elif "r" in data.__dict__.keys():
    try:
      angleIndex = np.where(angle==data.r["angles_deg"])[0][0]
    except IndexError: 
      raise ValueError("angle was not found")
    if type(frequency) in [int,np.int32,np.int64]:
      if frequency < len(data.set["freqs"]):
        frequencyIndex=frequency
      else:
        raise ValueError("frequencyIndex was not found. Use integer number to specify frequency index, use real to specify frequency directly.")
    elif type(frequency) in [float,np.float32,np.float64]:
      try:
        frequencyIndex = np.where(data.set["freqs"]==frequency)[0][0]
      except IndexError: 
        raise ValueError("frequency was not found. Use integer number to specify frequency index, use real to specify frequency directly.")
    else:
      raise ValueError("frequency type not understood! Use integer number to specify frequency index, use real to specify frequency directly.")
    
    
    #no get the data
    TB = data.r["tb"][:,:,levelIndex,angleIndex,frequencyIndex,stokesIndex]
    shapeTB = np.shape(TB)
    
    xgridArray = np.zeros(shapeTB)
    ygridArray = np.zeros(shapeTB)
    
    if xgrid == "index": #more options could be lat, lon, time....
      xgridArray.T[:] = range(shapeTB[0])
    else:
      raise ValueError("did not understand xgrid")
    
    if ygrid == "index":
      ygridArray[:] = range(shapeTB[1])
    else:
      raise ValueError("did not understand ygrid")
  else:
    raise ValueError("data type not understood!")
  
  if levels[0] == "min":
    TBmin = int(np.min(TB))
  else:
    TBmin = levels[0]
  if levels[1] == "max":
    TBmax = int(np.max(TB))+1
  else:
    TBmax = levels[1]   
  
  if fig == None:
    fig = plt.figure()
  sp = fig.add_subplot(111)
  #print zip(TB.ravel(),range(len(TB.ravel())))
  if (TB.shape[1] == 1):
    sp.plot(range(len(TB.ravel())),TB.ravel())
    sp.set_xlabel(xgrid)
    sp.set_ylabel("TB [K]")#, fontsize=7)
  elif (TB.shape[0] == 1):
    sp.plot(range(len(TB.ravel())),TB.ravel())
    sp.set_xlabel(ygrid)
    sp.set_ylabel("TB [K]")#, fontsize=7)
  else:
    #cf = sp.contourf(xgridArray,ygridArray,TB,extend="both",levels=np.linspace(TBmin,TBmax,101),cmap=cmap)
    cf = sp.pcolor(xgridArray,ygridArray,TB,vmin=TBmin,vmax=TBmax,cmap=cmap)
    cb = plt.colorbar(cf)
    sp.set_xlabel(xgrid)
    sp.set_ylabel(ygrid)  
    
  sp.set_title(title)
  
  return fig
  
  
def plotTBLine(data,polarisation="H",angle=180.,outlevel="space",frequency=0,cmap=None,title="",fig =None,xIndex = "all",yIndex =slice(0,1),xVar = None,freqIndices="all",legendLoc="auto",relative=False):
  """
  plot brightness temperatures from pamtra
  
  input:
    data: netcdf file OR pyPamtra object
    polarisatuion: H or V
    outlevel: space or ground
    frequency: if real number frequency, if int, index of frequency
    xgrid: index or?
    ygrid index or? 
    relative: relative to 0,0?
  """
  #first parse options for both netcdf and pyPamtra
  if polarisation=="H":
    stokesIndex =1
  elif polarisation=="V":
    stokesIndex =0
  else:
    raise ValueError("polarisation must be H or V")

  if outlevel=="space":
    levelIndex =0
  elif outlevel=="ground":
    levelIndex =1
  else:
    raise ValueError("outlevel must be space or ground")
    

  try:
    angleIndex = np.where(angle==data.r["angles_deg"])[0][0]
  except IndexError: 
    raise ValueError("angle was not found")

    
  if xIndex == "all":
    xIndex = slice(0,data._shape2D[0])
  if yIndex == "all":
    yIndex = slice(0,data._shape2D[1])
  #no get the data
  TB = data.r["tb"][xIndex,yIndex,levelIndex,angleIndex,:,stokesIndex]
  if fig == None:
    fig = plt.figure()
  sp = fig.add_subplot(111)
    
  if xVar == None:
    xData = np.arange(len(TB[...,0].ravel()))
  else:
    xData = data.p[xVar][xIndex,yIndex]
  
  if freqIndices == "all": freqIndices = range(len(data.set["freqs"]))
  
  cols = niceColors(len(freqIndices),cmap='hsv_r')
  
  rel = deepcopy(TB[0,0,:])
  if relative == False:
    rel[:] = 0.

  for fn, ff in enumerate(freqIndices):
    sp.plot(xData,TB[...,ff].ravel()-rel[ff],color=cols[fn],label="%.2f GHz"%data.set["freqs"][ff])
  
  sp.set_ylim(np.min(TB-rel),np.max(TB-rel))
  sp.set_xlim(np.min(xData),np.max(xData))
  
  if xVar == None:
      sp.set_xlabel("index")
  else:
      sp.set_xlabel(xVar)
  sp.set_ylabel("Brightness Temperature [K]")

  sp.legend(loc=legendLoc)
  sp.set_title(title)

  return fig