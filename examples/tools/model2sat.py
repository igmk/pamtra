import numpy as np
import satConv
from matplotlib import pyplot as plt
from matplotlib import colors, colorbar
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import pickle
import imp
import pdb

imp.reload(satConv)

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

def makePlot(fig, ax_pos_grid0, ax_pos_grid1,sc,scm,df,row, norm, cmap='jet',saveFig=False):

    freqs = ['50.3', '89.0', '157.0']
    olabels = ['(a)','(b)','(c)']
    mlabels = ['(d)','(e)','(f)']
    slabels = ['(g)','(h)','(i)']
    dlabels = ['(j)','(k)','(l)']

    ax = fig.add_subplot(ax_pos_grid0[row,0])

    map = Basemap(projection='merc',llcrnrlat=50.,urcrnrlat=60.,llcrnrlon=-42.5,urcrnrlon=-27.5,lat_ts=55.,resolution='c')

    map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,1])
    map.drawmeridians(np.arange(-180.,171.,10.),labels=[1,0,0,1])               

    X,Y = map(sc._reduce(sc.c['lon']),sc._reduce(sc.c['lat']))

    map.contourf(X,Y,sc._reduce(sc.c['obsTB']),25,norm=norm,cmap=cmap)
    if row == 0: plt.title('Observation')
    ax.text(-0.2, 0.5, freqs[row]+' GHz', fontsize=14, va='center', rotation=90,transform=ax.transAxes)
    ax.annotate(olabels[row],(0.01,0.9), xycoords='axes fraction')

    axm = fig.add_subplot(ax_pos_grid0[row,1])

    map = Basemap(projection='merc',llcrnrlat=50.,urcrnrlat=60.,llcrnrlon=-42.5,urcrnrlon=-27.5,lat_ts=55.,resolution='c')

    map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,1])
    map.drawmeridians(np.arange(-180.,171.,10.),labels=[1,0,0,1])               

    X,Y = map(sc._reduce(scm.c['lon']),sc._reduce(scm.c['lat']))

    map.contourf(X,Y,scm._reduce(scm.c['simTB']),25,norm=norm,cmap=cmap)
    if row == 0: plt.title('Mie')
    axm.annotate(mlabels[row],(0.01,0.9), xycoords='axes fraction')

    axs = fig.add_subplot(ax_pos_grid0[row,2])

    map = Basemap(projection='merc',llcrnrlat=50.,urcrnrlat=60.,llcrnrlon=-42.5,urcrnrlon=-27.5,lat_ts=55.,resolution='c')

    map.drawparallels(np.arange(-80.,81.,10.),labels=[1,0,0,1])
    map.drawmeridians(np.arange(-180.,171.,10.),labels=[1,0,0,1])               

    X,Y = map(sc._reduce(sc.c['lon']),sc._reduce(sc.c['lat']))

    map.contourf(X,Y,sc._reduce(sc.c['simTB']),25,norm=norm,cmap=cmap)
    if row == 0: plt.title('SSRGA')
    axs.annotate(slabels[row],(0.01,0.9), xycoords='axes fraction')

    axdf = fig.add_subplot(ax_pos_grid1[row,:])

    sns.violinplot(x='scat',y='tb',hue='type',split=True,data=df,inner="quartile")
    axdf.yaxis.set_label_position("right")
    axdf.yaxis.tick_right()
    axdf.set_xlabel('')
    axdf.set_ylabel('Brightness Temperature [K]')
    if row == 2:
        axdf.legend(loc='lower left')
        axdf.xaxis.set_ticklabels(['Mie', 'SSRGA'])
    if row != 2:
        axdf.xaxis.set_ticklabels(['',''])
        axdf.get_legend().set_visible(False)
   
    axdf.annotate(dlabels[row],(0.01,0.9), xycoords='axes fraction')
    
    if row == 2:
        axc_box = ax.get_position()
        axc1_box = axs.get_position()
        cax = fig.add_axes([axc_box.x0, axc_box.y0+0.04, axc1_box.x1-axc_box.x0, 0.02], rasterized=True)
        cb1 = colorbar.ColorbarBase(cax,norm=norm,orientation='horizontal',cmap='jet')

        cb1.set_label('Brightness Temperature [K]')

    if saveFig:
        fig.savefig('ecmwf-pamtra_mhs.png')

    return

def makeDataFrames(sc,scm):
    obsSer = pd.Series(sc.c['obsTB'],index=range(len(sc.c['obsTB'])))
    obsDF = pd.DataFrame(obsSer)
    obsDF.insert(loc = 1, column = 'type', value = np.repeat(['Obs'],len(sc.c['obsTB'])))

    simSerMie = pd.Series(scm.c['simTB'],index=range(len(scm.c['simTB'])))
    simMieDF = pd.DataFrame(simSerMie)
    simMieDF.insert(loc = 1, column = 'type', value = np.repeat(['Sim'],len(sc.c['simTB'])))

    osDF = pd.concat([obsDF,simMieDF])

    simSerSSRG = pd.Series(sc.c['simTB'],index=range(len(sc.c['simTB'])))
    simSSRGDF = pd.DataFrame(simSerSSRG)
    simSSRGDF.insert(loc = 1, column = 'type', value = np.repeat(['Sim'],len(sc.c['simTB'])))

    osDF = pd.concat([osDF, obsDF, simSSRGDF])

    osDF.insert(loc=2,column='scat',value = np.concatenate((np.repeat('Mie',len(sc.c['obsTB'])*2),np.repeat('SSRGA',len(sc.c['obsTB'])*2))))
    osDF.rename(columns={0:'tb'},inplace=True)

    return osDF


saveFig = [False,False,False]

makeConv = True

fig = plt.figure(figsize=(9,9))

gs1 = gridspec.GridSpec(3,3)
gs1.update(bottom=0.05, top=0.95,left=0.1, right=0.60, hspace=-0.4)#, wspace=0.1, hspace= 0.1)
gs2 = gridspec.GridSpec(3,2)
gs2.update(bottom=0.15, top=0.85,left=0.65, right=0.9, hspace=0.15)

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
        sc[i], scm[i] = performConvolution('data/ecmwf_ssrga.nc','data/ecmwf_mie.nc','data/'+satFile[i],sensor[i],simChan[i],satChan[i])
        # pdb.set_trace()

    with open('data/conv_ssrg_v2.pkl', 'wb') as handle:
        pickle.dump(sc, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open('data/conv_mie_v2.pkl', 'wb') as handle:
        pickle.dump(scm, handle, protocol=pickle.HIGHEST_PROTOCOL)

else:
    # Work on existing convolution
    with open('data/conv_ssrg_v2.pkl', 'rb') as handle:
        sc = pickle.load(handle)

    with open('data/conv_mie_v2.pkl', 'rb') as handle:
        scm = pickle.load(handle)

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

norm = colors.Normalize(vmin=np.nanmin(cmin), vmax=np.nanmax(cmax))
# cmap = sns.cubehelix_palette(light=1, as_cmap=True)
for i,fstr in enumerate(fstrs):
    df = makeDataFrames(sc[i],scm[i])
    makePlot(fig, gs1, gs2, sc[i],scm[i],df,i,norm,saveFig=saveFig[i])