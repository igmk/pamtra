# -*- coding: utf-8 -*-


from __future__ import division
import numpy as np
import datetime
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager 
import pyDia2vel

plt.figure()

d = np.linspace(1e-6,1e-4,50)
#name="pavlos_cloud"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)

#d = np.arange(1e-6,6e-3,1e-6)
#name="khvorostyanov01_drops"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)


#d = np.arange(1e-6,1e-3,1e-6)
#name="heymsfield10_particles_hexPlates"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)

#d = np.arange(1e-6,6e-3,1e-6)
#name="khvorostyanov01_particles_hexPlates"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)
#d = np.arange(1e-6,6e-3,1e-6)
#name="heymsfield10_particles_hexPlates"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)

v = 297*(d*100)**0.86
v = v/100.

plt.plot(d,v,label="fit")

d = np.linspace(1e-6,1e-4,10)
name="heymsfield10_particles_MPACE"
v = pyDia2vel.pydia2vel(name,d)
plt.plot(d,v,label=name)

d = np.linspace(1e-6,1e-4,10)
name="heymsfield10_particles_MPACE2"
v = pyDia2vel.pydia2vel(name,d)
plt.plot(d,v,label=name)
#d = np.arange(1e-6,1e-3,1e-6)
#name="heymsfield10_particles_rosettes"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)

#d = np.arange(1e-6,1e-3,1e-6)
#name="khvorostyanov01_particles_rosettes"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)


#d = np.arange(1e-6,6e-3,1e-6)
#name="khvorostyanov01_particles_rosettes"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)

##d = np.arange(1e-6,6e-3,1e-6)
##name="heymsfield10_particles_rosettes"
##v = pyDia2vel.pydia2vel(name,d)
##plt.plot(d,v,label=name)

#d = np.linspace(1e-6,1e-4,50)
#name="khvorostyanov01_particles_aggregates"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)

#d = np.arange(1e-6,6e-3,1e-6)
#name="heymsfield10_particles_aggregates"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)

#d = np.linspace(1e-6,1e-4,50)
#name="khvorostyanov01_particles_sector"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)


#d = np.linspace(1e-6,1e-4,50)
#name="heymsfield10_particles_sector"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)

#d = np.arange(1e-6,6e-3,1e-6)
#name="heymsfield10_particles_sector"
#v = pyDia2vel.pydia2vel(name,d)
##plt.plot(d,v,label=name)

#d = np.arange(1e-6,6e-3,1e-6)
#name="khvorostyanov01_particles_cosmoIce"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)

##d = np.arange(1e-6,6e-3,1e-6)
##name="heymsfield10_particles_cosmoIce"
##v = pyDia2vel.pydia2vel(name,d)
##plt.plot(d,v,label=name)

#d = np.arange(1e-6,6e-3,1e-6)
#name="khvorostyanov01_particles_cosmoSnow"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)

##d = np.arange(1e-6,6e-3,1e-6)
##name="khvorostyanov01_spheres"
##v = pyDia2vel.pydia2vel(name,d)
##plt.plot(d,v,label=name)

##d = np.arange(1e-6,6e-3,1e-6)
##name="khvorostyanov01_graupel"
##v = pyDia2vel.pydia2vel(name,d)
##plt.plot(d,v,label=name)

##d = np.arange(1e-5,6e-3,1e-6)
##name="khvorostyanov01_hail"
##v = pyDia2vel.pydia2vel(name,d)
##plt.plot(d,v,label=name)

##d = np.arange(1e-6,6e-3,1e-6)
##name="rogers_graupel"
##v = pyDia2vel.pydia2vel(name,d)
##plt.plot(d,v,label=name)

##d = np.arange(1e-4,5.8e-3,1e-5)
##name="foote69_rain"
##v = pyDia2vel.pydia2vel(name,d)
##plt.plot(d,v,label=name)

 
##d = np.arange(0.109e-3,6e-3,1e-5)
##name="metek_rain"
##v = pyDia2vel.pydia2vel(name,d)
##plt.plot(d,v,label=name)


##d = np.arange(0,6e-3,1e-6)
##name="rogers_drops"
#v = pyDia2vel.pydia2vel(name,d)
#plt.plot(d,v,label=name)


plt.legend(loc="upper left",prop={'size':8})

#plt.loglog()


plt.xlabel("dim")
plt.ylabel("V")


#d = np.arange(0,8.5,0.01)
#lam = 4.7
#theta = exp(-d/lam)+(1-exp(-d/lam))*(1/(1+(d/lam)))

#figure()
#plot(d,theta)
