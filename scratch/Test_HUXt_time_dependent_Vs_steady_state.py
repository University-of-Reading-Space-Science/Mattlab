#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 10:07:40 2024

@author: vy902033
"""

#a comparison of time-dependent and steady state HUXt solutions

import os
import datetime
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import pandas as pd



import huxt as H
import huxt_inputs as Hin
import huxt_analysis as HA


startdate = datetime.datetime(2000,1,1,0)
v_boundary = np.ones(128) * 400 * (u.km/u.s)
v_boundary[30:50] = 600 * (u.km/u.s)
v_boundary[95:125] = 700 * (u.km/u.s)



#steady-state solution
cr, cr_lon_init = Hin.datetime2huxtinputs(startdate)

model_ss = H.HUXt(v_boundary=v_boundary, cr_num=cr, cr_lon_init=cr_lon_init,
                simtime=7*u.day, r_min=21.5*u.solRad, r_max=250*u.solRad,
                frame='synodic', lon_start=0 * u.rad,
                lon_stop=0.1 * u.rad, dt_scale =4)
model_ss.solve([])
HA.plot_timeseries(model_ss, 215*u.solRad, lon=0.0*u.rad)




#time-dependent solution

vlongs=[]
mjds = []
for n in range(0,7):
    vlongs.append(v_boundary)
    mjds.append(n+51544.0)

vlongs = np.array(vlongs).T
mjds = np.array(mjds)

#set up a HUXt run with this boundary condition
model_td = Hin.set_time_dependent_boundary(vgrid_Carr=vlongs, time_grid=mjds,
                                        starttime=startdate, simtime=7*u.day, 
                                        r_min=21.5*u.solRad, r_max=250*u.solRad,
                                        frame='synodic', lon_start=0 * u.rad,
                                        lon_stop=0.1 * u.rad, dt_scale=4)
model_td.solve([])

HA.plot_timeseries(model_td, 215*u.solRad, lon=0.0*u.rad)


# <codecell>


import os
import datetime
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import astropy.units as u
import pandas as pd


import conversions.helio_time as htime
import conversions.helio_coords as hcoords
import plotting.huxt_plots as hplots
import system.system as system

import huxt as H
import huxt_inputs as Hin
import huxt_analysis as HA



#set plot defaults
system.plot_defaults()

datadir = os.path.join(os.getenv('DBOX'),'Data','WSA_MO_API')

ndays = 1  # download coronal solutions up to this many days prior to the forecast date

firstdate = datetime.datetime(2023,1,1,0)
finaldate = datetime.datetime(2023,2,1,0)

deacc = True #Whether to deaccelerate the WSA maps froms 215 to 21.5 rS
plot_map_frames = False

daysec = 24*60*60

bufferdays_td = 4 #number of days before the forecasttime to start the time-dependent runs

# <codecell>



thisdate = firstdate

count = 1
vlong_list = []
brlong_list = []
mjd_list = []
while thisdate <= finaldate:
    
    year = str(thisdate.year)
    month = str(thisdate.month)
    day = str(thisdate.day)
    
    #create the expected filename
    filename = 'models%2Fenlil%2F' + year +'%2F' + month +'%2F' + day + '%2F0%2Fwsa.gong.fits'
    filepath = os.path.join(datadir, filename)
    
    #get the MJD
    mjd = htime.datetime2mjd(thisdate)
    mjd_list.append(mjd)
    
    #get Earth lat
    Ecoords = hcoords.carringtonlatlong_earth(mjd)
    E_lat = np.pi/2 - Ecoords[0][0]
    
    
    
    if os.path.exists(filepath):
        vr_map, vr_longs, vr_lats, br_map, br_longs, br_lats, cr_fits \
            = Hin.get_WSA_maps(filepath)
            
        # deaccelerate the WSA map from 1-AU calibrated speeds to expected 21.5 rS values
        if deacc:
            vr_map_deacc = vr_map.copy()
            for nlat in range(1, len(vr_lats)):
                vr_map_deacc[nlat, :], lon_temp = Hin.map_v_inwards(vr_map[nlat, :], 215 * u.solRad,
                                                                    vr_longs, 21.5* u.solRad)
            vr_map = vr_map_deacc
            
        E_lats = E_lat * (vr_longs.value * 0 + 1)
        
            
        #get the Earth lat slice
        v_in = Hin.get_WSA_long_profile(filepath, lat=E_lat*u.rad)
        if deacc:
            # deaccelerate them?
            v_in, lon_temp = Hin.map_v_inwards(v_in, 215 * u.solRad, vr_longs,  21.5 * u.solRad)
        br_in = Hin.get_WSA_br_long_profile(filepath, lat=E_lat*u.rad)
         
        #store the data
        vlong_list.append(v_in)
        brlong_list.append(br_in)
        
    else:
        print(filepath + '; not found')
    
    
    
    #advance the date
    thisdate = thisdate + datetime.timedelta(days=ndays)
    count = count + 1
  
# <codecell>    
  
#convert to arrays    
vlongs = np.array(vlong_list).T
brlongs = np.array(brlong_list).T
mjds = np.array(mjd_list)
time = htime.mjd2datetime(mjds)

#plot the CarrLon - time Vin
fig = plt.figure(figsize = (10,10))
ax = plt.subplot(2,1,1)
pc = ax.pcolor(time, vr_longs.value*180/np.pi, vlongs, 
            shading='auto',vmin=250, vmax=650)
ax.set_ylabel('Carrington Longitude [deg]')
ax.axes.yaxis.set_ticks([0,90,180,270,360])
ax.text(0.15,1.05,r'$V_{SW}$ [km/s]' , 
        fontsize = 11, transform=ax.transAxes, backgroundcolor = 'w')
cbar = plt.colorbar(pc, ax=ax)

ax = plt.subplot(2,1,2)
pc = ax.pcolor(time, br_longs.value*180/np.pi, brlongs, 
            shading='auto')
ax.set_ylabel('Carrington Longitude [deg]')
ax.axes.yaxis.set_ticks([0,90,180,270,360])
ax.text(0.15,1.05,r'$Br [nT]' , 
        fontsize = 11, transform=ax.transAxes, backgroundcolor = 'w')
cbar = plt.colorbar(pc, ax=ax)


# <codecell> Steady-state model run

forecasttime = datetime.datetime(2023,1,17)

cr, cr_lon_init = Hin.datetime2huxtinputs(forecasttime)

#find the map with this date
id_t = np.argmin(np.abs(time - forecasttime))
#set up a HUXt run with this boundary condition
model_ss = H.HUXt(v_boundary=vlongs[:,id_t]*u.km/u.s, cr_num=cr, cr_lon_init=cr_lon_init,
               simtime=7*u.day, dt_scale=4, r_min = 21.5*u.solRad, lon_out=0.0*u.rad)
model_ss.solve([])

HA.plot_timeseries(model_ss, 215*u.solRad, lon=0.0*u.rad)


# <codecell> time-dependent model run

f_mjd = htime.datetime2mjd(forecasttime)
#runstart = forecasttime - datetime.timedelta(days = bufferdays_td)

#get the CR num and cr_lon_init that this corresponds to
cr, cr_lon_init = Hin.datetime2huxtinputs(forecasttime)

#find the map with this date
id_t = np.argmin(np.abs(time - forecasttime))

#create a input carr_v with the last few days and then the current value repeated into the future
vlongs_slice = []
mjds_slice = []

for i in range(0,10):
    vlongs_slice.append(vlongs[:,id_t])  
    mjds_slice.append(f_mjd + i)
    

vlongs_slice = np.array(vlongs_slice).T
mjds_slice = np.array(mjds_slice)


#set up a HUXt run with this boundary condition
model_td = Hin.set_time_dependent_boundary(vgrid_Carr=vlongs_slice, time_grid=mjds_slice,
                                        starttime=forecasttime, simtime=7*u.day, 
                                        r_min=21.5*u.solRad, r_max=250*u.solRad,
                                        frame='synodic', lon_start=0 * u.rad,
                                        lon_stop=0.1 * u.rad,
                                        dt_scale=4)
model_td.solve([])
HA.plot_timeseries(model_td, 215*u.solRad, lon=0.0*u.rad)