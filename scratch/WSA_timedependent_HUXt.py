# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 10:14:03 2024

@author: mathewjowens
"""


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
finaldate = datetime.datetime(2024,1,1,0)

deacc = True #Whether to deaccelerate the WSA maps froms 215 to 21.5 rS
plot_map_frames = False

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
        
        if plot_map_frames:
            #plot it
            fig, ax, axins = hplots.plotspeedmap(vr_map, vr_longs, vr_lats)
            
            ax.plot(vr_longs*180/np.pi,E_lats*180/np.pi,'k--',label = 'Earth')
            ax.plot(vr_longs*180/np.pi,E_lats*0,'k')
            ax.legend()
            ax.text(0.95,-0.13, thisdate.strftime("%Y-%m-%d"),  
                     fontsize = 14, transform=ax.transAxes)
            
            #save map image
            formatted_number = "{:03d}".format(count)
            map_image = os.path.join(datadir,'frames', 'frame' + formatted_number + '.png')
            plt.savefig(map_image)
            
            plt.close(fig)
            
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
vlongs = np.array(vlong_list)
brlongs = np.array(brlong_list)
mjds = np.array(mjd_list)
time = htime.mjd2datetime(mjds)

#plot the CarrLon - time Vin
fig = plt.figure(figsize = (10,10))
ax = plt.subplot(2,1,1)
pc = ax.pcolor(time, vr_longs.value*180/np.pi, vlongs.T, 
            shading='auto',vmin=250, vmax=650)
ax.set_ylabel('Carrington Longitude [deg]')
ax.axes.yaxis.set_ticks([0,90,180,270,360])
ax.text(0.15,1.05,r'$V_{SW}$ [km/s]' , 
        fontsize = 11, transform=ax.transAxes, backgroundcolor = 'w')
cbar = plt.colorbar(pc, ax=ax)

ax = plt.subplot(2,1,2)
pc = ax.pcolor(time, br_longs.value*180/np.pi, brlongs.T, 
            shading='auto')
ax.set_ylabel('Carrington Longitude [deg]')
ax.axes.yaxis.set_ticks([0,90,180,270,360])
ax.text(0.15,1.05,r'$Br [nT]' , 
        fontsize = 11, transform=ax.transAxes, backgroundcolor = 'w')
cbar = plt.colorbar(pc, ax=ax)


# <codecell> run HUXt with the time-dependent boundary condition

runstart = datetime.datetime(2023,1,1)
runend = datetime.datetime(2023,12,31)
simtime = (runend-runstart).days * u.day
r_min = 21.5 *u.solRad



#set up the model, with (optional) time-dependent bpol boundary conditions
model = Hin.set_time_dependent_boundary(vlongs.T, mjds, runstart, simtime, 
                                        r_min=r_min, r_max=250*u.solRad,
                                        #bgrid_Carr = brlongs.T, 
                                        dt_scale=10, latitude=0*u.deg,)

#trace a bunch of field lines from a range of evenly spaced Carrington longitudes
dlon = (20*u.deg).to(u.rad).value
lon_grid = np.arange(dlon/2, 2*np.pi-dlon/2 + 0.0001, dlon)*u.rad

#give the streakline footpoints (in Carr long) to the solve method
model.solve([], streak_carr = lon_grid)

HA.plot(model, 12*u.day)

HA.animate(model, tag='HUXt_WSA_2023_time_dependent', duration = 60, fps = 20) # This takes about two minutes.


# <codecell> put together 1 to 7-day advance forecasts

# Initialize lists for for1d to for7d
for_lists = [ [] for _ in range(1, 8) ]

# Initialize lists for tim1d to tim7d
tim_lists = [ [] for _ in range(1, 8) ]

startdate = datetime.datetime(2023,1,3,0)
stopdate = datetime.datetime(2024,1,3,0)

forecasttime = startdate
while forecasttime <=stopdate:

    #get the CR num and cr_lon_init that this corresponds to
    cr, cr_lon_init = Hin.datetime2huxtinputs(forecasttime)
    
    #find the map with this date
    id_t = np.argmin(np.abs(time - forecasttime))
    #set up a HUXt run with this boundary condition
    model = H.HUXt(v_boundary=vlongs[id_t,:]*u.km/u.s, cr_num=cr, cr_lon_init=cr_lon_init,
                   simtime=7*u.day, dt_scale=4, r_min = 21.5*u.solRad, lon_out=0.0*u.rad)
    model.solve([])
    
    daysec = 24*60*60
    #convert the model time to MJD
    tim_mjd = model.time_init.mjd + model.time_out.value/daysec
    
    #get conditions at Earth
    Earth_ts = HA.get_observer_timeseries(model, observer='Earth')
    
    #now hack out days of data for the various forecasts
    for d_advanced in range(0,7):
        mask = (model.time_out.value >= daysec*d_advanced) & (model.time_out.value < daysec*(d_advanced+1))
        tim_lists[d_advanced].extend(tim_mjd[mask])
        for_lists[d_advanced].extend(Earth_ts.loc[mask,'vsw'])
        
    #advance the date
    forecasttime = forecasttime + datetime.timedelta(days=1)
    

# <codecell> Now compare to 1-hour omni data

from sunpy.net import Fido
from sunpy.net import attrs
from sunpy.timeseries import TimeSeries

# Download the 1hr OMNI data from CDAweb
trange = attrs.Time(startdate, stopdate)
dataset = attrs.cdaweb.Dataset('OMNI2_H0_MRG1HR')
result = Fido.search(trange, dataset)
downloaded_files = Fido.fetch(result)

# Import the OMNI data
omni = TimeSeries(downloaded_files, concatenate=True)
data = omni.to_dataframe()

# Set invalid data points to NaN
id_bad = data['V'] == 9999.0
data.loc[id_bad, 'V'] = np.NaN

# Create a datetime column
data['datetime'] = data.index
data['mjd'] = htime.datetime2mjd(data['datetime'].to_numpy())


plt.figure()
plt.plot(htime.mjd2datetime(np.array(tim_lists[0])), for_lists[0])

plt.plot(data['datetime'], data['V'])


plt.figure()
plt.plot(htime.mjd2datetime(np.array(tim_lists[6])), for_lists[6])

plt.plot(data['datetime'], data['V'])

    
    

