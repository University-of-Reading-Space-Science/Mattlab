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


import conversions.helio_time as htime
import conversions.helio_coords as hcoords
import plotting.huxt_plots as hplots
import system.system as system


import huxt_inputs as Hin


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

