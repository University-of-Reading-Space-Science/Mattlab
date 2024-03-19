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
import imageio
from moviepy.editor import ImageSequenceClip


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

firstdate = datetime.datetime(2022,5,1,0)
finaldate = datetime.datetime(2024,3,1,0)

deacc = True #Whether to deaccelerate the WSA maps froms 215 to 21.5 rS
input_res_days = 0.1

daysec = 24*60*60


bufferdays_td = 4 #number of days before the forecasttime to start the time-dependent runs

single_run_now = False
single_anim_now = False
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
  
 
  
#convert to arrays    
vlongs_1d = np.array(vlong_list).T
brlongs_1d = np.array(brlong_list).T
mjds_1d = np.array(mjd_list)
time_1d = htime.mjd2datetime(mjds_1d)
n_longs = len(vlongs_1d[:,0])


#increase the time resolution of the vlongs
mjds = np.arange(mjds_1d[0], mjds_1d[-1], input_res_days)
time = htime.mjd2datetime(mjds)
vlongs = np.ones((n_longs, len(mjds)))
brlongs = np.ones((n_longs, len(mjds)))
for n in range(0, n_longs):
    vlongs[n,:] = np.interp(mjds, mjds_1d, vlongs_1d[n,:])
    brlongs[n,:] = np.interp(mjds, mjds_1d, brlongs_1d[n,:])


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


# <codecell> run HUXt with the time-dependent boundary condition
if single_run_now:
    runstart = firstdate#datetime.datetime(2023,1,1)
    smjd = htime.datetime2mjd(runstart)
    
    runend = finaldate#datetime.datetime(2023,12,31)
    fmjd = htime.datetime2mjd(runend)
    
    simtime = (runend-runstart).days * u.day
      
    
    #set up the model, with (optional) time-dependent bpol boundary conditions
    model_td = Hin.set_time_dependent_boundary(vlongs, mjds, runstart, simtime, 
                                            r_min=21.5 *u.solRad, r_max=250*u.solRad,
                                            #bgrid_Carr = brlongs.T, 
                                            dt_scale=10, latitude=0*u.deg,
                                            frame = 'synodic')
    
    #trace a bunch of field lines from a range of evenly spaced Carrington longitudes
    dlon = (20*u.deg).to(u.rad).value
    lon_grid = np.arange(dlon/2, 2*np.pi-dlon/2 + 0.0001, dlon)*u.rad
    
    #give the streakline footpoints (in Carr long) to the solve method
    model_td.solve([], streak_carr = lon_grid)
    
    HA.plot(model_td, 12*u.day)
    
    
    #get the Earth timeseries
    #get conditions at Earth
    td_Earth_ts = HA.get_observer_timeseries(model_td, observer='Earth')
    #convert the model time to MJD
    td_tim_mjd = model_td.time_init.mjd + model_td.time_out.value/daysec
    
    if single_anim_now:

        HA.animate(model_td, tag='HUXt_WSA_2023_time_dependent', duration = 60, fps = 20) # This takes about two minutes.


    


    #generate side-by-side images of steady-state and time-dependent runs
    #=====================================================================
    #delete existing images
    # List all files in the directory
    # Directory containing the PNG files
    png_dir = os.path.join(datadir,'frames')
    files = os.listdir(png_dir)
    
    # Iterate through the files and delete .png files
    for file in files:
        if file.endswith('.png'):
            os.remove(os.path.join(png_dir, file))
        
    
    count = 1
    fracs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    
    plot_mjd = smjd
    while plot_mjd <= fmjd - 1:
        plot_date = htime.mjd2datetime(plot_mjd).item()
        dmjd = plot_mjd - smjd
        
        dang = 2 * np.pi * dmjd / (model_td.rotation_period.value/daysec)
        plotlongrid = H._zerototwopi_(lon_grid + dang*u.rad) * u.rad
        
        #run the steady-state model
        #===========================
        #get the CR num and cr_lon_init that this corresponds to
        cr, cr_lon_init = Hin.datetime2huxtinputs(plot_date)
        #find the map with this date
        id_t = np.argmin(np.abs(time - plot_date))
        #set up a HUXt run with this boundary condition
        model_ss = H.HUXt(v_boundary=vlongs[:, id_t]*u.km/u.s, cr_num=cr, cr_lon_init=cr_lon_init,
                       simtime=1*u.day, dt_scale=4, r_min = 21.5*u.solRad, frame='synodic')
        model_ss.solve([], streak_carr = plotlongrid)
        
        
        for plot_day_frac in fracs:
        
            fig = plt.figure(figsize=(12, 7))
            fig.subplots_adjust(left=0.01, bottom=0.17, right=0.99, top=0.99)
            
            ax = plt.subplot(121,  projection='polar')
            HA.plot(model_ss, plot_day_frac*u.day, 
                    fighandle=fig, axhandle=ax, plotHCS=False, annotateplot = False)
            ax.set_title('Steady-state', fontsize =16)
            
            fig.legend(ncol=5, loc='lower center', frameon=False, handletextpad=0.2, columnspacing=1.0)
            
            ax = plt.subplot(122,  projection='polar')   
            HA.plot(model_td, (dmjd + plot_day_frac)*u.day, 
                    fighandle=fig, axhandle=ax, plotHCS=False,  annotateplot = False)
            ax.set_title('Time-dependent', fontsize =16)
            
            label = plot_date.strftime('%Y-%m-%d')
            fig.text(0.45, 0.85, label, fontsize=16)
            
            
            
            
            
            #save the plot
            formatted_number = "{:04d}".format(count)
            fig_image = os.path.join(datadir,'frames', 'huxt' + formatted_number + '.png')
            plt.savefig(fig_image)
            
            plt.close('all')
            
            count = count + 1
        
        plot_mjd = plot_mjd + 1.0
    
    
    
    # List the PNG files
    png_files = sorted([f for f in os.listdir(png_dir) if f.endswith('.png')])
    
    # Read PNG files into a list of images
    images = [imageio.imread(os.path.join(png_dir, f)) for f in png_files]
    
    # Create a video clip from the images
    clip = ImageSequenceClip(images, fps=24)  # Change fps as needed
    
    # Write the video clip to an MP4 file
    output_file = os.path.join(datadir,'frames',"output_detail.mp4")
    clip.write_videofile(output_file, codec='libx264', fps=10) 
# <codecell> put together 1 to 7-day advance forecasts for steady state models

# Initialize lists for for1d to for7d
for_lists = [ [] for _ in range(1, 8) ]

# Initialize lists for tim1d to tim7d
tim_lists = [ [] for _ in range(1, 8) ]

startdate = firstdate + datetime.timedelta(days=4) #datetime.datetime(2023,1,3,0)
stopdate = finaldate#datetime.datetime(2024,1,3,0)

forecasttime = startdate
while forecasttime <=stopdate:

    #get the CR num and cr_lon_init that this corresponds to
    cr, cr_lon_init = Hin.datetime2huxtinputs(forecasttime)
    
    #find the map with this date
    id_t = np.argmin(np.abs(time - forecasttime))
    #set up a HUXt run with this boundary condition
    model = H.HUXt(v_boundary=vlongs[:, id_t]*u.km/u.s, cr_num=cr, cr_lon_init=cr_lon_init,
                   simtime=7*u.day, dt_scale=4, r_min = 21.5*u.solRad, lon_out=0.0*u.rad)
    model.solve([])
    

    #convert the model time to MJD
    tim_mjd = model.time_init.mjd + model.time_out.value/daysec
    f_mjd = htime.datetime2mjd(forecasttime)
    
    #get conditions at Earth
    Earth_ts = HA.get_observer_timeseries(model, observer='Earth')
    
    #now hack out days of data for the various forecasts
    for d_advanced in range(0,7):
        # mask = (model.time_out.value >= daysec*d_advanced) & (model.time_out.value < daysec*(d_advanced+1))
        # tim_lists[d_advanced].extend(tim_mjd[mask])
        # for_lists[d_advanced].extend(Earth_ts.loc[mask,'vsw'])
        
        mask = (tim_mjd >= f_mjd + d_advanced) & (tim_mjd < f_mjd + d_advanced + 1)
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






#find the maximum period of overlap
smjd = tim_lists[6][0]
fmjd = tim_lists[0][-1]

mask_omni = (data['mjd'] >= smjd) & (data['mjd'] <= fmjd)
data = data[mask_omni]

mae_ss = []

for d_advanced in range(0,7):
    #interpolate the forecast onto the omni time step
    data['for'+str(d_advanced)] = np.interp(data['mjd'], tim_lists[d_advanced], 
                                            for_lists[d_advanced])
    
    mae_ss.append(np.nanmean(np.abs(data['for'+str(d_advanced)] - data['V'])))
    
print(mae_ss)

plt.figure()
plt.plot(htime.mjd2datetime(np.array(tim_lists[0])), for_lists[0])
plt.plot(data['datetime'], data['V'])


# <codecell> do a synodic time-dependent run

# #set up the model, with (optional) time-dependent bpol boundary conditions
# model = Hin.set_time_dependent_boundary(vlongs.T, mjds, runstart, simtime, 
#                                         r_min=r_min, r_max=250*u.solRad,
#                                         frame='synodic', lon_start=0 * u.rad,
#                                         lon_stop=0.1 * u.rad,
#                                         dt_scale=10, latitude=0*u.deg)



# #give the streakline footpoints (in Carr long) to the solve method
# model.solve([])

# #get the Earth timeseries
# #get conditions at Earth
# td_Earth_ts = HA.get_observer_timeseries(model, observer='Earth')
# #convert the model time to MJD
# td_tim_mjd = model.time_init.mjd + model.time_out.value/daysec


    
# data['fortd'] = np.interp(data['mjd'], td_tim_mjd, td_Earth_ts['vsw'])  
# mae = np.nanmean(np.abs(data['fortd'] - data['V']))
# print(mae)

# plt.figure()
# plt.plot(htime.mjd2datetime(np.array(td_tim_mjd)), td_Earth_ts['vsw'])
# plt.plot(data['datetime'], data['V'])



# <codecell> put together 1 to 7-day advance forecasts for time-dependent solar wind

# Initialize lists for for1d to for7d
for_lists_td = [ [] for _ in range(1, 8) ]

# Initialize lists for tim1d to tim7d
tim_lists_td = [ [] for _ in range(1, 8) ]

startdate = firstdate + datetime.timedelta(days = bufferdays_td)#datetime.datetime(2023,1,3,0)
stopdate = finaldate#datetime.datetime(2024,1,3,0)

forecasttime = startdate
while forecasttime <=stopdate:
    
    runstart = forecasttime - datetime.timedelta(days = bufferdays_td)

    #get the CR num and cr_lon_init that this corresponds to
    cr, cr_lon_init = Hin.datetime2huxtinputs(runstart)
    
    #find the map with this date
    id_t_start = np.argmin(np.abs(time - runstart))
    id_t_stop = np.argmin(np.abs(time - forecasttime))
    f_mjd = htime.datetime2mjd(forecasttime)
    
    #create a input carr_v with the last few days
    vlongs_slice =[]
    mjds_slice = []
    for n in range(id_t_start, id_t_stop):
        vlongs_slice.append(vlongs[:,n])
        mjds_slice.append(mjds[n])
    
    #then project forward using the current value
    for i in range(0,7):
        vlongs_slice.append(vlongs[:,id_t_stop])
        mjds_slice.append(f_mjd + i)
        
    vlongs_slice = np.array(vlongs_slice).T
    mjds_slice = np.array(mjds_slice)
    
    #set up a HUXt run with this boundary condition
    model = Hin.set_time_dependent_boundary(vgrid_Carr=vlongs_slice, time_grid=mjds_slice,
                                            starttime=runstart, simtime=(7+bufferdays_td)*u.day, 
                                            r_min=21.5*u.solRad, r_max=250*u.solRad,
                                            frame='synodic', lon_start=0 * u.rad,
                                            lon_stop=0.1 * u.rad,
                                            dt_scale=10, latitude=0*u.deg)
    model.solve([])
    

    #convert the model time to MJD
    tim_mjd_td = model.time_init.mjd + model.time_out.value/daysec
    f_mjd = htime.datetime2mjd(forecasttime)
    
    #get conditions at Earth
    Earth_ts = HA.get_observer_timeseries(model, observer='Earth')
    
    #now hack out days of data for the various forecasts
    for d_advanced in range(0,7):
        mask = (tim_mjd_td >= f_mjd + d_advanced) & (tim_mjd_td < f_mjd + d_advanced + 1)
        tim_lists_td[d_advanced].extend(tim_mjd_td[mask])
        for_lists_td[d_advanced].extend(Earth_ts.loc[mask,'vsw'])
        
    #advance the date
    forecasttime = forecasttime + datetime.timedelta(days=1)

print('=====================')    
mae_td=[]
for d_advanced in range(0,7):
    #interpolate the forecast onto the omni time step
    data['for_td'+str(d_advanced)] = np.interp(data['mjd'], tim_lists_td[d_advanced], 
                                            for_lists_td[d_advanced])
    
    mae_td.append(np.nanmean(np.abs(data['for_td'+str(d_advanced)] - data['V'])))
    
print(mae_td)

plt.figure()
plt.plot(htime.mjd2datetime(np.array(tim_lists_td[0])), for_lists_td[0])
plt.plot(data['datetime'], data['V'])


plt.figure()
plt.plot(np.arange(1,8,1), mae_ss, 'k',label =  'Steady state')
plt.plot(np.arange(1,8,1), mae_td, 'r',label =  'Time dependent')
plt.legend()
plt.xlabel('Forecast lead time [days]')
plt.ylabel('MAE [km/s]')
