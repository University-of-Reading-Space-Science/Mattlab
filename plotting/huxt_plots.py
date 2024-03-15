# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:06:51 2024

@author: mathewjowens
"""



import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plotspeedmap(vr_map, vr_longs, vr_lats, fig = None, ax = None):
    
    if fig == None:
        fig = plt.figure(figsize = (9,4.5))
    
    if ax == None:
        ax = plt.subplot(1,1,1)
    
    pc = ax.pcolor(vr_longs.value*180/np.pi, vr_lats.value*180/np.pi, vr_map.value, 
               shading='auto',vmin=250, vmax=650)

    ax.set_ylim([-90,90]); ax.set_xlim([0,360])
    ax.set_xlabel('Carrington Longitude [deg]')
    ax.set_ylabel('Latitude [deg]')


    ax.axes.xaxis.set_ticks([0,90,180,270,360])
    ax.axes.yaxis.set_ticks([-90,-45,0,45,90])
    #ax.axes.xaxis.set_ticklabels([])
    #ax.axes.yaxis.set_ticklabels([])
    plt.sca(ax)
    #colorbar
    axins = inset_axes(ax,
                        width="100%",  # width = 50% of parent_bbox width
                        height="10%",  # height : 5%
                        loc='upper right',
                        bbox_to_anchor=(0.28, 0.58, 0.72, 0.5),
                        bbox_transform=ax.transAxes,
                        borderpad=0,)
    cb = fig.colorbar(pc, cax = axins, orientation = 'horizontal',  pad = -0.1)
    cb.ax.xaxis.set_ticks_position("top")
    ax.text(0.15,1.05,r'$V_{SW}$ [km/s]' , 
            fontsize = 11, transform=ax.transAxes, backgroundcolor = 'w')
    
    return fig, ax, axins
