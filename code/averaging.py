#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 19:29:14 2024

@author: vy902033
"""

import numpy as np




def rolling(mjds, vals, window_days, min_points = 2, method = 'mean'):
    
    #a function to compute rolling estimates of a given parameter, as pandas sucks balls
    #
    #valid methods: mean, median, max, min, sum, var
    
    assert(len(mjds) == len(vals))
    
    smjd = mjds[0]
    fmjd = mjds[-1]
    dt = window_days/2
    
    assert(fmjd - smjd > window_days*2)
    
    rollingvals = vals * np.nan
    for i in range(0, len(mjds)):
        mjd = mjds[i]
        
        #check it's not the start or end of the data sequnce
        if mjd < smjd + dt:
            rollingvals[i] = np.nan
        elif mjd > fmjd - dt:
            rollingvals[i] = np.nan
        else:
            #check there's enough real data
            mask = (mjds > mjd - dt) & (mjds <= mjd +dt)
            npoints = np.count_nonzero(~np.isnan(vals[mask]))
            
            if npoints < min_points:
                rollingvals[i] = np.nan
            else: 
                if (method == 'mean'):
                    rollingvals[i] = np.nanmean(vals[mask])
                elif (method == 'median'):
                    rollingvals[i] = np.nanmedian(vals[mask])
                elif (method == 'max'):
                    rollingvals[i] = np.nanmax(vals[mask])
                elif (method == 'min'):
                    rollingvals[i] = np.nanmin(vals[mask])
                elif (method == 'sum'):
                    rollingvals[i] = np.nansum(vals[mask])
                elif (method == 'var'):
                    rollingvals[i] = np.nanvar(vals[mask])
                else:
                    print('unknown method for rolling')
                    
                    
    return rollingvals
            
    
        