#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 19:29:14 2024

@author: vy902033
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d




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
            
    
 
def periodic_interp(x, xp, fp, period = 2*np.pi, ignore_nans = True):
    """
    Perform interpolation with periodic boundaries.
    
    Parameters:
    - x: array-like, the x-coordinates where to interpolate.
    - xp: 1-D sequence of floats, the x-coordinates of the data points.
    - fp: 1-D sequence of floats, the y-coordinates of the data points.
    
    Returns:
    - array-like, the interpolated values.
    """
    # Remove NaNs from the data
    if ignore_nans:
        valid = ~np.isnan(fp)
        xp = xp[valid]
        fp = fp[valid]
            
    
    # Extend the data to handle periodic boundaries
    extended_xp = np.concatenate((xp - period, xp, xp + period))
    extended_fp = np.concatenate((fp, fp, fp))
    
    # Wrap the x values to be within the range [0, period)
    x_wrapped = np.mod(x, period)
    
    return np.interp(x_wrapped, extended_xp, extended_fp)


def apply_periodic_kernel(values, sigma, period=360):
    """
    Apply a Gaussian kernel with periodic boundaries.
    
    Parameters:
    - values: array-like, the values to which the kernel is applied.
    - sigma: float, the standard deviation of the Gaussian kernel.
    - period: int, the period of the data.
    
    Returns:
    - array-like, the values after applying the kernel.
    """
    # Remove NaNs from the data
    valid = ~np.isnan(values)
    values_cleaned = values[valid]
    
    # Extend the data to handle periodic boundaries
    extended_values = np.concatenate((values_cleaned, values_cleaned, values_cleaned))
    
    # Apply the Gaussian filter
    filtered_extended_values = gaussian_filter1d(extended_values, sigma)
    
    # Extract the middle part, corresponding to the original data length
    start = len(values_cleaned)
    end = 2 * len(values_cleaned)
    filtered_values = filtered_extended_values[start:end]
    
    # Reintroduce NaNs in their original positions
    full_filtered_values = np.full_like(values, np.nan)
    full_filtered_values[valid] = filtered_values
    
    return full_filtered_values
       