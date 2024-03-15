# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 09:01:36 2024

@author: mathewjowens
"""

import os
import datetime
import time

import huxt_tools.huxt_ensembles as Hens


datadir = os.path.join(os.getenv('DBOX'),'Data','WSA_MO_API')

ndays = 1  # download coronal solutions up to this many days prior to the forecast date
ddays = 1  # look for WSA solutions this many days apart
sleep_time = 2

firstdate = datetime.datetime(2024,1,1,0)
finaldate = datetime.datetime(2025,1,1,0)

thisdate = firstdate

requestdate_list = []
success_list = []
mapfilepath_list = []
modeltime_list = []
while thisdate <= finaldate:
    
    sdate = thisdate - datetime.timedelta(days=ndays)
    fdate = thisdate
    
    # get WSA solution
    success, mapfilepath, model_time = Hens.getMetOfficeWSA(sdate, fdate, datadir)
    
    requestdate_list.append(thisdate)
    success_list.append(success)
    mapfilepath_list.append(mapfilepath)
    modeltime_list.append(model_time)
    
    if success == False:
        print('No data for ' + thisdate.strftime("%Y-%m-%d"))
        
    time.sleep(sleep_time)
        
    #advance the date
    thisdate = thisdate + datetime.timedelta(days=ddays)