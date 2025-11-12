# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 16:04:46 2025

@author: mathe
"""

import requests

import os
from datetime import datetime

import sunpy.timeseries as ts
import pandas as pd


def get_SO_PAS_LL_data(StartTimeSel, EndTimeSel, folder = None):
    
    if folder is None:
        DSROOTFOLDER = os.path.join(os.getenv('DBOX'), 'Data', 'SolarOrbiter', 'LL')
    else:
        DSROOTFOLDER = folder
    
    
    # Selection criteria
    InstrumentSel = 'SWA'
    #StartTimeSel = datetime(2025, 2, 10)
    #EndTimeSel = datetime(2025, 3, 30)
    CreationTimeSel = datetime(2020, 2, 1)
    ProcLevelSel = 'LL02'
    FileNameSel = 'solo_LL02_swa-pas-mom'
    SOARProductType = 'LOW_LATENCY'
    
    # Download file list
    print('Downloading file list from SOAR...')
    URL = 'http://soar.esac.esa.int/soar-sl-tap/tap/sync'
    params = {
        'FORMAT': 'JSON',
        'REQUEST': 'doQuery',
        'LANG': 'ADQL',
        'QUERY': f"SELECT * FROM v_public_files WHERE instrument='{InstrumentSel}' AND archived_on >'2020-01-01'"
    }
    response = requests.get(URL, params=params)
    response.raise_for_status()
    data = response.json()
    
    # Filter files
    print('Files which qualify for download')
    FilePick = {'dataid': [], 'filename': []}
    NSel = 0
    
    for record in data['data']:
        CreationTimeHere = datetime.strptime(record[0], '%Y-%m-%dT%H:%M:%S.%f')
        StartTimeHere = datetime.strptime(record[1], '%Y-%m-%dT%H:%M:%S.%f')
    
        if (CreationTimeHere > CreationTimeSel and
            StartTimeSel < StartTimeHere < EndTimeSel and
            record[5] == InstrumentSel and
            FileNameSel in record[3] and
            record[8] == ProcLevelSel):
            
            print(f"{record[1]} | {record[2]} | {record[3]} | {record[5]} | {record[8]} made {record[0]}")
            NSel += 1
            FilePick['dataid'].append(record[6])
            FilePick['filename'].append(record[3])
    
    print(f'There are {NSel} files to download')
    
    # Download files
    SOARURL = 'http://soar.esac.esa.int/soar-sl-tap/data'
    
    for dataid, filename in zip(FilePick['dataid'], FilePick['filename']):
        print(f'Downloading {filename} from SOAR')
        params = {
            'product_type': SOARProductType,
            'retrieval_type': 'PRODUCT',
            'data_item_id': dataid
        }
        response = requests.get(SOARURL, params=params)
        response.raise_for_status()
        
        file_path = os.path.join(DSROOTFOLDER, filename)
        with open(file_path, 'wb') as file:
            file.write(response.content)
    
        print(f'Downloaded {file_path}')
        
        
    # Convert CDF files to DataFrame
    def cdf_to_dataframe(directory):
        dataframes = []
        for file in os.listdir(directory):
            if file.endswith('.cdf'):
                file_path = os.path.join(directory, file)
                try:
                    ts_data = ts.TimeSeries(file_path)
                    df = ts_data.to_dataframe().reset_index()  # Ensure time is a column
                    df['filename'] = file
                    dataframes.append(df)
                except Exception as e:
                    print(f'Failed to read {file}: {e}')
        return pd.concat(dataframes, ignore_index=True) if dataframes else pd.DataFrame()
    
    # Load all CDF files into a single DataFrame
    print('Converting CDF files to DataFrame...')
    dataframe = cdf_to_dataframe(DSROOTFOLDER)
    print('Combined DataFrame:')
    #print(dataframe.head())
    
    # Sort dataframe by EPOCH
    dataframe = dataframe.sort_values(by='EPOCH').reset_index(drop=True)
    
    # Resample data to 1-hour intervals, timestamp at the middle of the window
    numeric_cols = dataframe.select_dtypes(include=['number']).columns
    dataframe = dataframe.set_index('EPOCH').resample('1H', origin='epoch', closed='left', label='right')[numeric_cols].mean()
    dataframe.index = dataframe.index - pd.Timedelta(minutes=30)
    dataframe = dataframe.reset_index()
    #print(dataframe.head())
    
    
    df = dataframe[['EPOCH', 'SWA_PAS_VELOCITY_RTN_0', 'SWA_PAS_DENSITY']]  # Keep only the relevant columns
    
    # Rename the columns
    df.rename(columns={'EPOCH': 'time', 
                       'SWA_PAS_VELOCITY_RTN_0': 'vsw', 
                       'SWA_PAS_DENSITY': 'n'}, inplace=True)
    
    return df

def get_SO_MAG_LL_data(StartTimeSel, EndTimeSel, folder = None):
    
    if folder is None:
        DSROOTFOLDER = os.path.join(os.getenv('DBOX'), 'Data', 'SolarOrbiter', 'LL')
    else:
        DSROOTFOLDER = folder
    
    
    # Selection criteria
    InstrumentSel = 'MAG'
    #StartTimeSel = datetime(2025, 2, 10)
    #EndTimeSel = datetime(2025, 3, 30)
    CreationTimeSel = datetime(2020, 1, 1)
    ProcLevelSel = 'LL01'
    FileNameSel = 'solo_LL01_mag'
    SOARProductType = 'LOW_LATENCY'
    
    # Download file list
    print('Downloading file list from SOAR...')
    URL = 'http://soar.esac.esa.int/soar-sl-tap/tap/sync'
    params = {
        'FORMAT': 'JSON',
        'REQUEST': 'doQuery',
        'LANG': 'ADQL',
        'QUERY': f"SELECT * FROM v_public_files WHERE instrument='{InstrumentSel}' AND archived_on >'2020-01-01'"
  }
    response = requests.get(URL, params=params)
    response.raise_for_status()
    data = response.json()
    
    # Filter files
    print('Files which qualify for download')
    FilePick = {'dataid': [], 'filename': []}
    NSel = 0
    
    for record in data['data']:
        CreationTimeHere = datetime.strptime(record[0], '%Y-%m-%dT%H:%M:%S.%f')
        StartTimeHere = datetime.strptime(record[1], '%Y-%m-%dT%H:%M:%S.%f')
    
        if (CreationTimeHere > CreationTimeSel and
            StartTimeSel < StartTimeHere < EndTimeSel and
            record[5] == InstrumentSel and
            FileNameSel in record[3] and
            record[8] == ProcLevelSel):
            
            print(f"{record[1]} | {record[2]} | {record[3]} | {record[5]} | {record[8]} made {record[0]}")
            NSel += 1
            FilePick['dataid'].append(record[6])
            FilePick['filename'].append(record[3])
    
    print(f'There are {NSel} files to download')
    
    # Download files
    SOARURL = 'http://soar.esac.esa.int/soar-sl-tap/data'
    
    for dataid, filename in zip(FilePick['dataid'], FilePick['filename']):
        print(f'Downloading {filename} from SOAR')
        params = {
            'product_type': SOARProductType,
            'retrieval_type': 'PRODUCT',
            'data_item_id': dataid
        }
        response = requests.get(SOARURL, params=params)
        response.raise_for_status()
        
        file_path = os.path.join(DSROOTFOLDER, filename)
        with open(file_path, 'wb') as file:
            file.write(response.content)
    
        print(f'Downloaded {file_path}')
        
        
    # Convert CDF files to DataFrame
    def cdf_to_dataframe(directory):
        dataframes = []
        for file in os.listdir(directory):
            if file.endswith('.cdf'):
                file_path = os.path.join(directory, file)
                try:
                    ts_data = ts.TimeSeries(file_path)
                    df = ts_data.to_dataframe().reset_index()  # Ensure time is a column
                    df['filename'] = file
                    dataframes.append(df)
                except Exception as e:
                    print(f'Failed to read {file}: {e}')
        return pd.concat(dataframes, ignore_index=True) if dataframes else pd.DataFrame()
    
    # Load all CDF files into a single DataFrame
    print('Converting CDF files to DataFrame...')
    dataframe = cdf_to_dataframe(DSROOTFOLDER)
    print('Combined DataFrame:')
    #print(dataframe.head())
    
    # Sort dataframe by EPOCH
    dataframe = dataframe.sort_values(by='EPOCH').reset_index(drop=True)
    
    # Resample data to 1-hour intervals, timestamp at the middle of the window
    numeric_cols = dataframe.select_dtypes(include=['number']).columns
    dataframe = dataframe.set_index('EPOCH').resample('1H', origin='epoch', closed='left', label='right')[numeric_cols].mean()
    dataframe.index = dataframe.index - pd.Timedelta(minutes=30)
    dataframe = dataframe.reset_index()
    #print(dataframe.head())
    
    
    df = dataframe[['EPOCH', 'SWA_PAS_VELOCITY_RTN_0', 'SWA_PAS_DENSITY']]  # Keep only the relevant columns
    
    # Rename the columns
    df.rename(columns={'EPOCH': 'time', 
                       'SWA_PAS_VELOCITY_RTN_0': 'vsw', 
                       'SWA_PAS_DENSITY': 'n'}, inplace=True)
    
    return df

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    StartTimeSel = datetime(2025, 2, 10)
    EndTimeSel = datetime(2025, 2, 12)
    df_pas = get_SO_PAS_LL_data(StartTimeSel, EndTimeSel, folder = None)
    df_mag = get_SO_MAG_LL_data(StartTimeSel, EndTimeSel, folder = None)
    
    plt.figure()
    plt.plot(df_pas['time'], df_pas['vsw'])

