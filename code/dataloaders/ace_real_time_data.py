__author__ = 'Harriet Turner'
__email__ = 'h.m.turner@reading.ac.uk'

import os
import glob
import h5py
import numpy as np
import requests
import pandas as pd
from datetime import datetime
from astropy.time import Time


def scrape_ace_swepam_urls(date_start, date_stop):
    """
    Function to scrape the solarsoft website for the URLs of the daily SWEPAM data files from the ACE spacecraft.
    Parameters
    ----------
    date_start: a datetime object
    date_stop: a datetime object

    Returns
    -------
    url_list: a list of the URLs to the daily data files
    """

    # the page that contains all the links to the files
    parent_directory = 'https://sohoftp.nascom.nasa.gov/sdb/goes/ace/daily/'
    # earliest for when data is available
    data_start = datetime(2001, 8, 7)
    # latest for when data is available
    data_end = datetime.now()

    # checking to make sure the date is within the data range
    if date_start < data_start:
        raise ValueError("Date start is outside of data range.")

    if date_stop > data_end:
        raise ValueError("Date stop is in the future and outside of data range.")

    else:
        # list to append the URLs to
        url_list = []

        # list of daily dates between the start and stop times
        days = pd.date_range(date_start, date_stop, freq='24h')

        for day in days:
            str_day = day.strftime("%Y%m%d")
            link = parent_directory + str_day + '_ace_swepam_1m.txt'
            url_list.append(link)

    return url_list


def get_dirs():
    """
    Return a dictionary of relevant directories
    Returns:
    """
    root_dir = os.path.abspath(os.path.dirname(__file__))
    data_dir = os.path.join(root_dir, "data", "in_situ")
    
    #if the data/in_situ dir does not exist, create it
    if not os.path.exists(data_dir):
        # Create the directory
        os.makedirs(data_dir)
    
    dirs = {'root': root_dir, 'data': data_dir}
    return dirs

def update_ace_swepam_url_masterlist(url_list):
    """
    Function to update the master list of downloaded ACE SWEPAM URLs.
    Args:
        url_list: A list of URLs of the ACE SWEPAM Beacon data
    Returns:
    """
    dirs = get_dirs()
    # Write this list of URLs to file
    url_list_file = os.path.join(dirs['data'], 'ace_swepam_url_masterlist.txt')
    with open(url_list_file, 'a') as f:
        for url in url_list:
            f.write(f"{url}\n")

    return


def build_ace_swepam_url_masterlist( date_start = datetime(2019, 1, 1), 
                                     date_stop = datetime(2024, 1, 1) ):
    """
    Function to generate the URLs of the ACE SWEPAM 1m data from 2019-01-01 to present.
    These are saved in ace_swepam_url_masterlist.txt, which is used in forming the full archive of these data.
    Returns:
    """
    url_list = scrape_ace_swepam_urls(date_start, date_stop)
    update_ace_swepam_url_masterlist(url_list)

    return


def load_ace_swepam_url_masterlist():
    """
    Function to load in the list of ACE SWEPAM URLs already downloaded into archive.
    Returns:
    """
    dirs = get_dirs()
    url_list_file = os.path.join(dirs['data'], 'ace_swepam_url_masterlist.txt')
    with open(url_list_file, 'r') as f:
        url_list = f.readlines()

    url_list = [url.strip() for url in url_list]

    return url_list


def download_ace_swepam_url(url):
    """
    Function to download the ACE SWEPAM URL
    Args:
        url: A full URL to the data file to be downloaded
    Returns:
        data_file_path: The path to the downloaded data
    """
    dirs = get_dirs()

    # Download the URL
    print("Downloading: " + url)
    r = requests.get(url, stream='True')
    data_file_path = os.path.join(dirs['data'], url.split('/')[-1])
    with open(data_file_path, 'wb') as f:
        f.write(r.content)

    return data_file_path


def read_ace_swepam_txt_data(data_file_path):
    """
    Function to load the data from the ACE SWEPAM file into a pandas dataframe
    Args:
        data_file_path: Path to the txt file of the daily ACE SWEPAM data
    Returns:
        df: A dataframe of the times and solar wind speeds contained in the txt file. Invalid values set to NaN
        success_flag: A boolean, True if records were extracted from the cdf file.
    """
    # Open the txt file and export the time and solar wind speed
    txt_file = pd.read_fwf(data_file_path, skiprows=18,
                           colspecs=[(0, 4), (4, 7), (7, 10), (10, 14), (14, 16), (16, 24), (24, 32), (32, 37), (37, 48),
                                     (48, 60), (60, 72)],
                           names=['Year', 'Month', 'Day', 'Hour', 'Minute', 'MJD', 'Seconds', 'Status', 'Proton_density',
                                  'Bulk_speed', 'Ion_T'])

    try:
        # Combining the date parts into a datetime object
        dates = pd.to_datetime(txt_file[['Year', 'Month', 'Day', 'Hour', 'Minute']])
        # Extracting the bulk speed
        vsw = txt_file['Bulk_speed']

        # Combining into a new dataframe
        df = pd.DataFrame({'time': dates, 'vsw': vsw})
        # Replacing the invalid values to NaNs
        df['vsw'].mask(df['vsw'] < 0, np.NaN, inplace=True)

        success_flag = True

    except ValueError as ve:
        print(ve)
        print('Error loading {}. No records in file'.format(data_file_path))
        print('DataFrame of NaNs returned')
        df = pd.DataFrame({'time': [], 'vsw': []})
        success_flag = False

    return df, success_flag


def hourly_average_ace_swepam_data(df):
    """
    Function to compute hourly average values from the minute resolution ACE SWEPAM data.
    Args:
        df: A dataframe of the minute resolution solar wind data, as returned by read_ace_swepam_txt_data()
    Returns:
        df_avg: A dataframe of the hourly averaged solar wind data. Invalids set to NaN.
    """
    # Get new range of times to compute the averages over.
    time_min = df['time'].min().strftime("%Y-%m-%dT00:00")
    time_max = df['time'].max().strftime("%Y-%m-%dT23:00")
    times_hourly = pd.date_range(time_min, time_max, freq="1H") + pd.Timedelta(minutes=30)

    # Create the new data frame
    df_avg = pd.DataFrame({'time': times_hourly, 'vsw': np.NaN, 'n_samples': np.NaN})

    # Average the values
    for j in range(df_avg.shape[0]):
        id_vals = (df['time'] > df_avg.loc[j, 'time'] - pd.Timedelta(minutes=30)) & \
                  (df['time'] <= df_avg.loc[j, 'time'] + pd.Timedelta(minutes=30))

        # Only compute average if more than 10 mins of data in hourly mean
        n_samples = np.sum(id_vals)
        df_avg.loc[j, 'n_samples'] = n_samples
        if n_samples > 10:
            df_avg.loc[j, 'vsw'] = df.loc[id_vals, 'vsw'].mean()
        else:
            df_avg.loc[j, 'vsw'] = np.NaN

    # Only keep data in original range
    id_keep = (df_avg['time'] >= df['time'].min()) & (df_avg['time'] <= df['time'].max())
    df_avg = df_avg[id_keep]

    return df_avg


def download_ace_swepam_1min_archive(date_start = datetime(2019, 1, 1), 
                                     date_stop = datetime(2024, 1, 1)):
    """
    Function to download the archive of ACE SWEPAM SW data. Creates the master list of URLs and
    sequentially downloads the files and exports to the HDF5 archive of all these data. Data gaps are not interpolated
    over, and so time series is not necessarily continuous in HDF5 file.
    Returns:
    """
    # Get some needed dirs
    dirs = get_dirs()

    # First, scrape the list of URLs for the needed CDF files
    build_ace_swepam_url_masterlist(date_start, date_stop)

    # Open a HDF5 file to store the results in.
    out_file_name = os.path.join(dirs['data'], 'ace_swepam_1min_archive.h5')
    out_file = h5py.File(out_file_name, 'w')
    out_file.create_dataset('time', (1440,), maxshape=(None,), dtype='S48')  # 1440 minutes in a day hence shape
    out_file.create_dataset('vsw', (1440,), maxshape=(None,), dtype='float64')

    # Get list of all URLS of the STA Beacon data
    url_list = load_ace_swepam_url_masterlist()

    for i, url in enumerate(url_list):

        data_file_tmp = download_ace_swepam_url(url)
        df = read_ace_swepam_txt_data(data_file_tmp)[0]
        
        # Bin these data into the HDF5 file
        if i == 0:
            # First pass we just dump the data
            out_file['time'][:] = df['time'].dt.strftime("%Y-%m-%dT%H:%M:%S").values.astype('S')
            out_file['vsw'][:] = df['vsw'].values
        else:
            # Now we have to resize the data set and append to end
            out_file['time'].resize((out_file['time'].shape[0] + df.shape[0]), axis=0)
            out_file['time'][-df.shape[0]:] = df['time'].dt.strftime("%Y-%m-%dT%H:%M:%S").values.astype('S')

            out_file['vsw'].resize((out_file['vsw'].shape[0] + df.shape[0]), axis=0)
            out_file['vsw'][-df.shape[0]:] = df['vsw'].values

        # Update the file on disk
        out_file.flush()

    # Close out the HDF5 file
    out_file.close()

    # Clean up the not needed CDF files
    clean_up_ace_swepam_files()

    # Compute the 1hr archive while at it
    compute_ace_1hr_archive()
    return


def update_ace_1min_archive():
    """
    Function to update the archive of ACE 1min solar wind speed archive. Updates both the masterlist of downloaded
    URLs and the HDF5 file of data.  Data gaps are not interpolated over, and so time series is not necessarily
    continuous in HDF5 file.
    Returns:
    """
    # Get some needed dirs
    dirs = get_dirs()
    ace_data_path = os.path.join(dirs['data'], "ace_swepam_1min_archive.h5")

    # Load in the hdf5 of output data
    ace_1min = h5py.File(ace_data_path, 'r+')

    # Get time of last data entry
    time_last_entry = Time(ace_1min['time'][-1]).to_datetime()
    # Get start and stop times for URL search. Extend start by 32 days to prior month to ensure date_range behaves.
    date_start = time_last_entry - pd.Timedelta("32d")
    date_stop = datetime.now()

    # Get list of URLs for new epoch
    new_urls = scrape_ace_swepam_urls(date_start, date_stop)

    # Load in masterlist of URLs already processed.
    url_masterlist = load_ace_swepam_url_masterlist()
    # We will need final file for updating incomplete records
    final_file = url_masterlist[-1]

    # Trim the list of new URLs so that only those not in master list are processed.
    new_urls = [url for url in new_urls if url not in url_masterlist]
    # Add previous final file to new_urls, so previously incomplete data gets redownloaded.
    new_urls.insert(0, final_file)

    for url in new_urls:

        # Unzip and get dataframe of minute resolution data
        url_path = download_ace_swepam_url(url)
        df, success_flag = read_ace_swepam_txt_data(url_path)


        if success_flag:
            # Test if df has overlap with HDF5 file.
            t_start_new = df['time'][0].strftime("%Y-%m-%dT%H:%M:%S")
            if t_start_new < ace_1min['time'][-1].decode("utf8"):
                # There is overlap, so find index where it begins
                k = 0
                time_test = False
                while time_test is False:
                    k = k - 1
                    time_test = ace_1min['time'][k].decode("utf8") <= t_start_new

                # Update the HDF5 file with these data
                t_out = df['time'].dt.strftime("%Y-%m-%dT%H:%M:%S").values.astype('S')
                if k == -24:
                    ace_1min['time'][k:] = t_out
                    ace_1min['vsw'][k:] = df['vsw'].values
                else:  # Not sure this clause is strictly necessary? This case might be impossible?
                    ace_1min['time'][k:k + df.shape[0]] = t_out
                    ace_1min['vsw'][k:k + df.shape[0]] = df['vsw'].values

            else:
                # No overlap, so resize HDF5 data and append to end
                ace_1min['time'].resize((ace_1min['time'].shape[0] + df.shape[0]), axis=0)
                ace_1min['time'][-df.shape[0]:] = df['time'].dt.strftime("%Y-%m-%dT%H:%M:%S").values.astype('S')

                ace_1min['vsw'].resize((ace_1min['vsw'].shape[0] + df.shape[0]), axis=0)
                ace_1min['vsw'][-df.shape[0]:] = df['vsw'].values
        else:
            print("No data extracted from {}".format(url_path))


        # Update the file on disk
        ace_1min.flush()

    # Close out the HDF data
    ace_1min.close()

    # Update the masterlist of downloaded URLs
    new_urls.remove(final_file)
    update_ace_swepam_url_masterlist(new_urls)

    # Clean up the not needed CDF files
    clean_up_ace_swepam_files()

    return


def clean_up_ace_swepam_files():
    """
    Remove all the ACE SWEPAM txt files.
    Returns:
    """
    dirs = get_dirs()

    files_to_remove = glob.glob(os.path.join(dirs['data'], '*ace*swepam*1m*.txt'))
    for file in files_to_remove:
        print("Deleting file: " + file)
        os.remove(file)

    return


def compute_ace_1hr_archive():
    """
        Function to compute the 1HR mean archive of the 1min ACE SWEPAM data.
    """

    # Get some needed dirs
    dirs = get_dirs()
    ace_1min_path = os.path.join(dirs['data'], "ace_swepam_1min_archive.h5")
    ace_1hr_path = os.path.join(dirs['data'], "ace_swepam_1hr_archive.h5")

    # Load in the hdf5 of output data
    ace_1min = h5py.File(ace_1min_path, 'r')
    ace_1hr = h5py.File(ace_1hr_path, 'w')

    df = pd.DataFrame({'time': ace_1min['time'].asstr()[:], 'vsw': ace_1min['vsw'][:]})
    df['time'] = pd.to_datetime(df['time'])

    df_avg = hourly_average_ace_swepam_data(df)
    ace_1hr.create_dataset('time', (df_avg.shape[0],), maxshape=(None,), dtype='S48')
    ace_1hr.create_dataset('vsw', (df_avg.shape[0],), maxshape=(None,), dtype='float64')
    ace_1hr.create_dataset('n_samples', (df_avg.shape[0],), maxshape=(None,), dtype='float64')
    ace_1hr['time'][:] = df_avg['time'].dt.strftime("%Y-%m-%dT%H:%M:%S").values.astype('S')
    ace_1hr['vsw'][:] = df_avg['vsw'].values
    ace_1hr['n_samples'][:] = df_avg['n_samples'].values

    ace_1hr.flush()
    ace_1hr.close()
    ace_1min.close()

    return


def update_ace_1hr_archive():
    """
    Function to update the 1HR mean archive of the 1min ACE SWEPAM data.
    """

    # Get some needed dirs
    dirs = get_dirs()
    ace_1min_path = os.path.join(dirs['data'], "ace_swepam_1min_archive.h5")
    ace_1hr_path = os.path.join(dirs['data'], "ace_swepam_1hr_archive.h5")

    # Load in the hdf5 of output data
    ace_1min = h5py.File(ace_1min_path, 'r')
    ace_1hr = h5py.File(ace_1hr_path, 'r+')

    # Find where 1min data begins after final 1hr value
    time_last_1hr_val = Time(ace_1hr['time'][-1].decode("utf8")).to_datetime()
    i = 0
    time_test = False
    while time_test is False:
        i = i - 1
        time_test = Time(ace_1min['time'][i].decode("utf8")).to_datetime() <= time_last_1hr_val

    # Extract only subset of 1min data not averaged into 1hr archive
    df = pd.DataFrame({'time': ace_1min['time'].asstr()[i:], 'vsw': ace_1min['vsw'][i:]})
    df['time'] = pd.to_datetime(df['time'])
    # Average these values up to 1hr
    df_avg = hourly_average_ace_swepam_data(df)

    # Double check no overlap in time
    df_avg = df_avg[df_avg['time'] > time_last_1hr_val]
    if df_avg.shape[0] >= 1:
        # Resize the hdf5 datasets and append to end:
        ace_1hr['time'].resize((ace_1hr['time'].shape[0] + df_avg.shape[0]), axis=0)
        ace_1hr['time'][-df_avg.shape[0]:] = df_avg['time'].dt.strftime("%Y-%m-%dT%H:%M:%S").values.astype('S')

        ace_1hr['vsw'].resize((ace_1hr['vsw'].shape[0] + df_avg.shape[0]), axis=0)
        ace_1hr['vsw'][-df_avg.shape[0]:] = df_avg['vsw'].values

        ace_1hr['n_samples'].resize((ace_1hr['n_samples'].shape[0] + df_avg.shape[0]), axis=0)
        ace_1hr['n_samples'][-df_avg.shape[0]:] = df_avg['n_samples'].values
    else:
        print("ACE SWEPAM hourly archive already up to date")

    ace_1hr.flush()
    ace_1hr.close()
    ace_1min.close()

    return


def update_ace_archive():
    """
    Update the 1 min and 1 hr long ACE archvies
    Returns:
    """
    update_ace_1min_archive()
    update_ace_1hr_archive()
    return


def load_ace_interval(date_start, date_stop, res="hr"):
    """
    Function to load in the ACE real time archive data. Significant data gaps can exist in the series/HDF5
    data.
    Args:
        date_start: datetime of the beginning of the requested interval
        date_stop: datetime of the end of the requested interval
        res: String, "hr" or "min" to load the hourly or minute data. defaults to hourly
    Returns:
        df: A dataframe with the hourly ACE SWEPAM RT data.
    """
    # Get needed dirs
    dirs = get_dirs()

    # Open a HDF5 file to extract archive
    if res == "hr":
        data_file_name = 'ace_swepam_1hr_archive.h5'
    elif res == "min":
        data_file_name = 'ace_swepam_1min_archive.h5'
    else:
        raise ValueError('res should be hr or min, but {} was given'.format(res))

    data_file_name = os.path.join(dirs['data'], data_file_name)
    data_file = h5py.File(data_file_name, 'r')

    # Make pandas dataframe of output and sort out the times.
    if res == "hr":
        df = pd.DataFrame({'time': data_file['time'].asstr()[:], 'vsw': data_file['vsw'][:],
                           'n_samples': data_file['n_samples'][:]})
    elif res == "min":
        df = pd.DataFrame({'time': data_file['time'].asstr()[:], 'vsw': data_file['vsw'][:]})

    df['time'] = pd.to_datetime(df['time'])
    # Close the hdf5 file
    data_file.close()

    if (date_start < df['time'].min()) | (date_start > df['time'].max()):
        print("Warning. Requested interval start date is outside of local ACE SWEPAM archive.")

    if (date_stop < df['time'].min()) | (date_stop > df['time'].max()):
        print("Warning. Requested interval end date is outside of local ACE SWEPAM archive.")

    # Only return data between requested limits
    df = df[(df['time'] >= date_start) & (df['time'] <= date_stop)]

    return df


def bravda_export(date_start, date_stop, output_path=None):
    """
    Function to export ACE SWEPAM RT data in the format required by BRaVDA.
    Args:
        date_start: A date time object giving the start of the window.
        date_stop: A date time object giving the end of the window.
    Returns:
    """

    df = load_ace_interval(date_start, date_stop, res='hr')
    df['year'] = df['time'].dt.year
    df['doy'] = df['time'].dt.dayofyear
    df['hour'] = df['time'].dt.hour
    df.drop(columns='time', inplace=True)
    cols = ['year', 'doy', 'hour', 'vsw']

    output_basename = 'ACE_rt_observations.txt'
    if output_path is None:
        output_name = output_basename
    else:
        output_name = os.path.join(output_path, output_basename)

    formatters = {"year": "{:d}".format, "doy": "{:d}".format, "hour": "{:d}".format, "vsw": "{:3.2f}".format}
    with open(output_name, 'w') as out_file:
        out_file.write(df.to_string(columns=cols, header=False, index=False, na_rep='nan', formatters=formatters))

    return


def load_ace_rt_vsw(date_start = datetime(2024, 2, 1), 
                    date_stop = datetime(2024,3, 1)):
    """
    a wrapper to download the ACE real time data, average to 1-hour and return
    a dataframe with datetime and vsw.

    """

    #delete any existing url masterlist
    dirs = get_dirs()
    url_list_file_name = os.path.join(dirs['data'], "ace_swepam_url_masterlist.txt")
    if os.path.exists(url_list_file_name):
        # Delete the file
        os.remove(url_list_file_name)
    
    #download the data for this interval
    download_ace_swepam_1min_archive( date_start,  date_stop )
    #average to 1 hour.
    compute_ace_1hr_archive() 
    
    #return the dataframe for this interval
    df = load_ace_interval(date_start, date_stop, res="hr")
    return df