#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script calculates the wind speed and the annual mean wind speed for a given year.

    It reads the u and v components of the wind from the ERA5 dataset, calculates the variables, and saves the results in a NetCDF file.
'''

import numpy as np
import xarray as xr
import sys

import settings as settings
import utilities_for_data_processing as utilities


# Read the year to be processed
year = int(sys.argv[1])

if settings.use_dask_mpi:
    # Create a Dask cluster and client
    dask_workers_folder_name = 'dask_workers_step_1_'+str(year)
    client = utilities.create_dask_cluster_and_client(dask_workers_folder_name)

# Define the log file name
log_file_name = 'step_1_'+str(year)

utilities.write_to_log_file(log_file_name, 'Year '+str(year)+'\n\n', new_file=True, write_time=False)
utilities.write_to_log_file(log_file_name, 'Task started\n\n')

# Load the climate data for the year
u_component_time_series = xr.open_dataarray(settings.climate_data_directory+'/ERA5-{:d}-100m_u_component_of_wind.nc'.format(year), chunks=settings.chunks)
v_component_time_series = xr.open_dataarray(settings.climate_data_directory+'/ERA5-{:d}-100m_v_component_of_wind.nc'.format(year), chunks=settings.chunks)

# Calculate the wind speed
wind_speed_time_series = np.sqrt(np.power(u_component_time_series,2)+np.power(v_component_time_series,2)).rename('100m wind speed')
utilities.save(settings.climate_data_directory+'/{:d}-wind_speed_time_series.nc'.format(year), wind_speed_time_series)
utilities.write_to_log_file(log_file_name, 'Wind speed calculated\n\n')

# Calculate the annual mean wind speed
annual_mean_wind_speed = wind_speed_time_series.mean(dim='time').rename('100m annual mean wind speed')
annual_mean_wind_speed = annual_mean_wind_speed.expand_dims(dim={'year': 1})
annual_mean_wind_speed['year'] = np.array([year])
utilities.save(settings.result_directory+'/{:d}-annual_mean_wind_speed.nc'.format(year), annual_mean_wind_speed)
utilities.write_to_log_file(log_file_name, 'Annual wind speed calculated\n\n')

utilities.write_to_log_file(log_file_name, 'Task completed\n\n')

if settings.use_dask_mpi:
    utilities.close_dask_cluster_and_client(dask_workers_folder_name, client)
