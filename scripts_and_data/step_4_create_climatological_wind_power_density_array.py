#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script creates an empty array for the climatological wind power density time series that will be filled in step 5.
'''

import numpy as np
import xarray as xr

import settings as settings
import utilities_for_data_processing as utilities


if settings.use_dask_mpi:
    # Create a Dask cluster and client
    dask_workers_folder_name = 'dask_workers_step_4'
    client = utilities.create_dask_cluster_and_client(dask_workers_folder_name)

# Define the log file name
log_file_name = 'step_4'

utilities.write_to_log_file(log_file_name, 'Task started\n\n', new_file=True)

# Load the wind power density time series of a leap year (1980)
wind_power_density_time_series = xr.open_dataarray(settings.climate_data_directory+'/1980-wind_power_density_time_series.nc', chunks=settings.chunks)

# Change the wind power density time series to make it consistent with the climatological wind power density time series
wind_power_density_time_series['time'] = list(np.arange(366*24))
wind_power_density_time_series = wind_power_density_time_series.rename({'time': 'hour of year'})

# Create an empty array of the climatological wind power density time series
climatological_wind_power_density_time_series = xr.zeros_like(wind_power_density_time_series).rename('Climatological 100m wind power density')
utilities.save(settings.result_directory+'/climatological_wind_power_density_time_series.nc', climatological_wind_power_density_time_series)
utilities.write_to_log_file(log_file_name, 'Climatological wind power calculated\n\n')

utilities.write_to_log_file(log_file_name, 'Task completed\n\n')

if settings.use_dask_mpi:
    utilities.close_dask_cluster_and_client(dask_workers_folder_name, client)
