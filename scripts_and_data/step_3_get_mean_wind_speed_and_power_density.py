#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script combines the annual mean wind speed and annual mean wind power density for all years into a single file.

    It also calculates the minimum annual mean wind power density, the mean wind speed, and the mean wind power density.
'''

import numpy as np
import xarray as xr

import settings as settings
import utilities_for_data_processing as utilities


if settings.use_dask_mpi:
    # Create a Dask cluster and client
    dask_workers_folder_name = 'dask_workers_step_3'
    client = utilities.create_dask_cluster_and_client(dask_workers_folder_name)

# Define the log file name
log_file_name = 'step_3'

utilities.write_to_log_file(log_file_name, 'Task started\n\n', new_file=True)

# Load the climate data for all years
annual_mean_wind_speed = xr.open_mfdataset(settings.result_directory+'/*-annual_mean_wind_speed.nc', chunks=settings.chunks, parallel=True)
annual_mean_wind_power_density = xr.open_mfdataset(settings.result_directory+'/*-annual_mean_wind_power_density.nc', chunks=settings.chunks, parallel=True)

# Save the annual mean wind speeds and annual mean wind power densities in a single file
utilities.save(settings.result_directory+'/annual_mean_wind_speed.nc', annual_mean_wind_speed)
utilities.save(settings.result_directory+'/annual_mean_wind_power_density.nc', annual_mean_wind_power_density)
utilities.write_to_log_file(log_file_name, 'Annual values calculated\n\n')

# Calculate the minimum annual mean wind speed
min_annual_mean_wind_power_density = annual_mean_wind_power_density.min(dim='year')
utilities.save(settings.result_directory+'/min_annual_mean_wind_power_density.nc', min_annual_mean_wind_power_density)
utilities.write_to_log_file(log_file_name, 'Minimum values calculated\n\n')

# Calculate the start and end year
start_year = int(annual_mean_wind_power_density['year'][0])
end_year = int(annual_mean_wind_power_density['year'][-1])

# Calculate a series of weights where each weight is the number of hours in a year
weights = []
for year in range(start_year,end_year+1):
    if year%4 == 0:
        weights.append(366*24)
    else:
        weights.append(365*24) 
weights = xr.DataArray(data=weights/np.sum(weights), dims=['year'], coords={'year': annual_mean_wind_speed['year']}, name='weights')

# Calculate the mean wind speed
mean_wind_speed = (annual_mean_wind_speed['100m annual mean wind speed']*weights).sum(dim='year').rename('100m mean wind speed')
utilities.save(settings.result_directory+'/mean_wind_speed.nc', mean_wind_speed)
utilities.write_to_log_file(log_file_name, 'Mean wind speed calculated\n\n')

# Calculate the mean wind power density
mean_wind_power_density = (annual_mean_wind_power_density['100m annual mean wind power density']*weights).sum(dim='year').rename('100m mean wind power density')
utilities.save(settings.result_directory+'/mean_wind_power_density.nc', mean_wind_power_density)
utilities.write_to_log_file(log_file_name, 'Mean wind power density calculated\n\n')

utilities.write_to_log_file(log_file_name, 'Task completed\n\n')

if settings.use_dask_mpi:
    utilities.close_dask_cluster_and_client(dask_workers_folder_name, client)
