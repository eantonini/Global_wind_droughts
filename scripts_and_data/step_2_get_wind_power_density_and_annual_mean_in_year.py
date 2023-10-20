#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script calculates the wind power density and the annual mean wind power density for a given year.

    It reads the wind speed calculated in step 1, the surface pressure, and the 2m temperature from the ERA5 dataset, calculates the variables, and saves the results in a NetCDF file.
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
    dask_workers_folder_name = 'dask_workers_step_2_'+str(year)
    client = utilities.create_dask_cluster_and_client(dask_workers_folder_name)

# Define the log file name
log_file_name = 'step_2_'+str(year)

utilities.write_to_log_file(log_file_name, 'Year '+str(year)+'\n\n', new_file=True, write_time=False)
utilities.write_to_log_file(log_file_name, 'Task started\n\n')

# Load the climate data for the year
wind_speed_time_series = xr.open_dataarray(settings.climate_data_directory+'/{:d}-wind_speed_time_series.nc'.format(year), chunks=settings.chunks)
pressure_time_series = xr.open_dataarray(settings.climate_data_directory+'/ERA5-{:d}-surface_pressure.nc'.format(year), chunks=settings.chunks)
temperature_time_series = xr.open_dataarray(settings.climate_data_directory+'/ERA5-{:d}-2m_temperature.nc'.format(year), chunks=settings.chunks)

# Calculate the wind power density
R = 287.05 # Specific gas constant for dry air
wind_power_density_time_series = (0.5*np.power(wind_speed_time_series,3)*pressure_time_series/(R*temperature_time_series)).rename('100m wind power density')
utilities.save(settings.climate_data_directory+'/{:d}-wind_power_density_time_series.nc'.format(year), wind_power_density_time_series)
utilities.write_to_log_file(log_file_name, 'Wind power density calculated\n\n')

# Calculate the annual mean wind power density
annual_mean_wind_power_density = wind_power_density_time_series.mean(dim='time').rename('100m annual mean wind power density')
annual_mean_wind_power_density = annual_mean_wind_power_density.expand_dims(dim={'year': 1})
annual_mean_wind_power_density['year'] = np.array([year])
utilities.save(settings.result_directory+'/{:d}-annual_mean_wind_power_density.nc'.format(year), annual_mean_wind_power_density)
utilities.write_to_log_file(log_file_name, 'Annual wind power density calculated\n\n')

utilities.write_to_log_file(log_file_name, 'Task completed\n\n')

if settings.use_dask_mpi:
    utilities.close_dask_cluster_and_client(dask_workers_folder_name, client)
