#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script calculates the energy deficit for the seasonal variability.
'''

import xarray as xr

import settings as settings
import utilities_for_data_processing as utilities


if settings.use_dask_mpi:
    # Create a Dask cluster and client
    dask_workers_folder_name = 'dask_workers_step_7'
    client = utilities.create_dask_cluster_and_client(dask_workers_folder_name)

# Define the log file name
log_file_name = 'step_7'

utilities.write_to_log_file(log_file_name, 'Task started\n\n', new_file=True)

# Load the climate data
climatological_wind_power_density_time_series = xr.open_dataarray(settings.result_directory+'/climatological_wind_power_density_time_series.nc', chunks=settings.chunks_climatological)
climatological_wind_power_density_time_series = climatological_wind_power_density_time_series.rename({'hour of year': 'time'})

# Calculate the seasonal variability
energy_deficit_for_seasonal_variability = utilities.calc_energy_deficit(climatological_wind_power_density_time_series, climatological_wind_power_density_time_series, 'constant_target')
utilities.write_to_log_file(log_file_name, 'energy_deficit_for_seasonal_variability '+str(energy_deficit_for_seasonal_variability.nbites/1e6)+'\n\n')
utilities.save(settings.result_directory+'/energy_deficit_for_seasonal_variability.nc', energy_deficit_for_seasonal_variability)
utilities.write_to_log_file(log_file_name, 'Seasonal variability calculated\n\n')

utilities.write_to_log_file(log_file_name, 'Task completed\n\n')

if settings.use_dask_mpi:
    utilities.close_dask_cluster_and_client(dask_workers_folder_name, client)
