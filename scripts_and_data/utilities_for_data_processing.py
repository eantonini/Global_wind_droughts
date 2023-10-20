#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains some utility functions used to process the climate data such as calculate the energy deficits, write to log files, save data, and create and close Dask clusters and clients.
'''

import numpy as np
import xarray as xr
import os
import shutil
from datetime import datetime

import settings as settings


def write_to_log_file(filename, message, new_file=False, write_time=True):
    '''
    Write a message to a log file.

    Parameters
    ----------
    filename : str
        Name of the log file without the extension
    message : str
        Message to write to the log file
    new_file : bool, optional
        If True, the log file is created
    write_time : bool, optional
        If True, the current time is written before the message
    '''
    
    # Create the log file if it does not exist.
    if not os.path.exists(settings.working_directory+'/log_files'):
        os.makedirs(settings.working_directory+'/log_files')
    
    # Determine whether to append or overwrite the log file.
    mode = 'w' if new_file else 'a'

    # Write the message to the log file.
    with open(settings.working_directory+'/log_files/'+filename+'.log', mode) as output_file:
        if write_time:
            # Write the current time to the log file.
            now = datetime.now()
            prefix_time = now.strftime('%H:%M:%S') + ' - '
            output_file.write(prefix_time + message)
        else:
            output_file.write(message)


def save(output_data_filename, data_to_save, compute=True):
    '''
    Save data to a NetCDF file.
    
    Parameters
    ----------
    output_data_filename : str
        Name of the output NetCDF file
    data_to_save : xarray.DataArray or xarray.Dataset
        Data to save
    compute : bool, optional
        If True, the dask array is computed before saving the data
    '''
    
    # Create the directory where the results will be saved if it does not exist.
    if not os.path.exists(settings.result_directory):
        os.mkdir(settings.result_directory)
    
    # Append the data to the existing file if it exists.
    if os.path.exists(output_data_filename):
        with xr.open_dataarray(output_data_filename, chunks=settings.chunks) as original_data:
            dim_of_interest = ('time' if 'time' in original_data.dims else 'year')
            data_to_save = xr.concat([original_data,data_to_save], dim=dim_of_interest)
            data_to_save = data_to_save.sortby(dim_of_interest)

    # Save the data to a NetCDF file.
    data_to_save.to_netcdf(output_data_filename, compute=compute)


def create_dask_cluster_and_client(dask_workers_folder_name):
    '''
    Create a Dask cluster and client.

    Parameters
    ----------
    dask_workers_folder_name : str
        Name of the folder where the Dask workers are stored
    
    Returns
    -------
    client : dask.distributed.Client
        Dask client
    '''
    
    from dask_mpi import initialize
    from dask.distributed import Client

    # Initialize the Dask cluster and client.
    initialize(local_directory=settings.working_directory+'/'+dask_workers_folder_name, memory_limit=1)

    # Create the Dask client.
    client = Client()
    
    return client


def close_dask_cluster_and_client(dask_workers_folder_name, client):
    '''
    Close a Dask cluster and client.
    
    Potential issue: https://github.com/dask/distributed/issues/7192

    Parameters
    ----------
    dask_workers_folder_name : str
        Name of the folder where the Dask workers are stored
    client : dask.distributed.Client
        Dask client
    '''
    
    # Close the Dask client and cluster.
    client.retire_workers()
    client.shutdown()
    client.close()

    # Delete the Dask workers folder.
    shutil.rmtree(settings.working_directory+'/'+dask_workers_folder_name)


def calc_energy_deficit_value(generation, demand):
    '''
    Calculate the energy deficit value given the generation and demand time series.

    Parameters
    ----------
    generation : xarray.DataArray
        Generation time series (longitude x latitude x time)
    demand : xarray.DataArray
        Demand time series (longitude x latitude x time)

    Returns
    -------
    energy_deficit : xarray.DataArray
        Energy deficit time series (longitude x latitude)
    '''

    # Calculate the time series of net generation and change the time dimension to a range of integers. The actual values of the time dimension are not important.
    net_generation = (generation - demand)
    net_generation['time'] = np.arange(len(net_generation['time']))

    # Create a copy of the net generation time series and append it to the original time series.
    second_net_generation = net_generation.copy()
    second_net_generation['time'] = len(net_generation['time']) + np.arange(len(net_generation['time']))
    doubled_net_generation = xr.concat([net_generation,second_net_generation], dim='time').chunk(settings.chunks)
    
    # Calculate the stored energy time series by integrating the net generation time series along the time dimension.
    stored_energy = doubled_net_generation.cumsum(dim='time')

    # Calculate the previous and next values of the stored energy time series.
    stored_energy_next_value = stored_energy.roll(time=-1)
    stored_energy_previous_value = stored_energy.roll(time=1)
    
    # Calculate the maximum and minimum values of the stored energy time series.
    max_indexes = np.logical_and(stored_energy>stored_energy_previous_value,stored_energy>stored_energy_next_value)
    max_values = xr.where(max_indexes, stored_energy, stored_energy*np.nan)
    min_indexes = np.logical_and(stored_energy<stored_energy_previous_value,stored_energy<stored_energy_next_value)
    min_values = xr.where(min_indexes, stored_energy, stored_energy*np.nan)
    
    # Find the maximum value of the stored energy time series that precedes any value in the stored energy time series.
    previous_max_values = max_values.shift(time=1).rolling(time=len(max_values.time), min_periods=1).max()

    # Calculate the energy deficit value as the largest difference between any minumum value of the stored energy time series and its preceding maximum value.
    energy_deficit = (previous_max_values - min_values).max(dim='time')
   
    return energy_deficit


def calc_energy_deficit(time_series, normalization_value, target_equal_to_constant_or_climatological, year=None, climatological_time_series=None, overbuild_factor=1.0, demand_reduction_factor=1.0):
    '''
    Calculate the energy deficit by setting up the generation and demand time series according to the given time series and normalization values.

    Parameters
    ----------
    time_series : xarray.DataArray
        Time series that can be either the wind power density time series of the given year or the climatological wind power density time series (longitude x latitude x time)
    normalization_value : xarray.DataArray
        Normalization value that can be the wind power density time series of the given year, the climatological wind power density time series (longitude x latitude x time), or the minumun annual mean wind power density (longitude x latitude)
    target_equal_to_constant_or_climatological : str
        If 'constant_target', the demand time series is set to a constant value; if 'climatological_target', the demand time series is set to the normalized climatological wind power density time series
    year : int, optional
        Year of the wind power density time series
    climatological_time_series : xarray.DataArray, optional
        Climatological wind power density time series (longitude x latitude x time)
    overbuild_factor : float, optional
        Overbuild factor.
    demand_reduction_factor : float, optional
        Demand reduction factor

    Returns
    -------
    energy_deficit : xarray.DataArray
        Energy deficit values (longitude x latitude)
    '''
    
    # Calculate the normalization value. If the normalization value is a time series, the mean along the time dimension is calculated.
    if 'time' in normalization_value.dims:
        normalization_value = normalization_value.mean(dim='time')
    normalization_value = xr.where(normalization_value > 0, normalization_value, 1.0)
    
    # Calculate the generation time series.
    generation = time_series/normalization_value*overbuild_factor
    
    # Calculate the demand time series.
    if target_equal_to_constant_or_climatological == 'constant_target':
        # Set the demand time series to a constant value.
        demand = xr.ones_like(generation)*demand_reduction_factor
    elif target_equal_to_constant_or_climatological == 'climatological_target':
        # Set the demand time series to the normalized climatological wind power density time series.
        mean_climatological_time_series = climatological_time_series.mean(dim='time')
        demand = xr.where(mean_climatological_time_series > 0, climatological_time_series/mean_climatological_time_series*demand_reduction_factor, demand_reduction_factor)
    
    # Calculate the energy deficit.
    energy_deficit = calc_energy_deficit_value(generation, demand)
    energy_deficit = xr.where(normalization_value > 0, energy_deficit, len(generation['time']))
    
    # Add the year dimension if it is not None.
    if year is not None:
        energy_deficit = energy_deficit.expand_dims(dim={'year': 1})
        energy_deficit['year'] = np.array([year])

    return energy_deficit
