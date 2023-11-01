#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script downloads the ERA5 data.
'''

import cdsapi

import settings as settings

# Create a CDS API client.
c = cdsapi.Client()

# Define the variables to be downloaded.
ERA5_variable_names = ['2m_temperature', 'surface_pressure', '100m_u_component_of_wind', '100m_v_component_of_wind']

# Loop over the years and variables and download the data.
for year in range(1979,2023):
    for ERA5_variable_name in ERA5_variable_names:
        data_path = settings.climate_data_directory+'/ERA5-{:d}-'.format(year)+ERA5_variable_name+'.nc'
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': ERA5_variable_name,
                'year': str(int(year)),
                'month': [str(int(x)) for x in range(1,13)],
                'day': [str(int(x)) for x in range(1,32)],
                'time': [str(int(x))+':00' if x>=10 else '0'+str(int(x))+':00' for x in range(0,24)],
                'format': 'netcdf',
            },
            data_path)
