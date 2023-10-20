#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the settings used to process the climate data.
'''

import os


on_cluster = False

if on_cluster:
    working_directory = 'your_working_directory_on_your_hpc_cluster'
    climate_data_directory = 'your_climate_data_directory_on_your_hpc_cluster'
else:
    working_directory = os.getcwd()
    climate_data_directory = working_directory + '/ERA5_data'

general_data_directory = working_directory + '/data'
result_directory = working_directory + '/postprocessed_results'
figures_directory = working_directory + '/figures'

use_dask_mpi = False
use_dask = (True or use_dask_mpi)

if use_dask:
    chunks = {'latitude': 10, 'longitude': 10}
    chunks_climatological = {'latitude': 10, 'longitude': 10}
else:
    chunks = None
    chunks_climatological = None
