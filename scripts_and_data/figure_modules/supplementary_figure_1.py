#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the function to plot and save sample time series of wind power density, generation and demand.

    This is Supplementary Figure 1 of the paper.
'''

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

import settings as settings


def plot_supplementary_figure_1():
    '''
    Plot the sample time series of wind power density, generation and demand.
    '''

    # Read the wind power density time series.
    wind_speed_time_series = xr.open_dataset(settings.general_data_directory + '/wind_speed_time_series_at_lat_53.00_lon_3.00.nc', engine='netcdf4')

    # Remove every 29th of February from the time series.
    feb_29 = np.logical_and(wind_speed_time_series['time'].dt.month == 2, wind_speed_time_series['time'].dt.day == 29)
    wind_speed_time_series = wind_speed_time_series['100m wind speed'].loc[~feb_29.values] # type: ignore

    # Calculate the wind power density time series.
    wind_power_density_time_series = (0.5*wind_speed_time_series**3).rename('100m wind power')

    # Calculate the annual mean wind power density.
    annual_mean_wind_power_density = wind_power_density_time_series.resample(time='AS').mean('time')
    annual_mean_wind_power_density['time'] = annual_mean_wind_power_density['time'].dt.year
    annual_mean_wind_power_density = annual_mean_wind_power_density.rename({'time': 'year'})

    # Calculate the climatological wind power density time series.
    hour_of_year = xr.DataArray(data=np.tile(np.arange(365*24), 43), dims=['time'],
                                coords=dict(time=wind_power_density_time_series['time']), name='hour of year')
    climatological_wind_power_density_time_series = wind_power_density_time_series.groupby(hour_of_year).mean()

    # Select a representative year.
    selected_year = pd.date_range('2010-01-01', '2011-01-01', freq='H')[:-1]

    # Extract the wind power density time series for the selected year.
    selected_wind_power_density_time_series = wind_power_density_time_series.sel(time=selected_year)

    # Calculate the mean wind power density of the selected year and the climatological mean wind power density.
    mean_selected_wind_power_density_time_series = selected_wind_power_density_time_series.mean()
    climatological_mean_wind_power_density = climatological_wind_power_density_time_series.mean()
    
    # Define the variables to plot.
    variables_to_plot = [selected_wind_power_density_time_series, climatological_wind_power_density_time_series,
                         selected_wind_power_density_time_series/mean_selected_wind_power_density_time_series,
                         climatological_wind_power_density_time_series/climatological_mean_wind_power_density,
                         np.array([1,1]), climatological_wind_power_density_time_series/climatological_mean_wind_power_density]
    
    # Define the scale factors, y labels, descriptions, colors, plot titles, titles and panel letters of each figure.
    scale_factors = [10000, 2000, 15, 3, 3, 3]
    y_labels = ['Power density', 'Power density', 'Generation', 'Generation', 'Target generation', 'Target generation']
    descriptions = ['Mean power density', 'Mean power density', 'Mean unit generation', 'Mean unit generation', 'Constant unit generation', 'Mean unit generation']
    colors = ['green', 'green', 'orange', 'orange', 'blue', 'blue']
    plot_titles = ['Individual year', 'Climatological year']
    panel_letter = ['b', 'a', 'd', 'c', 'e', 'f']
    titles = ['/Supplementary Figure 1 - Source from individual year capacity factor time series',
              '/Supplementary Figure 1 - Source from climatological capacity factor time series',
              '/Supplementary Figure 1 - Generation time series from individual year capacity factor',
              '/Supplementary Figure 1 - Generation time series from climatological capacity factor',
              '/Supplementary Figure 1 - Demand time series of constant value',
              '/Supplementary Figure 1 - Demand time series from climatological capacity factor']

    # Plot the sample time series of wind power density, generation and demand.
    for ii in range(len(variables_to_plot)):

        # Initialize the figure and set the font size.
        plt.figure(figsize=(8,4))
        plt.rc('font', size=20)

        # Plot the time series.
        plt.plot(np.linspace(0,1,len(variables_to_plot[ii])),variables_to_plot[ii], color=colors[ii], linestyle='-', alpha=0.6)

        # Plot the mean value.
        plt.plot(np.array([0,1]), np.array([np.mean(variables_to_plot[ii]),np.mean(variables_to_plot[ii])]), color='k', linestyle='--')

        # Plot the x and y labels.
        plt.xlabel('Fraction of year')
        plt.ylabel(y_labels[ii])

        # Set the limits of the x and y axes.
        plt.xlim([-0.05,1.05])
        plt.ylim([0-scale_factors[ii]*0.05,scale_factors[ii]*1.05])

        # Annotate the description and the panel letter.
        plt.text(0.5, np.mean(variables_to_plot[ii])+scale_factors[ii]*0.1, descriptions[ii], horizontalalignment='center', weight='bold')
        plt.annotate(panel_letter[ii], xy=(0.95, 0.90), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26)

        # Add the title.
        if ii < 2:
            plt.title(plot_titles[ii])

        # Save the figure.
        plt.savefig(settings.figures_directory + titles[ii]+'.png', bbox_inches = 'tight', dpi = 300)
