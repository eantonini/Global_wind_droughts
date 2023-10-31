#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the function to plot and save the maps of percentile ranks of wind power density, seasonal variability and weather variability.

    This is Figure 1 of the paper.
'''

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import settings as settings


def plot_figure_1(lons_grid, lats_grid, percentile_rank_of_wind_power_density, percentile_rank_of_seasonal_variability, percentile_rank_of_weather_variability, wind_resource, energy_deficits, land_and_coast_mask):
    '''
    Plot and save the maps of percentile ranks of wind power density, seasonal variability and weather variability.
    
    Parameters
    ----------
    lons_grid : numpy.ndarray
        Longitudes of the grid cells
    lats_grid : numpy.ndarray
        Latitudes of the grid cells
    percentile_rank_of_wind_power_density : xarray.DataArray
        Percentile rank of the mean wind power density
    percentile_rank_of_seasonal_variability : xarray.DataArray
        Percentile rank of the energy deficits for the seasonal variability
    percentile_rank_of_weather_variability : xarray.DataArray
        Percentile rank of the mean energy deficits for the weather variability
    wind_resource : xarray.Dataset
        Wind resource data
    energy_deficits : xarray.Dataset
        Energy deficits data
    '''
    
    # Define the map levels, colormaps, variable names, and panel letters.
    map_levels = np.linspace(0, 100, 11)
    map_colormaps = ['plasma','viridis','cividis']
    map_variable_name = ['Mean power\ndensity\n[W m$^{-2}$]', 'Seasonal\nvariability\nmetric', 'Weather\nvariability\nmetric']
    panel_letter = ['a', 'b', 'c']

    # Define the variables to plot and the original variables.
    variables_to_plot = [percentile_rank_of_wind_power_density, percentile_rank_of_seasonal_variability, percentile_rank_of_weather_variability]
    original_variables = [wind_resource['mean_wind_power_density'], energy_deficits['normalized_seasonal_variability'], energy_deficits['normalized_weather_variability'].mean(dim='year')]
    
    # Define the map projection and the coordinate reference system of the data to plot.
    map_projection = ccrs.Robinson(central_longitude=0, globe=None)
    data_crs = ccrs.PlateCarree()

    # Initialize the figure and the axes of the subplots, and set the font size.
    fig, axs = plt.subplots(3, 1, figsize=(15, 20), subplot_kw={'projection': map_projection})
    plt.rc('font', size=20)

    # Plot the maps of percentile ranks of wind power density, seasonal variability and weather variability.
    for ii in range(3):
        
        # Select the axis of the subplot.
        ax = axs[ii]

        # Plot the coastlines.
        ax.coastlines()

        # Plot the percentile rank.
        im = ax.contourf(lons_grid, lats_grid, variables_to_plot[ii], levels=map_levels, cmap=map_colormaps[ii], transform=data_crs)

        # Plot the contour lines.
        ax.contour(lons_grid, lats_grid, variables_to_plot[ii], levels=map_levels, colors='k', linewidths=0.1, transform=data_crs)

        # Annotate the panel letter.
        ax.annotate(panel_letter[ii], xy=(0.05, 0.95), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26)
        
        # Add a colorbar and remove the ticks.
        cax = ax.inset_axes([1.08, 0.15, 0.02, 0.6])
        cbar = plt.colorbar(im, ax=ax, cax=cax)
        cbar.set_ticks([])

        # Define new tick values.
        tick_values = map_levels[::2]
        
        # For each tick value of the percentile rank, calculate the corresponding value of the original variable and annotate it.
        for kk in range(len(tick_values)):

            # Remove the values of the original variable that are not on land or coast.
            original_variables[ii] = original_variables[ii].where(land_and_coast_mask, np.nan)

            # Calculate the value of the original variable. The energy deficits have a descending percentile rank, so the tick values are reversed.
            if ii > 0:
                tick_percentile = original_variables[ii].quantile(tick_values[-1-kk]/100, skipna=True)
            else:
                tick_percentile = original_variables[ii].quantile(tick_values[kk]/100, skipna=True)
            
            # Set the starting y position and the y range of the ticks according to the dimension of the colorbar.
            starting_y_position = 0.150
            y_range = 0.5998

            # Annotate the tick value of the percentile rank and the corresponding value of the original variable.
            ax.annotate('{:d}'.format(int(tick_values[kk])), xy=(1.110, starting_y_position+y_range*kk/(len(tick_values)-1)), xycoords='axes fraction', horizontalalignment='left', verticalalignment='center', fontsize=20)
            if ii > 0:
                ax.annotate('{:.2f}'.format(tick_percentile), xy=(1.070, starting_y_position+y_range*kk/(len(tick_values)-1)), xycoords='axes fraction', horizontalalignment='right', verticalalignment='center', fontsize=20)
            else:
                ax.annotate('{:d}'.format(int(tick_percentile)), xy=(1.070, starting_y_position+y_range*kk/(len(tick_values)-1)), xycoords='axes fraction', horizontalalignment='right', verticalalignment='center', fontsize=20)
            
            # Annotate the ticks on both sides of the colorbar.
            ax.annotate('-', xy=(1.102, starting_y_position+y_range*kk/(len(tick_values)-1)), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=12)
            ax.annotate('-', xy=(1.078, starting_y_position+y_range*kk/(len(tick_values)-1)), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=12)
        
        # Annotate the variable names on both sides of the colorbar.
        ax.annotate(map_variable_name[ii], xy=(1.07, 0.90), xycoords='axes fraction', horizontalalignment='right', verticalalignment='center', weight='bold', fontsize=20)
        if ii == 0:
            ax.annotate('Ascending\npercentile\nrank', xy=(1.11, 0.90), xycoords='axes fraction', horizontalalignment='left', verticalalignment='center', weight='bold', fontsize=20)
        else:
            ax.annotate('Descending\npercentile\nrank', xy=(1.11, 0.90), xycoords='axes fraction', horizontalalignment='left', verticalalignment='center', weight='bold', fontsize=20)
    
    # Adjusts subplots to fit figure area.
    plt.tight_layout()

    # Set the title and save the figure.
    title = '/Figure 1 - Maps of individual percentiles'
    fig.savefig(settings.figures_directory + title+'.png', bbox_inches = 'tight', dpi = 300)
    fig.savefig(settings.figures_directory + title+'.eps', format='eps', bbox_inches = 'tight', dpi = 300)