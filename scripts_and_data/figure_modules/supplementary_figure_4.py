#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the function to plot and save the maps of trends in the wind power density, seasonal variability and weather variability.

    This is Supplementary Figure 4 of the paper.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs

import settings as settings


def plot_supplementary_figure_4(lons_grid, lats_grid, wind_resource, energy_deficits):
    '''
    Plot the maps of wind power density and variability.

    Parameters
    ----------
    lons_grid : numpy.ndarray
        Longitudes of the grid cells
    lats_grid : numpy.ndarray
        Latitudes of the grid cells
    wind_resource : xarray.Dataset
        Wind resource data
    energy_deficits : xarray.Dataset
        Energy deficits data
    '''

    # Define the map levels, colormaps, variable names, and panel letters.
    exponent_1 = -2.5
    exponent_2 = -1.15
    map_level_1 = np.insert(np.logspace(exponent_1, 0, 10),0,0)
    map_level_2 = np.round(np.insert(np.logspace(exponent_2, 0, 10),0,0),2)
    map_levels = [map_level_1, map_level_2, map_level_2]
    map_linthresh = [10**(exponent_1), 10**(exponent_2), 10**(exponent_2)]
    map_colormaps = ['plasma','viridis_r','cividis_r']
    map_variable_name = ['Mean power density', 'Seasonal variability', 'Weather variability']
    map_variable_secondary_name = ['[W m$^{-2}$]', 'Energy deficit for seasonal profile\nto generate constant power\n[fraction of year]', 'Mean energy deficit for individual year\nprofile to generate seasonal power\n[fraction of year]']
    panel_letter = ['a', 'b', 'c']
    
    # Define the variables to plot.
    variables_to_plot = [wind_resource['mean_wind_power_density'], energy_deficits['normalized_seasonal_variability'], energy_deficits['normalized_weather_variability'].mean(dim='year')]
    
    # Define the normalization factor used only for the wind power density.
    rounded_max_wind_power_density = np.ceil(wind_resource['mean_wind_power_density'].max().values/10)*10
    normalization_factor = [rounded_max_wind_power_density, 1, 1]

    # Define the map projection and the coordinate reference system of the data to plot.
    map_projection = ccrs.Robinson(central_longitude=0, globe=None)
    data_crs = ccrs.PlateCarree()
    
    # Initialize the figure and the axes of the subplots, and set the font size.
    fig, axs = plt.subplots(3, 1, figsize=(15, 20), subplot_kw={'projection': map_projection})
    plt.rc('font', size=20)

    # Plot the maps of wind power density and variability.
    for ii in range(3):

        # Select the axis of the subplot.
        ax = axs[ii]

        # Plot the coastlines.
        ax.coastlines()

        # Plot the variable.
        im = ax.contourf(lons_grid, lats_grid, variables_to_plot[ii]/normalization_factor[ii], levels=map_levels[ii], cmap=map_colormaps[ii], norm=colors.SymLogNorm(linthresh=map_linthresh[ii], linscale=map_linthresh[ii]), transform=data_crs) # type: ignore

        # Plot the contour lines.
        ax.contour(lons_grid, lats_grid, variables_to_plot[ii], levels=map_levels[ii], colors='k', linewidths=0.1, transform=data_crs)

        # Add the variable name and the panel letter.
        if ii == 0:
            ax.annotate(map_variable_name[ii], xy=(1.205, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26, rotation=270) # 1.16
            ax.annotate(map_variable_secondary_name[ii], xy=(1.16, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=16, rotation=270)
        else:
            ax.annotate(map_variable_name[ii], xy=(1.25, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26, rotation=270)
            ax.annotate(map_variable_secondary_name[ii], xy=(1.18, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=16, rotation=270)
        ax.annotate(panel_letter[ii], xy=(0.05, 0.95), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26)

        # Add a colorbar.
        cbar = plt.colorbar(im, shrink=0.6, ax=ax, pad=0.03)

        # Set the ticks of the colorbar.
        cbar.set_ticks(map_levels[ii][::2])
        if ii == 0:
            cbar.ax.set_yticklabels((map_levels[ii][::2]*normalization_factor[ii]).astype(int))
    
    # Adjusts subplots to fit figure area.
    plt.tight_layout()

    # Set the title and save the figure.
    title = '/Supplementary Figure 4 - Maps of wind power density and variability'
    fig.savefig(settings.figures_directory + title+'.png', bbox_inches = 'tight', dpi = 300)
    fig.savefig(settings.figures_directory + title+'.eps', format='eps', bbox_inches = 'tight', dpi = 300)
