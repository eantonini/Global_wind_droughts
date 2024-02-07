#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the function to plot and save the maps of the year of occurrence of the most severe wind drought, the intensity of the most severe wind drought, and the probability of a wind drought causing more than 400 hours of energy deficit.

    This is Supplementary Figure 10 of the paper.
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import settings as settings


# Set the font size.
plt.rc('font', size=20)


def plot_supplementary_figure_10(lons_grid, lats_grid, geographic_masks, wind_resource, energy_deficits, maximum_energy_deficits_for_wind_droughts):
    '''
    Plot and save the maps of the year of occurrence of the most severe wind drought, the intensity of the most severe wind drought, and the probability of a wind drought causing more than 400 hours of energy deficit.

    Parameters
    ----------
    lons_grid : numpy.ndarray
        Longitudes of the grid cells
    lats_grid : numpy.ndarray
        Latitudes of the grid cells
    geographic_masks : xarray.Dataset
        Dataset containing the masks of land, coast and sea
    wind_resource : xarray.Dataset
        Wind resource data
    energy_deficits : xarray.Dataset
        Energy deficits data
    maximum_energy_deficits_for_wind_droughts : xarray.DataArray
        Maximum energy deficits for wind droughts
    '''

    # Combine the land and coast masks.
    combined_geographic_masks = np.logical_or(geographic_masks['land'], geographic_masks['coast'])

    # Select the cells of interest considering regions with a mean power density greater than 150 W/m2 and the land and coast masks.
    cells_of_interest = np.logical_and(combined_geographic_masks, wind_resource['mean_wind_power_density'] > 150)
    
    # Select the maximum energy deficits for wind droughts, the year of occurrence, and the frequency of long deficits for the cells of interest.
    maximum_energy_deficits_of_interest = maximum_energy_deficits_for_wind_droughts.where(cells_of_interest, np.nan)
    year_of_max_deficit_of_interest = maximum_energy_deficits_for_wind_droughts['year'].where(cells_of_interest, np.nan)
    frequency_of_long_deficits = energy_deficits['wind_droughts'].where(energy_deficits['wind_droughts'] >= 400, np.nan).count(dim='year')/len(energy_deficits['year'])
    frequency_of_long_deficits_of_interest = frequency_of_long_deficits.where(cells_of_interest, np.nan)
    
    # Define a custom colormap.
    N = 256
    vals = np.ones((N, 4))
    base_colors = np.array([[255, 0, 0], [255, 255, 0], [0, 255, 0], [0, 0, 255], [255, 0, 255]])/255
    for ii in range(len(base_colors)-1):
        vals[int(N*ii/(len(base_colors)-1)):int(N*(ii+1)/(len(base_colors)-1)), 0] = np.linspace(base_colors[ii][0], base_colors[ii+1][0], int(N/(len(base_colors)-1)))
        vals[int(N*ii/(len(base_colors)-1)):int(N*(ii+1)/(len(base_colors)-1)), 1] = np.linspace(base_colors[ii][1], base_colors[ii+1][1], int(N/(len(base_colors)-1)))
        vals[int(N*ii/(len(base_colors)-1)):int(N*(ii+1)/(len(base_colors)-1)), 2] = np.linspace(base_colors[ii][2], base_colors[ii+1][2], int(N/(len(base_colors)-1)))
    custom_colormap = mpl.colors.ListedColormap(vals)

    # Define the map levels, colormaps, variable names, panel letters and ticks.
    map_levels = [energy_deficits['year'].values, np.linspace(0, 2000, 21), np.linspace(0,1,21)]
    map_colormaps = [custom_colormap, 'YlOrRd', 'GnBu']
    map_variable_name = ['Year of most\nsevere wind drought', 'Intensity of most\nsevere wind drought', 'Probability of a wind\ndrought causing more\nthan 400 hours\nof energy deficit']
    map_variable_secondary_name = ['[hours of energy deficit]']
    panel_letter = ['a', 'b', 'c']
    map_ticks = [np.linspace(1980,2020,5), map_levels[1][::4], np.linspace(0,1,6)]

    # Define the variables to plot.
    variables_to_plot = [year_of_max_deficit_of_interest, maximum_energy_deficits_of_interest, frequency_of_long_deficits_of_interest]

    # Define the map projection and the coordinate reference system of the data to plot.
    map_projection = ccrs.Robinson(central_longitude=0, globe=None)
    data_crs = ccrs.PlateCarree()
    
    # Initialize the figure and the axes of the subplots.
    fig, axs = plt.subplots(3, 1, figsize=(16, 20), subplot_kw={'projection': map_projection})

    # Plot the maps of the year of occurrence of the most severe wind drought, the intensity of the most severe wind drought, and the probability of a wind drought causing more than 400 hours of energy deficit.
    for ii in range(3):

        # Select the axis of the subplot.
        ax = axs[ii]

        # Plot the coastlines.
        ax.coastlines()

        # Plot the combined land and coast masks with a grey color.
        ax.contourf(lons_grid, lats_grid, combined_geographic_masks.where(combined_geographic_masks, np.nan), levels=np.linspace(0, 3, 10), cmap='Greys', transform=data_crs)

        # Plot the variable
        if ii == 1:
            im = ax.contourf(lons_grid, lats_grid, variables_to_plot[ii], levels=map_levels[ii], cmap=map_colormaps[ii], extend='max', transform=data_crs)
        else:
            im = ax.contourf(lons_grid, lats_grid, variables_to_plot[ii], levels=map_levels[ii], cmap=map_colormaps[ii], transform=data_crs)

        # Annotate the maps with the variable names.
        if ii == 0:
            ax.annotate(map_variable_name[ii], xy=(1.18, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26, rotation=270)
        elif ii == 1:
            ax.annotate(map_variable_name[ii], xy=(1.21, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26, rotation=270)
            ax.annotate(map_variable_secondary_name[ii-1], xy=(1.155, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=16, rotation=270)
        elif ii == 2:
             ax.annotate(map_variable_name[ii], xy=(1.215, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26, rotation=270)

        # Annotate the panel letter.
        ax.annotate(panel_letter[ii], xy=(0.05, 0.95), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26)

        # Annotate the description of the grey regions.
        if ii == 2:
            cmap = mpl.colormaps.get_cmap('Greys') # type: ignore
            ax.annotate('Grey regions have a mean power density\nless than 150 W m$^{-2}$', xy=(0.5, -0.10), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=20, color=cmap(1/3))
        
        # Add a colorbar. 
        cbar = plt.colorbar(im, shrink=0.6, ax=ax, pad=0.03)
        cbar.set_ticks(map_ticks[ii])
    
    # Adjust the space between the subplots.
    plt.tight_layout()

    # Set the title and save the figure.
    title = '/Supplementary Figure 10 - Maps of drought analysis'
    fig.savefig(settings.figures_directory + title+'.png', bbox_inches = 'tight', dpi = 300)
    fig.savefig(settings.figures_directory + title+'.eps', format='eps', bbox_inches = 'tight', dpi = 300)