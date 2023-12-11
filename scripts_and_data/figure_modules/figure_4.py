#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the function to plot and save the maps of trends in the wind power density, weather variability and wind droughts.

    This is Figure 4 of the paper.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs

import settings as settings


# Set the font size.
plt.rc('font', size=20)


def plot_figure_4(lons_grid, lats_grid, regression_of_wind_power_density, regression_of_weather_variability, regression_of_wind_droughts):
    '''
    Plot the maps of the trends in the wind power density, weather variability and wind droughts.
    
    Parameters
    ----------
    lons_grid : numpy.ndarray
        Longitudes of the grid
    lats_grid : numpy.ndarray
        Latitudes of the grid
    regression_of_wind_power_density : xarray.DataArray
        Linear regression of the wind power density
    regression_of_weather_variability : xarray.DataArray
        Linear regression of the weather variability
    regression_of_wind_droughts : xarray.DataArray
        Linear regression of the wind droughts
    '''

    # Define the map levels, colormaps, variable names, and panel letters.
    map_levels = [-10,-3,-1,-0.3,-0.1,-0.03,-0.01,0,0.01,0.03,0.1,0.3,1,3,10]
    map_colormaps = 'RdBu_r'
    map_variable_name = ['Annual\npercentage change in\npower density',
                         'Annual\npercentage change in\nweather variability',
                         'Annual\npercentage change in\nwind drought severity']
    panel_letter = ['a', 'b', 'c']
    
    # Define the variables to plot and their p-values.
    variables_to_plot = [regression_of_wind_power_density['slope'].copy()*100,
                         regression_of_weather_variability['slope'].copy()*100,
                         regression_of_wind_droughts['slope'].copy()*100]
    p_value_of_variables_to_plot = [regression_of_wind_power_density['p_value'].copy(),
                                    regression_of_weather_variability['p_value'].copy(),
                                    regression_of_wind_droughts['p_value'].copy()]
    
    # Define the map projection and the coordinate reference system of the data to plot.
    map_projection = ccrs.Robinson(central_longitude=0, globe=None)
    data_crs = ccrs.PlateCarree()

    # Initialize the figure and the axes of the subplots.
    fig, axs = plt.subplots(3, 1, figsize=(16, 20), subplot_kw={'projection': map_projection})

    # Set the hatch linewidth.
    plt.rcParams['hatch.linewidth'] = 0.5
    
    # Plot the maps of trends in the wind power density, weather variability and wind droughts.
    for ii in range(3):

        # Select the axis of the subplot.
        ax = axs[ii]

        # Plot the coastlines.
        ax.coastlines()

        # Plot the trends.
        variables_to_plot[ii] = variables_to_plot[ii].where(p_value_of_variables_to_plot[ii] <= 0.05, np.nan)
        im = ax.contourf(lons_grid, lats_grid, variables_to_plot[ii], levels=map_levels, cmap=map_colormaps, norm=colors.SymLogNorm(linthresh=0.01, linscale=0.01), transform=data_crs) # type: ignore

        # Plot the p-values.
        p_value_of_variables_to_plot[ii] = p_value_of_variables_to_plot[ii].where(p_value_of_variables_to_plot[ii] > 0.05, -1).where(p_value_of_variables_to_plot[ii] <= 0.05, 1)
        hatch = ax.contourf(lons_grid, lats_grid, p_value_of_variables_to_plot[ii], 2, colors='none', hatches=[None,'//'], extend='lower', transform=data_crs)

        # Set the color of the hatch lines.
        for collection in hatch.collections:
            collection.set_edgecolor('black')
            collection.set_linewidth(0)
        
        # Annotate the variable name and the panel letter.
        ax.annotate(map_variable_name[ii], xy=(1.21, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26, rotation=270)
        ax.annotate(panel_letter[ii], xy=(0.05, 0.95), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26)

        # Annotate the definition of the hatch lines.
        if ii == 2:
            ax.annotate('Hatches indicate regions\nwithout statistically significant trends', xy=(0.5, -0.10), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=20)

        # Add a colorbar.
        cbar = plt.colorbar(im, shrink=0.6, ax=ax, pad=0.03)
        cbar.set_ticks(map_levels[::2])
    
    # Adjusts subplots to fit figure area.
    plt.tight_layout()

    # Set the title and save the figure.
    title = '/Figure 4 - Maps of trends'
    fig.savefig(settings.figures_directory + title+'.png', bbox_inches = 'tight', dpi = 300)
    fig.savefig(settings.figures_directory + title+'.eps', format='eps', bbox_inches = 'tight', dpi = 300)