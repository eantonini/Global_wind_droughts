#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the function to plot and save the map of the spatial extension of the most severe wind droughts in each continent and in each year.

    This is Figure 3 of the paper.
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs

import settings as settings


# Set the font size.
plt.rc('font', size=20)


def plot_figure_3(lons_grid, lats_grid, geographic_masks, continent_masks, wind_resource, drought_extension):
    '''
    Plot and save the map of the spatial extension of the most severe wind droughts in each continent and in each year.

    Parameters
    ----------
    lons_grid : numpy.ndarray
        Longitudes of the grid cells
    lats_grid : numpy.ndarray
        Latitudes of the grid cells
    geographic_masks : xarray.Dataset
        Dataset containing the land, coast, and sea masks
    continent_masks : xarray.Dataset
        Continent masks
    wind_resource : xarray.Dataset
        Wind resource data
    drought_extension : xarray.Dataset
        Spatial extension of the most severe wind droughts in each continent and in each year
    '''
    
    # Define the names and colors of each continent.
    continent_names = list(drought_extension.data_vars)
    continent_names = [continent.replace(' ', '\n') for continent in continent_names]
    continent_color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:cyan']
    
    # Define the start and end year for the analysis of the maximum drought extension.
    start_year = drought_extension['year'].values[0]
    end_year = drought_extension['year'].values[-1]

    # Calculate the maximum drought extension.
    max_drought_extension = drought_extension.max().to_array().max().values

    # Select the area with power density ≥ 150 W/m^2.
    wind_resource_of_interest = wind_resource['mean_wind_power_density'].where(wind_resource['mean_wind_power_density'] >= 150, np.nan)

    # Select the cells that are land or coast.
    wind_resource_of_interest = wind_resource_of_interest.where(np.logical_or(geographic_masks['land'], geographic_masks['coast']), np.nan)

    # Define the map projection and the coordinate reference system of the data to plot.
    map_projection = ccrs.Robinson(central_longitude=0, globe=None)
    data_crs = ccrs.PlateCarree()

    # Initialize the figure.
    fig, ax = plt.subplots(figsize=(12, 9))

    # Plot the spatial extension of the most severe wind droughts in each continent and in each year.
    for ii, continent in enumerate(list(drought_extension.data_vars)):

        for year in range(start_year, end_year+1):

            # Plot a rectangle proportional to the drought extension.
            ax.add_patch(Rectangle((ii+1-drought_extension[continent].sel(year=year).values / (2*max_drought_extension), year-0.5),
                                   drought_extension[continent].sel(year=year).values / max_drought_extension, 1, facecolor=continent_color[ii], fill=True))
    
    # Plot the x and y ticks and labels.
    ax.set_xticks(np.linspace(1,len(continent_names),len(continent_names)))
    ax.set_xticklabels(continent_names, verticalalignment='center',position=(0,-0.05))
    for ii in range(len(continent_color)):
        ax.get_xticklabels()[ii].set_color(continent_color[ii])
    ax.set_yticks(np.linspace(np.round(start_year/10)*10,np.round(end_year/10)*10,int((np.round(end_year/10)*10-np.round(start_year/10)*10)/10+1)))
    ax.set_ylabel('Year')

    # Set the limits of the x and y axes.
    ax.set_xlim(0.4,6.6)
    ax.set_ylim(start_year-9,end_year+17)

    # Annotate the title and the legend.
    ax.annotate('Distribution across years\nof the fraction of continents affected\nby the most severe wind droughts',
                xy=(0.04, 0.90), xycoords='axes fraction', horizontalalignment='left', verticalalignment='center', weight='bold', fontsize=20)
    ax.add_patch(Rectangle((0.57, 1973), 0.1/max_drought_extension, 1, facecolor='k', fill=True))
    ax.annotate('10% of continent area with power density ≥ 150 W m$^{-2}$',
                xy=(0.108, 0.049), xycoords='axes fraction', horizontalalignment='left', verticalalignment='center', weight='bold', fontsize=20)

    # In the top right corner, add a subplot of the world map with the regions with power density ≥ 150 W/m^2
    ax_subplot = fig.add_axes([0.67,0.745,0.25,0.25], projection=map_projection)

    # Plot the coastlines.
    ax_subplot.coastlines() # type: ignore

    # Loop over each continent.
    for ii, continent in enumerate(list(continent_masks.data_vars)):

        # Select the area with power density ≥ 150 W/m^2 in the continent.
        region_of_interest = wind_resource_of_interest.where(continent_masks[continent], np.nan)
        region_of_interest = region_of_interest.where(region_of_interest.isnull(), 1)

        # Plot the area with power density ≥ 150 W/m^2 in the continent.
        ax_subplot.contourf(lons_grid, lats_grid, region_of_interest, transform=data_crs, colors=continent_color[ii])

    # Adjusts subplots to fit figure area.
    plt.tight_layout()

    # Set the title and save the figure.
    title = '/Figure 3 - Drought analysis per continent'
    fig.savefig(settings.figures_directory + title+'.png', bbox_inches = 'tight', dpi = 300)
    fig.savefig(settings.figures_directory + title+'.eps', format='eps', bbox_inches = 'tight', dpi = 300)