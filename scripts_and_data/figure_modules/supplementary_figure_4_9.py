#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the function to plot and save the map of minimum percentile rank across wind power density, seasonal variability and weather variability for each continent.

    These are Supplementary Figures 4 - 9 of the paper.
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import settings as settings


# Set the font size.
plt.rc('font', size=20)


def plot_supplementary_figure_4_9(lons_grid, lats_grid, minimum_percentile_rank, continent):
    '''
    Plot and save the map of minimum percentile rank across wind power density, seasonal variability and weather variability.

    Parameters
    ----------
    lons_grid : numpy.ndarray
        Longitudes of the grid cells
    lats_grid : numpy.ndarray
        Latitudes of the grid cells
    minimum_percentile_rank : xarray.DataArray
        Minimum percentile rank across wind power density, seasonal variability and weather variability
    continent : str
        Continent of interest
    '''

    # Define the extent of the map.
    if continent == 'Europe':
        ii = 4
        map_extent = [-20, 35, 30, 75]
    elif continent == 'North America':
        ii = 5
        map_extent = [-135, -65, 10, 82]
    elif continent == 'South America':
        ii = 6
        map_extent = [-90, -30, -60, 20]
    elif continent == 'Africa':
        ii = 7
        map_extent = [-20, 60, -35, 40]
    elif continent == 'Asia':
        ii = 8
        map_extent = [65, 140, 3, 80]
    elif continent == 'Oceania':
        ii = 9
        map_extent = [110, 175, -53, 0]

    # Define the map levels and colormap.
    map_levels = np.linspace(0, 100, 11)

    # Get the original colormap.
    original_colormap = mpl.colormaps['turbo']

    # Extract the base colors of the original colormap.
    base_colors = original_colormap(np.linspace(0, 1, 11))
    
    # Modify some base colors to soften some green tones.
    base_colors[4,:] = np.array([163, 247, 131, 255])/255 # original [70.4, 247.6, 131.7, 255]
    base_colors[5,:] = np.array([221, 252, 60, 255])/255 # original [164.1, 252.4, 59.6, 255]
    base_colors[6,:] = np.array([255, 244, 126, 255])/255 # original [225.2, 220.7, 55.3, 255]

    # Create a custom colormap.
    N = (len(base_colors) - 1) * 20
    vals = np.ones((N, 4))    
    for jj in range(len(base_colors)-1):
        vals[int(N*jj/(len(base_colors)-1)):int(N*(jj+1)/(len(base_colors)-1)), 0] = np.linspace(base_colors[jj][0], base_colors[jj+1][0], int(N/(len(base_colors)-1)))
        vals[int(N*jj/(len(base_colors)-1)):int(N*(jj+1)/(len(base_colors)-1)), 1] = np.linspace(base_colors[jj][1], base_colors[jj+1][1], int(N/(len(base_colors)-1)))
        vals[int(N*jj/(len(base_colors)-1)):int(N*(jj+1)/(len(base_colors)-1)), 2] = np.linspace(base_colors[jj][2], base_colors[jj+1][2], int(N/(len(base_colors)-1)))
    custom_colormap = mpl.colors.ListedColormap(vals)

    # Define the map projection and the coordinate reference system of the data to plot.
    map_projection = ccrs.Stereographic(central_longitude=map_extent[0]+(map_extent[1]-map_extent[0])/2, central_latitude=map_extent[2]+(map_extent[3]-map_extent[2])/2)
    data_crs = ccrs.PlateCarree()

    # Initialize the figure and the axes of the subplots.
    fig, ax = plt.subplots(figsize=(21, 14), subplot_kw={'projection': map_projection})

    ax.set_extent(map_extent, crs=data_crs)
    
    # Plot the coastlines.
    ax.coastlines(resolution='50m')

    # Plot the minimum percentile rank.
    im = ax.contourf(lons_grid, lats_grid, minimum_percentile_rank, levels=map_levels, cmap=custom_colormap, transform=data_crs)

    # Plot the contour lines.
    ax.contour(lons_grid, lats_grid, minimum_percentile_rank, levels=map_levels, colors='k', zorder=6, linewidths=0.1, transform=data_crs)

    # Add a colorbar.
    cbar = plt.colorbar(im, shrink=0.5, location='bottom', pad=0.04)
    cbar.set_ticks(map_levels[::2])

    # Add the variable names.
    plt.annotate('Areas with abundant and reliable wind power', xy=(0.5, -0.17), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26)
    plt.annotate('Minimum percentile rank across power density,\nseasonal variability, and weather variability', xy=(0.5, -0.23), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=20)

    # Adjusts subplots to fit figure area.
    plt.tight_layout()

    # Set the title and save the figure.
    title = '/Supplementary Figure {:d} - Map of minimum percentile in '.format(ii) + continent
    fig.savefig(settings.figures_directory + title+'.png', bbox_inches = 'tight', dpi = 300)