#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the function to plot and save the map of minimum percentile rank across wind power density, seasonal variability and weather variability.

    This is Figure 2 of the paper.
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import settings as settings


# Set the font size.
plt.rc('font', size=20)


def plot_figure_2(lons_grid, lats_grid, minimum_percentile_rank):
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
    '''

    # Define the map levels.
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
    for ii in range(len(base_colors)-1):
        vals[int(N*ii/(len(base_colors)-1)):int(N*(ii+1)/(len(base_colors)-1)), 0] = np.linspace(base_colors[ii][0], base_colors[ii+1][0], int(N/(len(base_colors)-1)))
        vals[int(N*ii/(len(base_colors)-1)):int(N*(ii+1)/(len(base_colors)-1)), 1] = np.linspace(base_colors[ii][1], base_colors[ii+1][1], int(N/(len(base_colors)-1)))
        vals[int(N*ii/(len(base_colors)-1)):int(N*(ii+1)/(len(base_colors)-1)), 2] = np.linspace(base_colors[ii][2], base_colors[ii+1][2], int(N/(len(base_colors)-1)))
    custom_colormap = mpl.colors.ListedColormap(vals)

    # Define the map projection and the coordinate reference system of the data to plot.
    map_projection = ccrs.Robinson(central_longitude=0, globe=None)
    data_crs = ccrs.PlateCarree()

    # Initialize the figure and the axes of the subplots.
    fig, ax = plt.subplots(figsize=(15, 10), subplot_kw={'projection': map_projection})
    
    # Plot the coastlines.
    ax.coastlines()

    # Plot the minimum percentile rank.
    im = ax.contourf(lons_grid, lats_grid, minimum_percentile_rank, levels=map_levels, cmap=custom_colormap, transform=data_crs)

    # Plot the contour lines.
    ax.contour(lons_grid, lats_grid, minimum_percentile_rank, levels=map_levels, colors='k', zorder=6, linewidths=0.1, transform=data_crs)

    # Add a colorbar.
    cbar = plt.colorbar(im, shrink=0.5, location='bottom', pad=0.04)
    cbar.set_ticks(map_levels[::2])

    # Add the variable names.
    plt.annotate('Areas with abundant and reliable wind power', xy=(0.5, -0.20), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26)
    plt.annotate('Minimum percentile rank across power density,\nseasonal variability, and weather variability', xy=(0.5, -0.28), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=20)

    # Adjusts subplots to fit figure area.
    plt.tight_layout()

    # Set the title and save the figure.
    title = '/Figure 2 - Map of minimum percentile'
    fig.savefig(settings.figures_directory + title+'.png', bbox_inches = 'tight', dpi = 300)
    fig.savefig(settings.figures_directory + title+'.eps', format='eps', bbox_inches = 'tight', dpi = 300)