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

    # Define the map levels and colormap.
    map_levels = np.linspace(0, 100, 11)
    map_colormaps = 'turbo'

    # Define the map projection and the coordinate reference system of the data to plot.
    map_projection = ccrs.Robinson(central_longitude=0, globe=None)
    data_crs = ccrs.PlateCarree()

    # Initialize the figure and the axes of the subplots.
    fig, ax = plt.subplots(figsize=(15, 10), subplot_kw={'projection': map_projection})
    
    # Plot the coastlines.
    ax.coastlines()

    # Plot the minimum percentile rank.
    im = ax.contourf(lons_grid, lats_grid, minimum_percentile_rank, levels=map_levels, cmap=map_colormaps, transform=data_crs)

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