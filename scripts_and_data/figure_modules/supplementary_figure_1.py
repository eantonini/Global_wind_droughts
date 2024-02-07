#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the function to plot and save the maps of the land, coast and sea masks.

    This is Supplementary Figure 1 of the paper.
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import settings as settings


# Set the font size.
plt.rc('font', size=20)


def plot_supplementary_figure_1(lons_grid, lats_grid, geographic_masks):
    '''
    Plot and save the maps of the land, coast and sea masks.

    Parameters
    ----------
    lons_grid : numpy.ndarray
        Longitudes of the grid cells
    lats_grid : numpy.ndarray
        Latitudes of the grid cells
    geographic_masks : xarray.Dataset
        Dataset containing the land, coast, and sea masks
    '''

    # Define color maps, variable names, and panel letters.
    color_maps = ['Greens', 'Oranges', 'Blues']
    variable_name = ['Land regions', 'Coastal regions', 'Sea and ice sheets']
    panel_letter = ['a', 'b', 'c']
    
    # Define the map projection and the coordinate reference system of the data to plot.
    map_projection = ccrs.Robinson(central_longitude=0, globe=None)
    data_crs = ccrs.PlateCarree()

    # Initialize the figure and the axes of the subplots.
    fig, axs = plt.subplots(3, 1, figsize=(15, 20), subplot_kw={'projection': map_projection})

    # Plot the maps of the land, coast and sea masks.
    for ii, mask in enumerate(list(geographic_masks)):

        # Select the axis of the subplot.
        ax = axs[ii]

        # Plot the coastlines.
        ax.coastlines()

        # Plot the mask.
        ax.contourf(lons_grid, lats_grid, geographic_masks[mask].where(geographic_masks[mask], np.nan), levels=np.linspace(0, 1, 3), cmap=color_maps[ii], transform=data_crs)

        # Add the variable name and the panel letter.
        ax.annotate(variable_name[ii], xy=(1.08, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', color=mpl.colormaps.get_cmap(color_maps[ii])(0.8), weight='bold', fontsize=26, rotation=270) # type: ignore
        ax.annotate(panel_letter[ii], xy=(0.05, 0.95), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', fontsize=26)
    
    # Adjusts subplots to fit figure area.
    plt.tight_layout()

    # Set the title and save the figure.
    title = '/Supplementary Figure 1 - Geographic masks'
    fig.savefig(settings.figures_directory + title+'.png', bbox_inches = 'tight', dpi = 300)
    fig.savefig(settings.figures_directory + title+'.eps', format='eps', bbox_inches = 'tight', dpi = 300)
