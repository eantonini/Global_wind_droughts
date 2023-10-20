#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains the function to plot and save the 2D histograms of the joint distribution of wind power density and energy deficits.

    This is Supplementary Figure 5 of the paper.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import settings as settings


def plot_histogram(fig, x_bin_edges, y_bin_edges, wind_power_density_ticks, median_storage_requirement, histogram,
                   main_histogram_colormap, main_histogram_norm, line_color, y_variable_name, y_variable_secondary_name, secondary_histogram_color,
                   secondary_histogram_x_levels, secondary_histogram_y_levels, y_axis_label, histogram_title, plot_colorbar, plot_median_label,
                   median_label_y, plot_secondary_histogram_x, plot_subplot_title, panel_letter, rounded_max_wind_power_density):
    '''
    Plot the 2D histogram of the joint distribution of wind power density and energy deficits.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure where to plot the histogram
    x_bin_edges : numpy.ndarray
        Bin edges of the wind power density
    y_bin_edges : numpy.ndarray
        Bin edges of the energy deficits
    wind_power_density_ticks : numpy.ndarray
        Ticks of the wind power density
    median_storage_requirement : numpy.ndarray
        Median energy deficit of the joint distribution of wind power density and energy deficits
    histogram : numpy.ndarray
        2D histogram of the joint distribution of wind power density and energy deficits
    main_histogram_colormap : matplotlib.colors.Colormap
        Colormap of the main histogram
    main_histogram_norm : matplotlib.colors.Normalize
        Normalization of the main histogram
    line_color : str
        Color of the line of the median energy deficit
    y_variable_name : str
        Name of the y variable
    y_variable_secondary_name : str
        Name of the y variable for the secondary histogram
    secondary_histogram_color : str
        Color of the secondary histogram
    secondary_histogram_x_levels : numpy.ndarray
        Levels of the secondary histogram for the x axis
    secondary_histogram_y_levels : numpy.ndarray
        Levels of the secondary histogram for the y axis
    y_axis_label : bool
        Whether to plot the y axis label
    histogram_title : str
        Title of the histogram
    plot_colorbar : bool
        Whether to plot the colorbar
    plot_median_label : bool
        Whether to plot the median label
    median_label_y : float
        y coordinate of the median label
    plot_secondary_histogram_x : bool
        Whether to plot the secondary histogram for the x axis
    plot_subplot_title : bool
        Whether to plot the subplot title
    panel_letter : str
        Panel letter
    rounded_max_wind_power_density : float
        Maximum wind power density
    '''

    # Define the bin edges of the 2D histogram.
    bin_edges_mesh_x, bin_edges_mesh_y = np.meshgrid(x_bin_edges, y_bin_edges)

    # Define the parameters of the histograms.
    main_histogram_left, main_histogram_width = 0.1, 0.55
    main_histogram_bottom, main_histogram_height = 0.27, 0.50
    spacing = 0.08
    secondary_histogram_height, secondary_histogram_width = 0.10, 0.15
    rect_2D_histogram = [main_histogram_left, main_histogram_bottom, main_histogram_width, main_histogram_height]
    rect_histogram_x = [main_histogram_left, main_histogram_bottom + main_histogram_height + spacing, main_histogram_width, secondary_histogram_height]
    rect_histogram_y = [main_histogram_left + main_histogram_width + spacing, main_histogram_bottom, secondary_histogram_width, main_histogram_height]
    x_bar_width = 1/25*rounded_max_wind_power_density
    y_bar_width = 1/25
    
    # Add a subplot to the figure.
    ax = fig.add_axes(rect_2D_histogram)

    # Plot the 2D histogram of the joint distribution of wind power density and energy deficits.
    im = ax.pcolormesh(bin_edges_mesh_x, bin_edges_mesh_y, histogram, cmap=main_histogram_colormap, norm=main_histogram_norm)

    # Set the ticks of the wind power density.
    ax.set_xticks(np.array([0,0.333,0.666,1])*rounded_max_wind_power_density)
    ax.set_xticklabels((np.round(np.insert(np.logspace(-1.7, 0, 3), 0, 0)*rounded_max_wind_power_density/10)*10).astype(int))

    # Plot the median energy deficit.
    ax.plot(wind_power_density_ticks, median_storage_requirement, '--', color=line_color)
    
    # Plot the y variable name.
    if y_axis_label:
        ax.annotate(y_variable_name, xy=(-0.61, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', rotation=90, weight='bold') # -0.4
        ax.annotate(y_variable_secondary_name, xy=(-0.36, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', rotation=90, fontsize=14) # -0.24

    # Plot the median label.
    if plot_median_label:
        ax.annotate('Median value', xy=(0.4, median_label_y), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', color=line_color)
    
    # Plot the subplot title.
    if plot_subplot_title:
        ax.annotate(histogram_title, xy=(0.65, 1.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold', color=secondary_histogram_color)
    
    # Plot the colorbar.
    if plot_colorbar:
        ax.set_xlabel('Mean power density', weight='bold')
        cax = fig.add_axes([0.1, 0.1, 0.75, 0.03])
        cbar = fig.colorbar(im, cax=cax, shrink=0.8, orientation='horizontal')
        cbar.ax.set_xticklabels(['0','','10$^{-5}$','','10$^{-4}$','','10$^{-3}$','','10$^{-2}$','','10$^{-1}$'])
        cbar.ax.set_xlabel('Probability', weight='bold')
    
    # Plot the panel letter.
    ax.annotate(panel_letter, xy=(0.9, 0.9), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', weight='bold')
    
    # Plot the secondary histograms on the x axis.
    if plot_secondary_histogram_x:
        ax_histogram_x = fig.add_axes(rect_histogram_x, sharex=ax)
        ax_histogram_x.bar(x_bin_edges[:-1], np.sum(histogram,axis=0), align='edge', width=x_bar_width, color=secondary_histogram_color, alpha=0.7)
        ax_histogram_x.tick_params(axis='x', labelbottom=False)
        ax_histogram_x.set_yticks(secondary_histogram_x_levels)
    
    # Plot the secondary histograms on the y axis.
    ax_histogram_y = fig.add_axes(rect_histogram_y, sharey=ax)
    ax_histogram_y.barh(y_bin_edges[:-1], np.sum(histogram, axis=1), align='edge', height=y_bar_width, color=secondary_histogram_color, alpha=0.7)
    ax_histogram_y.tick_params(axis='y', labelleft=False)
    ax_histogram_y.set_xticks(secondary_histogram_y_levels)
    ax_histogram_y.set_xticklabels(secondary_histogram_y_levels,rotation=270)
    ax_histogram_y.xaxis.tick_top()


def plot_supplementary_figure_5(wind_resource, histogram_of_seasonal_variability, histogram_of_weather_variability, median_of_seasonal_variability, median_of_weather_variability):
    '''
    Plot and save the 2D histograms of the joint distribution of wind power density and seasonal variability, and of wind power density and weather variability.

    Parameters
    ----------
    histogram_of_seasonal_variability : xarray.Dataset
        Dataset containing the 2D histograms of the joint distribution of wind power density and seasonal variability
    histogram_of_weather_variability : xarray.Dataset
        Dataset containing the 2D histograms of the joint distribution of wind power density and weather variability
    median_of_seasonal_variability : xarray.Dataset
        Dataset containing the median energy deficit of the joint distribution of wind power density and seasonal variability
    median_of_weather_variability : xarray.Dataset
        Dataset containing the median energy deficit of the joint distribution of wind power density and weather variability
    '''

    # Get the max wind power density.
    rounded_max_wind_power_density = np.ceil(wind_resource['mean_wind_power_density'].max().values/10)*10

    # Define the levels of the main and secondary histograms.
    main_histogram_levels = np.insert(np.logspace(-5.5,-1,10),0,0)
    secondary_histogram_x_levels = np.linspace(0,0.3,2)
    secondary_histogram_y_levels = [np.linspace(0,0.3,2), np.linspace(0,0.6,2)]
    
    # Define the parameters of the histograms and whether to plot the colorbar, the median label, the secondary histogram, the subplot title and the panel letter.
    plot_colorbar = [False, True]
    plot_median_label = [True, False]
    median_label_y =[[0.25, 0.25, 0.25], [0.25, 0.25, 0.25]]
    plot_secondary_histogram_x = [True, False]
    plot_subplot_title = [True, False]
    panel_letter = [['a', 'b', 'c'], ['d', 'e', 'f']]
    y_variable_name = ['Seasonal variability', 'Weather variability']
    y_variable_secondary_name = ['Energy deficit\nfor seasonal profile\nto generate constant power\n[fraction of maximum]',
                                 'Mean energy deficit\nfor individual year profile\nto generate seasonal power\n[fraction of maximum]']
    
    # Define the colormaps, the line colors, the titles and whether to plot the y axis label.
    main_histogram_colormaps = [plt.colormaps['Greens'],  plt.colormaps['Oranges'], plt.colormaps['Blues']] # type: ignore
    main_histogram_norms = [colors.BoundaryNorm(main_histogram_levels, ncolors=main_histogram_colormaps[0].N, clip=True),
                            colors.BoundaryNorm(main_histogram_levels, ncolors=main_histogram_colormaps[1].N, clip=True),
                            colors.BoundaryNorm(main_histogram_levels, ncolors=main_histogram_colormaps[2].N, clip=True)]
    line_colors = ['black', 'blue', 'red']
    secondary_histogram_colors = ['green', 'orange', 'blue']
    histogram_title = ['Grid cells over land', 'Grid cells over coast', 'Grid cells over sea\nand ice sheets']
    histogram_y_axis_label = [True, False, False]

    # Define the bin edges and the ticks of the wind power density.
    x_bin_edges = np.linspace(0,1,26)*rounded_max_wind_power_density
    y_bin_edges = np.linspace(0,1,26)
    wind_power_density_ticks = np.linspace(0.02,0.98,25)*rounded_max_wind_power_density
    
    # Initialize the first figure and the axes of the subplots, and set the font size.
    fig = plt.figure(figsize=(16, 7))
    plt.rc('font', size=20)

    # Add three subfigures to the figure.
    subfigsnest_histograms = fig.subfigures(1, 3)

    # Plot the histograms of the joint distribution of wind power density and seasonal variability.
    for ii, case in enumerate(list(histogram_of_seasonal_variability.data_vars)):
        plot_histogram(subfigsnest_histograms[ii], x_bin_edges, y_bin_edges, wind_power_density_ticks, median_of_seasonal_variability[case].values,
                       histogram_of_seasonal_variability[case].values, main_histogram_colormaps[ii], main_histogram_norms[ii], line_colors[ii], y_variable_name[0],
                       y_variable_secondary_name[0], secondary_histogram_colors[ii], secondary_histogram_x_levels, secondary_histogram_y_levels[0], histogram_y_axis_label[ii],
                       histogram_title[ii], plot_colorbar[0], plot_median_label[0], median_label_y[0][ii], plot_secondary_histogram_x[0], plot_subplot_title[0], panel_letter[0][ii], rounded_max_wind_power_density)
    
    # Set the title and save the first figure.
    title = '/Supplementary Figure 4 - Histogram of wind resource and seasonal variability'
    fig.savefig(settings.figures_directory + title+'.png', bbox_inches = 'tight', dpi = 300)

    # Initialize the second figure and the axes of the subplots, and set the font size.
    fig = plt.figure(figsize=(16, 7))
    plt.rc('font', size=20)

    # Add three subfigures to the figure.
    subfigsnest_histograms = fig.subfigures(1, 3)

    # Plot the histograms of the joint distribution of wind power density and weather variability.
    for ii, case in enumerate(list(histogram_of_weather_variability.data_vars)):
        plot_histogram(subfigsnest_histograms[ii], x_bin_edges, y_bin_edges, wind_power_density_ticks, median_of_weather_variability[case].values,
                       histogram_of_weather_variability[case].values, main_histogram_colormaps[ii], main_histogram_norms[ii], line_colors[ii], y_variable_name[1],
                       y_variable_secondary_name[1], secondary_histogram_colors[ii], secondary_histogram_x_levels, secondary_histogram_y_levels[1], histogram_y_axis_label[ii],
                       histogram_title[ii], plot_colorbar[1], plot_median_label[1], median_label_y[1][ii], plot_secondary_histogram_x[1], plot_subplot_title[1], panel_letter[1][ii], rounded_max_wind_power_density)

    # Set the title and save the second figure.
    title = '/Supplementary Figure 5 - Histogram of wind resource and weather variability'
    fig.savefig(settings.figures_directory + title+'.png', bbox_inches = 'tight', dpi = 300)
