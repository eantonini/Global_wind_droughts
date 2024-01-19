#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script analyzes the results of the wind resource and energy deficits and produces the figures of the paper.
'''

import numpy as np
import os

import settings as settings
import utilities_for_data_analysis as utilities
from figure_modules import *


# Create the directory where the figures will be stored if it does not exist.
if not os.path.exists(settings.figures_directory):
    os.mkdir(settings.figures_directory)

# Define the longitude and latitude grid.
lons = np.array(np.arange(0,360.25,0.25))
lats = np.array(np.arange(90,-90.25,-0.25))
lons_grid, lats_grid = np.meshgrid(lons,lats)

# Calculate the grid cell areas.
cell_areas = utilities.get_grid_cell_area(lons, lats)

# Define the land, coast, sea and continent masks.
geographic_masks = utilities.get_land_coast_sea_masks()
continent_masks = utilities.get_continent_masks(lons, lats)
land_and_coast_mask = np.logical_or(geographic_masks['land'], geographic_masks['coast'])

# Read the wind resource and energy deficits.
wind_resource = utilities.read_wind_speed_and_power_density()
energy_deficits = utilities.read_energy_deficits()

# Calculate the 2D histograms of the joint distribution of wind power density and energy deficits for land, coast and sea.
histogram_of_seasonal_variability = utilities.get_histogram(geographic_masks, wind_resource['annual_mean_wind_power_density'].max(), wind_resource['mean_wind_power_density'], energy_deficits['normalized_seasonal_variability'])
histogram_of_weather_variability = utilities.get_histogram(geographic_masks, wind_resource['annual_mean_wind_power_density'].max(), wind_resource['annual_mean_wind_power_density'], energy_deficits['normalized_weather_variability'])

# Calculate the median energy deficit of the joint distribution of wind power density and energy deficits for land, coast and sea.
median_of_seasonal_variability = utilities.get_meadian_energy_deficit(histogram_of_seasonal_variability)
median_of_weather_variability = utilities.get_meadian_energy_deficit(histogram_of_weather_variability)

# Calculate a linear regression of the wind power density, weather variability and wind droughts, normalized by their mean values to get the relative change.
regression_of_wind_power_density = utilities.linear_regression(wind_resource['annual_mean_wind_power_density']/wind_resource['annual_mean_wind_power_density'].mean(dim='year'))
regression_of_weather_variability = utilities.linear_regression(energy_deficits['weather_variability']/energy_deficits['weather_variability'].mean(dim='year'))
regression_of_wind_droughts = utilities.linear_regression(energy_deficits['wind_droughts']/energy_deficits['wind_droughts'].mean(dim='year'))

# Calculate the percentile rank of the wind power density, seasonal variability and weather variability.
percentile_rank_of_wind_power_density = utilities.get_percentile_rank(wind_resource['mean_wind_power_density'], land_and_coast_mask)
percentile_rank_of_seasonal_variability = utilities.get_percentile_rank(-energy_deficits['seasonal_variability'], land_and_coast_mask)
percentile_rank_of_weather_variability = utilities.get_percentile_rank(-energy_deficits['weather_variability'].mean(dim='year'), land_and_coast_mask)

# Calculate the minimum percentile rank across the wind power density, seasonal variability and weather variability.
minimum_percentile_rank = np.minimum(np.minimum(percentile_rank_of_wind_power_density, percentile_rank_of_seasonal_variability), percentile_rank_of_weather_variability)

# Calculate the maximum energy deficit for wind droughts and the year of occurrence.
maximum_energy_deficits_for_wind_droughts = energy_deficits['wind_droughts'].isel(year=energy_deficits['wind_droughts'].argmax(dim='year'))

# Calculate the spatial extension of the most severe wind droughts in each continent and in each year.
drought_extension = utilities.get_drought_extension(cell_areas, geographic_masks, continent_masks, wind_resource, maximum_energy_deficits_for_wind_droughts)

# Calculate the percentage of land and coast areas with a mean wind power density equal to or greater than 150 W/m2.
land_and_coast_surface = cell_areas*land_and_coast_mask
land_and_coast_surface_with_high_wind_power = (wind_resource['mean_wind_power_density']>=150)*land_and_coast_surface
percentage_of_land_and_coast_surface_with_with_high_wind_power = land_and_coast_surface_with_high_wind_power.sum()/land_and_coast_surface.sum()*100

# Calculate the percentage of land and coast areas with a mean wind power density equal to or greater than 150 W/m2 and with a significant change in the wind power density, weather variability and wind droughts.
land_and_coast_surface_with_high_wind_power_and_significant_wind_power_change = (regression_of_wind_power_density['p_value'] <= 0.05)*land_and_coast_surface_with_high_wind_power
land_and_coast_surface_with_high_wind_power_and_significant_weather_variability_change = (regression_of_weather_variability['p_value'] <= 0.05)*land_and_coast_surface_with_high_wind_power
land_and_coast_surface_with_high_wind_power_and_significant_wind_drought_change = (regression_of_wind_droughts['p_value'] <= 0.05)*land_and_coast_surface_with_high_wind_power
percentage_of_land_and_coast_surface_with_with_high_wind_power_and_significant_wind_power_change = land_and_coast_surface_with_high_wind_power_and_significant_wind_power_change.sum()/land_and_coast_surface.sum()*100
percentage_of_land_and_coast_surface_with_with_high_wind_power_and_significant_weather_variability_change = land_and_coast_surface_with_high_wind_power_and_significant_weather_variability_change.sum()/land_and_coast_surface.sum()*100
percentage_of_land_and_coast_surface_with_with_high_wind_power_and_significant_wind_drought_change = land_and_coast_surface_with_high_wind_power_and_significant_wind_drought_change.sum()/land_and_coast_surface.sum()*100


# FIGURE 1
figure_1.plot_figure_1(lons_grid, lats_grid, percentile_rank_of_wind_power_density, percentile_rank_of_seasonal_variability, percentile_rank_of_weather_variability, wind_resource, energy_deficits, land_and_coast_mask)

# FIGURE 2
figure_2.plot_figure_2(lons_grid, lats_grid, minimum_percentile_rank)

# FIGURE 3
figure_3.plot_figure_3(lons_grid, lats_grid, geographic_masks, continent_masks, wind_resource, drought_extension)

# FIGURE 4
figure_4.plot_figure_4(lons_grid, lats_grid, regression_of_wind_power_density, regression_of_weather_variability, regression_of_wind_droughts)


# SUPPLEMENTARY FIGURE 1
supplementary_figure_1.plot_supplementary_figure_1(lons_grid, lats_grid, geographic_masks)

# SUPPLEMENTARY FIGURE 2
supplementary_figure_2.plot_supplementary_figure_2(lons_grid, lats_grid, wind_resource, energy_deficits)

# SUPPLEMENTARY FIGURE 3
supplementary_figure_3.plot_supplementary_figure_3(wind_resource, histogram_of_seasonal_variability, histogram_of_weather_variability, median_of_seasonal_variability, median_of_weather_variability)

# SUPPLEMENTARY FIGURES 4 TO 9
continents = ['Europe', 'North America', 'South America', 'Africa', 'Asia', 'Oceania']
for continent in continents:
    supplementary_figure_4_9.plot_supplementary_figure_4_9(lons_grid, lats_grid, minimum_percentile_rank, continent)

# SUPPLEMENTARY FIGURE 10
supplementary_figure_10.plot_supplementary_figure_10(lons_grid, lats_grid, geographic_masks, wind_resource, energy_deficits, maximum_energy_deficits_for_wind_droughts)

# SUPPLEMENTARY FIGURE 11
supplementary_figure_11.plot_supplementary_figure_11()