#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Author: Enrico Antonini

Created on: 2023-11-01

License: GNU General Public License v3.0

Description:

    This script contains some utility functions used to analyze the results of the wind resource and energy deficits.
'''

import numpy as np
import xarray as xr
import geopandas as gp
import regionmask
from scipy import stats
from xhistogram.xarray import histogram

import settings as settings


def get_grid_cell_area(lons, lats):
    '''
    Calculate the area of each cell defined by the lat/lon grid.
    
    Source: https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2004/msg00023.html
    Source: https://en.wikipedia.org/wiki/Spherical_sector

    Parameters
    ----------
    lons : numpy.ndarray
        Array containing the longitude values of the grid cell midpoints
    lats : numpy.ndarray
        Array containing the latitude values of the grid cell midpoints

    Returns
    -------
    cell_areas : xarray.DataArray
        DataArray containing the area of each cell defined by the lat/lon grid
    '''

    # Define the latitude values of the grid cell boundaries.
    bounds_lat = np.insert(0.5*(lats[1:] + lats[:-1]), 0, lats[0])
    bounds_lat = np.insert(bounds_lat, len(bounds_lat), lats[-1])
    bounds_lat = xr.Dataset(data_vars={'upper_lat': (['latitude', 'longitude'], np.tile(bounds_lat[1:].T, (len(lons), 1)).T),
                                       'lower_lat': (['latitude', 'longitude'], np.tile(bounds_lat[:-1], (len(lons), 1)).T)},
                            coords={'longitude': lons, 'latitude': lats})

    # Define the grid resolution.
    delta_lon = 0.25

    # Define the Earth's radius.
    R_earth = 6.371*10**6

    # Calculate the area of each cell.
    cell_areas = 2*np.pi*R_earth**2*np.absolute(np.sin(bounds_lat['upper_lat']*np.pi/180) - np.sin(bounds_lat['lower_lat']*np.pi/180))*np.absolute(delta_lon)/360
    
    return cell_areas


def read_wind_speed_and_power_density():
    '''
    Read the wind speed and power density data from the results directory.

    Returns
    -------
    wind_resource : xarray.Dataset
        Dataset containing the wind speed and power density data
    '''

    # Initialize the wind resource dataset.
    wind_resource = xr.Dataset()

    # Read the wind speed and power density data.
    wind_resource['mean_wind_speed'] = xr.open_dataarray(settings.results_directory + '/mean_wind_speed.nc')
    wind_resource['annual_mean_wind_speed'] = xr.open_dataarray(settings.results_directory + '/annual_mean_wind_speed.nc')
    wind_resource['mean_wind_power_density'] = xr.open_dataarray(settings.results_directory + '/mean_wind_power_density.nc')
    wind_resource['annual_mean_wind_power_density'] = xr.open_dataarray(settings.results_directory + '/annual_mean_wind_power_density.nc')

    # Create additional element to be placed at the end of the original dataset along the longitude dimension.
    extend_longitude = wind_resource.sel(longitude=wind_resource['longitude'][0].values)
    extend_longitude['longitude'] = [wind_resource['longitude'][-1].values + 0.25]
                
    # Add the additional element to the original dataset. 
    wind_resource = xr.combine_by_coords([wind_resource, extend_longitude])

    return wind_resource


def read_energy_deficits():
    '''
    Read the energy deficit data from the results directory.

    Returns
    -------
    energy_deficits : xarray.Dataset
        Dataset containing the energy deficit data
    '''

    # Initialize the energy deficit dataset.
    energy_deficits = xr.Dataset()

    # Read the energy deficit data.
    energy_deficits['seasonal_variability'] = xr.open_dataarray(settings.results_directory + '/energy_deficit_for_seasonal_variability.nc')
    energy_deficits['weather_variability'] = xr.open_dataarray(settings.results_directory + '/energy_deficit_for_weather_variability.nc')
    energy_deficits['wind_droughts'] = xr.open_dataarray(settings.results_directory + '/energy_deficit_for_wind_droughts.nc')

    # Calculate the number of hours in each year to normalize the energy deficit data.
    start_year = energy_deficits['year'][0].values
    end_year = energy_deficits['year'][-1].values
    hours_in_year = []
    for year in range(start_year,end_year+1):
        if year%4 == 0:
            hours_in_year.append(366*24)
        else:
            hours_in_year.append(365*24) 
    hours_in_year = xr.DataArray(data=hours_in_year, dims=['year'], coords={'year': energy_deficits['year']}, name='hours_in_year')

    # Add normalized versions of the energy deficit data to the dataset.
    energy_deficits['normalized_seasonal_variability'] = energy_deficits['seasonal_variability']/(366*24)
    energy_deficits['normalized_weather_variability'] = energy_deficits['weather_variability']/hours_in_year
    energy_deficits['normalized_wind_droughts'] = energy_deficits['wind_droughts']/hours_in_year

    # Create additional element to be placed at the end of the original dataset along the longitude dimension.
    extend_longitude = energy_deficits.sel(longitude=energy_deficits['longitude'][0].values)
    extend_longitude['longitude'] = [energy_deficits['longitude'][-1].values + 0.25]
                
    # Add the additional element to the original dataset. 
    energy_deficits = xr.combine_by_coords([energy_deficits, extend_longitude])

    return energy_deficits


def add_buffer_layers(array, n_layers):
    '''
    Add buffer layers (extra longitude and latitude layers) to the array to avoid edge effects in the averaging step.

    Parameters
    ----------
    array : xarray.DataArray
        Array to which buffer layers should be added
    n_layers : int
        Number of buffer layers to add

    Returns
    -------
    array : xarray.DataArray
        Array with buffer layers added
    '''
        
    # Define the grid resolution.
    grid_resolution = 0.25

    # Add buffer layers to the array.
    for ii in range(1, n_layers+1):

        # Add longitude values to the right of the array.
        original_longitude, new_longitude = 0.00+grid_resolution*ii, 360.00+grid_resolution*ii
        column_to_add = xr.DataArray(data=np.atleast_2d(array.sel(longitude=original_longitude).values).T,
                                     coords=[('latitude', array['latitude'].values), ('longitude', [new_longitude])])
        array = array.combine_first(column_to_add)

        # Add longitude values to the left of the array.
        original_longitude, new_longitude = 360.00-grid_resolution*ii, 0.00-grid_resolution*ii
        column_to_add = xr.DataArray(data=np.atleast_2d(array.sel(longitude=original_longitude).values).T,
                                     coords=[('latitude', array['latitude'].values), ('longitude', [new_longitude])])
        array = array.combine_first(column_to_add)
        
        # Add latitude values to the top of the array.
        original_latitude, new_latitude = 90.00-grid_resolution*ii, 90.00+grid_resolution*ii
        row_to_add = xr.DataArray(data=np.atleast_2d(array.sel(latitude=original_latitude).values),
                                  coords=[('latitude', [new_latitude]), ('longitude', array['longitude'].values)])
        array = array.combine_first(row_to_add)

        # Add latitude values to the bottom of the array.
        original_latitude, new_latitude = -90.00+grid_resolution*ii, -90.00-grid_resolution*ii
        row_to_add = xr.DataArray(data=np.atleast_2d(array.sel(latitude=original_latitude).values), 
                                  coords=[('latitude', [new_latitude]), ('longitude', array['longitude'].values)])
        array = array.combine_first(row_to_add)
    
    return array


def get_land_coast_sea_masks():
    '''
    Calculate the land, coast, and sea masks (containing True or False values). Greenlad and Antarctica are included in the sea mask.

    Returns
    -------
    geographic_masks : xarray.Dataset
        Dataset containing the land, coast, and sea masks
    '''
    
    # Read the land-sea index.
    land_sea_index = xr.open_dataarray(settings.general_data_directory + '/land_sea_mask.nc')

    # Define the latitude and longitude values of the wanted grid.
    target_longitude = np.insert(land_sea_index['longitude'].values, len(land_sea_index['longitude'].values), land_sea_index['longitude'][-1].values + 0.25)
    target_latitude = land_sea_index['latitude'].values
    
    # Add columns (longitude values) to the land-sea index to match the original grid.
    column_to_add = xr.DataArray(data=np.atleast_2d(land_sea_index.sel(longitude=target_longitude[0]).values).T,
                                 coords=[('latitude', target_latitude), ('longitude', [target_longitude[-1]])])
    land_sea_index = land_sea_index.combine_first(column_to_add)

    # Calculate the land mask
    land_mask = land_sea_index > 0.7

    # Add columns (longitude values) and rows (latitude values) to the land mask to avoid edge effects in the following averaging step.
    n_layers = 1
    land_sea_index = add_buffer_layers(land_sea_index, n_layers)

    # Calculate a smoothed land-sea index by averaging over a 3x3 grid and then extract only the values for the original grid.
    # This is to exclude more cells around land when calculating the sea mask.
    land_sea_index_smooth = land_sea_index.rolling(longitude=1+2*n_layers, latitude=1+2*n_layers, center=True).mean()
    land_sea_index_smooth = land_sea_index_smooth.sel(latitude=target_latitude, longitude=target_longitude)

    # Calculate the sea mask.
    sea_mask = land_sea_index_smooth == 0

    # Calculate the coast mask.
    coast_mask = np.logical_and(land_sea_index_smooth>0, land_sea_index<=0.7)

    # Get a mask for all contries in the world, where each country is assigned a unique number.
    countries = regionmask.defined_regions.natural_earth_v5_0_0.countries_50.mask(target_longitude[:-1],target_latitude) # type: ignore
    countries = countries.rename({'lon': 'longitude', 'lat': 'latitude'})

    # Add columns (longitude values) to the countries mask to match the original grid.
    column_to_add = xr.DataArray(data=np.atleast_2d(countries.sel(longitude=target_longitude[0]).values).T,
                                 coords=[('latitude', countries['latitude'].values), ('longitude', [target_longitude[-1]])])
    countries = countries.combine_first(column_to_add)

    # Get a mask for Greenland, where Greenland is assigned the number 181.
    greenland_mask = countries==181

    # Add columns (longitude values) and rows (latitude values) to the Greenland mask to avoid edge effects in the following averaging step.
    n_layers = 8
    greenland_mask = add_buffer_layers(greenland_mask, n_layers)

    # Calculate a smoothed Greenland mask by averaging over a 17x17 grid and then extract only the values for the original grid.
    # This is to exclude more cells around Greenland when calculating the resulting mask.
    greenland_mask_smooth = greenland_mask.rolling(longitude=1+2*n_layers, latitude=1+2*n_layers, center=True).mean()
    greenland_mask_smooth = greenland_mask_smooth.sel(latitude=target_latitude, longitude=target_longitude)

    # Calculate the Greenland mask
    greenland_mask = greenland_mask_smooth>0

    # Get a mask for Antarctica, where all values south of 55 degrees south are excluded.
    # Alternatively, Antarctica could be defined as all cells with the country number 239.
    antartica_mask = xr.DataArray(data=np.full((len(target_latitude),len(target_longitude)), False),
                                  coords=[('latitude', target_latitude), ('longitude', target_longitude)])
    antartica_mask = antartica_mask.where(antartica_mask['latitude']>=-55, True)

    # Combine the Greenland and Antarctica masks to get a mask for all regions to exclude.
    regions_to_exclude = np.logical_and(np.invert(antartica_mask), np.invert(greenland_mask))

    # Combine the land, coast, and sea masks with the regions to exclude to get the final masks.
    land_mask = np.logical_and(land_mask, regions_to_exclude)
    coast_mask = np.logical_and(coast_mask, regions_to_exclude)
    sea_mask = np.logical_or(sea_mask, ~regions_to_exclude)

    # Create a dataset containing the resulting land, coast, and sea masks.
    geographic_masks = xr.Dataset({'land': land_mask, 'coast': coast_mask, 'sea': sea_mask})

    return geographic_masks


def get_continent_masks(lons, lats):
    '''
    Calculate the continent masks (containing True or False values).

    Parameters
    ----------
    lons : numpy.ndarray
        Array containing the longitude values of the grid cell midpoints
    lats : numpy.ndarray
        Array containing the latitude values of the grid cell midpoints
    
    Returns
    -------
    continent_masks : xarray.Dataset
        Dataset containing masks for each continent
    '''
    
    # Read the continents shapefile.
    continents = gp.read_file(settings.general_data_directory + '/continents/continent.shp')
    
    # Create a mask for the continents, where each continent is assigned a unique number.
    continent_masks = regionmask.mask_geopandas(continents, lons[:-1], lats)
    continent_masks = continent_masks.rename({'lon': 'longitude', 'lat': 'latitude'})
    
    # Add columns (longitude values) to the continent mask to match the original grid
    column_to_add = xr.DataArray(data=np.atleast_2d(continent_masks.sel(longitude=0.00).values).T,
                                 coords=[('latitude', continent_masks['latitude'].values), ('longitude', [360.00])])
    continent_masks = continent_masks.combine_first(column_to_add)

    # Add a buffer layer to the continent mask to avoid edge effects in the following averaging step.
    n_layers = 2
    continent_masks = add_buffer_layers(continent_masks, n_layers)

    # Extract the continent mask for each continent
    continent_masks_dataset = xr.Dataset(coords={'latitude': continent_masks['latitude'].values, 'longitude': continent_masks['longitude'].values})

    # Loop over the continents
    for ii in range(6):

        # Assign the number 1 where the continent mask is equal to the continent number and 0 otherwise.
        # For the continents 5 and 6, which are Oceania and Australia, assign the same number.
        if ii < 5:
            continent_masks_dataset[continents['CONTINENT'][ii]] = xr.ones_like(continent_masks).where(continent_masks==ii, 0)
        else:
            continent_masks_dataset[continents['CONTINENT'][ii]] = xr.ones_like(continent_masks).where(np.logical_or(continent_masks==5, continent_masks==6), 0)

    # Calculate a smoothed continent mask by averaging over a 5x5 grid and then extract only the values for the original grid.
    # This is to exclude more cells around the continents when calculating the resulting mask.
    continent_masks_dataset_smooth = continent_masks_dataset.rolling(longitude=1+2*n_layers, latitude=1+2*n_layers, center=True).mean()
    continent_masks_dataset_smooth = continent_masks_dataset_smooth.sel(latitude=lats, longitude=lons)

    # Calculate the continent mask.
    continent_masks = continent_masks_dataset_smooth > 0

    # Reoder the continent variables.
    continent_masks = continent_masks[['North America', 'South America', 'Europe', 'Africa', 'Asia', 'Oceania']]

    return continent_masks


def get_histogram(geographic_masks, wind_power_density_max, wind_power_density, energy_deficits):
    '''
    Calcualte the 2D histogram showing the joint probability distribution of energy deficit and wind power density.

    Parameters
    ----------
    geographic_masks : xarray.Dataset
        Dataset containing the land, coast, and sea masks
    wind_power_density : xarray.DataArray
        Dataset containing the wind power density data
    energy_deficits : xarray.DataArray
        Dataset containing the energy deficit data
    
    Returns
    -------
    histogram_2d : xarray.Dataset
        Dataset containing the 2D histogram showing the joint probability distribution of energy deficit and wind power density
    '''

    # Get the max wind power density. Keep it the same for both seasonal and weather variability.
    rounded_max_wind_power_density = np.ceil(wind_power_density_max.values/10)*10

    # Define the bin edges of the 2D histogram.
    x_bin_edges = np.insert(np.logspace(-2.5, 0, 25),0,0)*rounded_max_wind_power_density
    y_bin_edges = np.linspace(0,1,26)

    # Get the latitude and longitude values of the wind resource.
    lons = wind_power_density['longitude'].values
    lats = wind_power_density['latitude'].values

    # Define an initial mask where all cells are included.
    non_redundant_cells = xr.DataArray(data=np.full((len(lats), len(lons)), True), dims=['latitude', 'longitude'],
                                       coords={'latitude': lats, 'longitude': lons}, name='cells_of_interest')

    # At latitude 90 and -90, keep only one value of longitude because the other values are redundant.
    non_redundant_cells.loc[dict(latitude=90, longitude=slice(lons[1], lons[-1]))] = False
    non_redundant_cells.loc[dict(latitude=-90, longitude=slice(lons[1], lons[-1]))] = False

    # At longitude 360, do not keep any values of latitude because they are redundant.
    non_redundant_cells.loc[dict(longitude=360)] = False

    # Initialize the 2D histogram dataset.
    histogram_2d = xr.Dataset()

    # Add the 2D histogram for each mask to the dataset.
    for mask in list(geographic_masks.data_vars):

        # Combine the mask with the land, coast, and sea masks.
        cells_of_interest = np.logical_and(non_redundant_cells, geographic_masks[mask])

        # Assign nan values to the cells that are not of interest.
        wind_power_density_of_interest = wind_power_density.where(cells_of_interest, np.nan)
        energy_deficits_of_interest = energy_deficits.where(cells_of_interest, np.nan)

        # Calculate the 2D histogram.
        if 'year' in energy_deficits_of_interest.dims:

            # Calculate the 2D histogram for each year and add the resulting occurrences together.
            create_new = True
            for year in energy_deficits['year']:
                if create_new:
                    histogram_2d_temp = histogram(energy_deficits_of_interest.sel(dict(year=year)),
                                                  wind_power_density_of_interest.sel(dict(year=year)),
                                                  bins=[y_bin_edges, x_bin_edges])
                    create_new = False
                else:
                    histogram_2d_temp = histogram_2d_temp + histogram(energy_deficits_of_interest.sel(dict(year=year)), # type: ignore
                                                                      wind_power_density_of_interest.sel(dict(year=year)),
                                                                      bins=[y_bin_edges, x_bin_edges])
        else:

            # Calculate the 2D histogram.
            histogram_2d_temp = histogram(energy_deficits_of_interest, wind_power_density_of_interest, bins=[y_bin_edges, x_bin_edges])

        # Normalize the 2D histogram.
        histogram_2d[mask] = histogram_2d_temp/histogram_2d_temp.sum() # type: ignore
    
    # Rename the dimensions of the 2D histogram.
    histogram_2d = histogram_2d.rename({energy_deficits.name+'_bin': 'energy_deficit', wind_power_density.name+'_bin': 'mean_wind_power_density'})
    
    return histogram_2d


def get_meadian_energy_deficit(histogram_2d):
    '''
    Calculate the median energy deficit for each wind power density bin.

    Source: https://math.stackexchange.com/questions/2591946/how-to-find-median-from-a-histogram

    Parameters
    ----------
    histogram_2d : xarray.Dataset
        Dataset containing the 2D histogram showing the joint probability distribution of energy deficit and wind power density

    Returns
    -------
    median_energy_deficit : xarray.Dataset
        Dataset containing the median energy deficit for each wind power density bin
    '''

    # Define the y-bin edges of the 2D histogram.
    y_bin_edges = np.linspace(0,1,26)

    # Initialize the median energy deficit dataset.
    median_energy_deficit = xr.Dataset()

    # Calculate the median energy deficit for each mask.
    for mask in list(histogram_2d.data_vars):

        # Calculate the total number of occurrences in each wind power density bin.
        total_occurrences_in_bin = histogram_2d[mask].sum(dim='energy_deficit')

        # Calculate the cumulative number of occurrences in each wind power density bin.
        cumulative_occurrences_in_bin = histogram_2d[mask].cumsum(dim='energy_deficit')

        # Find the bin that contains the median energy deficit, i.e., the bin where half of the total occurrences falls in between the lower and upper bin edges.
        median_group_in_bin = xr.full_like(histogram_2d[mask], False).where(cumulative_occurrences_in_bin < total_occurrences_in_bin/2, True).argmax(dim='energy_deficit')

        # Calculate the median energy deficit according to the formula in the description.
        median_energy_deficit[mask] = y_bin_edges[median_group_in_bin] + (((total_occurrences_in_bin/2 - cumulative_occurrences_in_bin[median_group_in_bin-1,:])/
                                                                           histogram_2d[mask][median_group_in_bin,:])*
                                                                           (y_bin_edges[median_group_in_bin+1] - y_bin_edges[median_group_in_bin]))

    return median_energy_deficit
    

def linear_regression(array):
    '''
    Calculate the linear regression of the array along the year dimension.
    
    Parameters
    ----------
    array : xarray.DataArray
        Array to calculate the linear regression of

    Returns
    -------
    regression : xarray.Dataset
        Dataset containing the slope, intercept, r-value, p-value, and standard error of the linear regression
    '''
    
    # Define a wrapper around scipy linregress to use in apply_ufunc
    def local_linear_regression(years, values_along_years):
        
        # Perform the linear regression.
        slope, intercept, r_value, p_value, std_err = stats.linregress(years, values_along_years)
        
        return slope, intercept, r_value, p_value, std_err
    
    # Apply a vectorized version of the scipy linregress function to the array.
    slope, intercept, r_value, p_value, std_err = xr.apply_ufunc(local_linear_regression, array['year'], array,
                                                                 input_core_dims=[['year'], ['year']],
                                                                 output_core_dims=[[], [], [], [], []],
                                                                 vectorize=True)
    
    # Create a dataset containing the results of the linear regression.
    regression = xr.Dataset({'slope': slope, 'intercept': intercept, 'r_value': r_value, 'p_value': p_value, 'std_err': std_err})

    return regression


def get_percentile_rank(array, mask):
    '''
    Calculate the percentile rank of each value in the array.

    Parameters
    ----------
    array : xarray.DataArray
        Array to calculate the percentile rank of
    mask : xarray.DataArray
        Array containing the mask to use when calculating the percentile rank

    Returns
    -------
    percentile_rank : xarray.DataArray
        Array containing the percentile rank of each value in the array
    '''

    # Create a temporary array where the values outside the mask are set to nan.
    array_temp = array.where(mask, np.nan)

    # Create a 1D array containing all the values in the temporary array.
    values_to_rank = array_temp.values.reshape(len(array_temp['latitude'])*len(array_temp['longitude']))

    # Calculate the percentile rank of each value in the temporary array.
    percentile_rank = stats.rankdata(values_to_rank, nan_policy='omit')*100/array_temp.count().values

    # Reshape the percentile rank array to match the original array.
    percentile_rank = xr.DataArray(data=percentile_rank.reshape(len(array_temp['latitude']),len(array_temp['longitude'])),
                                   coords=[('latitude', array_temp['latitude'].values), ('longitude', array_temp['longitude'].values)])

    return percentile_rank


def get_drought_extension(cell_areas, geographic_masks, continent_masks, wind_resource, maximum_energy_deficits_for_wind_droughts):
    '''
    Calculate the spatial extension of the most severe wind droughts in each continent and in each year.

    Parameters
    ----------
    cell_areas : xarray.DataArray
        Array containing the area of each cell defined by the lat/lon grid
    geographic_masks : xarray.Dataset
        Dataset containing the land, coast, and sea masks
    continent_masks : xarray.Dataset
        Dataset containing masks for each continent
    wind_resource : xarray.Dataset
        Dataset containing the wind speed and power density data
    maximum_energy_deficits_for_wind_droughts : xarray.DataArray
        Array containing the maximum energy deficit of wind droughts in each cell and in each year

    Returns
    -------
    drought_extension : xarray.Dataset
        Dataset containing the spatial extension of the most severe wind droughts in each continent
    '''

    # Define the wind resource of interest, i.e., the cells with a mean wind power density of at least 150 W/m^2.
    wind_resource_of_interest = wind_resource['mean_wind_power_density'].where(wind_resource['mean_wind_power_density'] >= 150, np.nan)

    # Select the cells that are land or coast.
    wind_resource_of_interest = wind_resource_of_interest.where(np.logical_or(geographic_masks['land'], geographic_masks['coast']), np.nan)

    # Initialize the drought extension dataset.
    drought_extension = xr.Dataset()

    # Loop over the continents.
    for continent in list(continent_masks.data_vars):

        # Select the cells in the continent of interest.
        region_of_interest = wind_resource_of_interest.where(continent_masks[continent], np.nan)

        # Select the energy deficits in the cells of interest and assign them a value of 1.
        energy_deficits_of_interest = xr.ones_like(maximum_energy_deficits_for_wind_droughts).where(region_of_interest.notnull(), np.nan)

        # Calculate the total extension of the most severe wind droughts in the continent of interest for each year.
        total_extension = (cell_areas*energy_deficits_of_interest).groupby('year').sum(skipna=True)

        # Calculate the total area of the region of interest.
        total_continent_area = cell_areas.where(region_of_interest.notnull(), np.nan).sum(skipna=True)

        # Calculate the fraction of the continent of interest that is affected by the most severe wind droughts in each year.
        drought_extension[continent] = total_extension / total_continent_area
    
    return drought_extension
