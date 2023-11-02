# Global wind droughts
This repository contains scripts to reproduce the results of the following paper:

Enrico G. A. Antonini, Edgar Virgüez, Sara Ashfaq, Lei Duan, Tyler H. Ruggles, Ken Caldeira, "Global distribution of and trends in wind droughts", Communications Earth & Environment, 2023.

If you need a copy of the paper or have any questions about the scripts, please email me at enrico.antonini@eiee.org.

## Prerequisites

To run the scripts in this repository, you need to have:

* [Conda](https://docs.conda.io/en/latest/), which is an open source package management system and environment management system that comes with [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/).
* [CDS API key](https://cds.climate.copernicus.eu/api-how-to), which is needed to download the climate data from the Copernicus' Climate Data Store.
* 6 TB of free disk space. 3 TB are needed to download the climate data and 3 TB are needed to process the climate data and calculate the wind power densities and energy deficits.
* Sufficient RAM memory to process the climate data. Downloaded climate files each containing a single variable for a single year are about 17 GB. Other derived files have dimensions up to 68 GB. I would recommend at least 128 GB of RAM memory. Some steps were intentionally made supre basic to minimize the RAM memory usage.

## Conda environments

The folder [environments](https://github.com/eantonini/Global_wind_droughts/tree/main/environments) contains yaml files to create conda environments for

* downloading climate data from the Copernicus' Climate Data Store,
* processing the climate data and calculate the wind power densities and energy deficits,
* analyzing the wind power densities and energy deficits and generating the figures of the paper.

A conda environment can be created by running the following command:
```
conda env create --file selected_environment.yml
```
where you need to specify the actual yaml file name.

## Climate data download

[Required conda environemt: `download_climate_data`]

To download the climate data from the Copernicus' Climate Data Store, you need to run the following command from within the [script_and_data](https://github.com/eantonini/Global_wind_droughts/tree/main/scripts_and_data) folder:
```
python step0_download_data.py
```
The download will take several hours, depending on your internet connection, and will require about 3 TB of disk space.

## Climate data processing

[Required conda environemt: `wind_droughts_processing`]

To process the climate data and calculate the wind power densities and energy deficits, you need to run the following commands from within the [script_and_data](https://github.com/eantonini/Global_wind_droughts/tree/main/scripts_and_data) folder.

### Step 1

Get the wind speed time series in a given year and its annual mean:
```
python step_1_get_wind_speed_and_annual_mean_in_year.py year
```
where you need to specify the actual year.

### Step 2

Get the wind power density time series in a given year and its annual mean:
```
python step_2_get_wind_power_density_and_annual_mean_in_year.py year
```
where you need to specify the actual year.

### Step 3

Get the mean wind speed and mean wind power density across all years:
```
python step_3_get_mean_wind_speed_and_power_density.py
```

### Step 4

Initialize the climatological wind power density time series file for subsequent steps:
```
python step_4_create_climatological_wind_power_density_array.py
```

### Step 5

Add the wind power density time series in a given year to the climatological wind power density time series:
```
python step_5_add_wind_power_density_in_year.py year
```
where you need to specify the actual year.

### Step 6

Divide the climatological wind power density time series by the number of years:
```
python step_6_divide_climatological_wind_power_density.py
```

### Step 7

Get the energy deficit for the seasonal variability:
```
python step_7_get_energy_deficit_for_seasonal_variability.py
```

### Step 8

Get the energy deficit for the weather variability in a given year:
```
python step_8_get_energy_deficit_for_weather_variability_in_year.py year
```
where you need to specify the actual year.

### Step 9

Get the energy deficit for wind droughts:
```
python step_9_get_energy_deficit_for_wind_droughts_in_year.py year
```
where you need to specify the actual year.

Note that step 8 and 9 generate energy deficit for a specific year and save it to a single file. All the generated files should then be combined into a single file.

## Results analysis

[Required conda environemt: `wind_droughts_analysis`]

To analyze the wind power densities and energy deficits and generate the figures of the paper, you need to run the following command from within the [script_and_data](https://github.com/eantonini/Global_wind_droughts/tree/main/scripts_and_data) folder:
```
python step_10_analyze_results.py
```