# Global wind droughts
This repository contains scripts to reproduce the results of the following paper:

Enrico G. A. Antonini, Edgar Virgüez, Sara Ashfaq, Lei Duan, Tyler H. Ruggles, Ken Caldeira, "Global distribution of and trends in wind droughts", Communications Earth & Environment, 2023.

If you need a copy of the paper or have any questions about the scripts, please email me at enrico.antonini@eiee.org.

## How to use the files in this repository

### Prerequisites

To run the scripts in this repository, you need to have:

* [Conda](https://docs.conda.io/en/latest/), which is an open source package management system and environment management system that comes with [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/).
* [CDS API key](https://cds.climate.copernicus.eu/api-how-to), which is needed to download the climate data from the Copernicus' Climate Data Store.
* 6 TB of free disk space. 3 TB are needed to download the climate data and 3 TB are needed to process the climate data and calculate the wind power densities and energy deficits.
* Sufficient RAM memory to process the climate data, where each file containing a single variable for a single year is about 18 GB.

### Conda environments

The folder [environments](https://github.com/eantonini/Global_wind_droughts/tree/main/environments) contains yaml files to create conda environments for
* downloading climate data from the Copernicus' Climate Data Store,
* processing the climate data and calculate the wind power densities and energy deficits,
* analyzing the wind power densities and energy deficits and generating the figures of the paper.

A conda environment can be created by running the following command:
```
conda env create --file selected_environment.yml
```
where you need to specify the actual yaml file name.

### Climate data download

To download the climate data from the Copernicus' Climate Data Store, you need to run the following command from within the [script_and_data](https://github.com/eantonini/Global_wind_droughts/tree/main/scripts_and_data) folder:

```
python step0_download_data.py
```

The download will take several hours, depending on your internet connection, and will require about 3 TB of disk space.

### Climate data processing

To process the climate data and calculate the wind power densities and energy deficits, you need to run the following command from within the [script_and_data](https://github.com/eantonini/Global_wind_droughts/tree/main/scripts_and_data) folder:


