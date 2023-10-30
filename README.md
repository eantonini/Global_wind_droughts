# Global wind droughts
This repository contains scripts to reproduce the results of the following paper:

Enrico G. A. Antonini, Edgar Virguez, Sara Ashfaq, Lei Duan, Tyler H. Ruggles, Ken Caldeira, "Global distribution of and trends in wind droughts", Communications Earth & Environment, 2023.

If you need a copy of the paper or have any questions about the scripts, please email me at enrico.antonini@eiee.org.

## How to use the files in this repository

### Conda environments

The folder [environments](https://github.com/eantonini/Global_wind_droughts/tree/main/environments) contains yaml files to create conda environments for
* downloading climate data from the Copernicus's Climate Data Store,
* processing the climate data and calculate the wind power densities and energy deficits,
* analyzing the wind power densities and energy deficits and generating the figures of the paper.

A conda environment can be created by running the following command:
```
conda env create --file selected_environment.yml
```
where you need to specify the actual yaml file name.

