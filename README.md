# [pfra-hydromet](https://dewberry.github.io/pfra-hydromet/)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Dewberry/pfra-hydromet/master)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

---

# Description

__pfra-hydromet__ is a collection of tools for developing excess rainfall scenarios
for input to hydraulic models using:
  1. Meteorological data
  2. Random sampling
  3. Hydrologic transform
  4. Convolution algorithm for  grouping

These tools ([jupyter notebooks](https://jupyter.org/) ) ingest data from the NOAA Hydrometeorological Design Studies Center ([HDSC](https://www.nws.noaa.gov/oh/hdsc/index.html)) and return unique, weighted runoff events suitable for use in 2D hydraulic *rain-on-grid* models. Executed notebooks should be saved as documentation of the inputs, outputs, and results for a given project location.

__NOTE__: [EventsTable](EventsTable.ipynb) is currently the primary notebook for developing excess precipitation scenarios, and is managed and called by [*Papermill*](https://pypi.org/project/papermill/). Within this repo, papermill is designed to act as a *manager* to maintain consistency in computation, and ensure cells are executed in order. *Manager* notebooks are designated with the `PM-` prefix.

---

## Contents

##### Notebooks:

- [__PrecipTable__](PrecipTable.ipynb): Retrieve NOAA Atlas 14 precipitation statisics for an Area of Interest (AOI).

- [__PM-EventsTable__](PM-EventsTable.ipynb): Manager notebook that executes `EventsTable` and/or `reEventsTable`.

- [__EventsTable__](EventsTable.ipynb): Calculates excess rainfall using area-averaged NOAA Atlas 14 precipitation data, temporal distributions, and the curve number (CN)* transform. The output is a set of unique, weighted excess rainfall time series.

- [__reEventsTable__](reEventsTable.ipynb): Calculates the reduced excess rainfall given a user-specified stormwater removal rate and capacity. 

- [__distalEventsTable__](distalEventsTable.ipynb): Calculates excess rainfall using updated randomized curve numbers and the original precipitation events calculated in `EventsTable.ipynb`. The events are combined using the groups determined during the convolution steps in `EventsTable.ipynb`. The `reEventsTable` notebook can be then be executed in order to calculate the reduced excess rainfall.

- [__MetadataExplorer__](MetadataExplorer.ipynb): Explores the metadata file created by `PM-EventsTable` or `distalEventsTable` during the excess rainfall calculations.

- [__Convolution_Parameters__](Convolution_Parameters.ipynb): Describes the test statistic and parameters used during the convolution step in the `EventsTable` notebook.


##### DataRepository:

- __Temporal_Distributions__: Folder containing csv files of temporal distributions of observed rainfall patterns broken down by volume, region, duration, and quartile [NOAA Published](https://hdsc.nws.noaa.gov/hdsc/pfds/pfds_temporal.html). Note that the original data were compiled into csv's for uniform formatting.

- __Temporal_Distributions_Plots__: Folder containing a Jupyter Notebook for each NOAA Atlas 14 volume with the plotted temporal distributions for each region, duration, and quartile.

- `NEH630_Table_10_1.json`: A formatted copy of Table 10-1 from the National Engineering Handbook [Chapter 10].(https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/NEHhydrology/ch10.pdf.) which lists the CN values for dry and wet antecedent moisture conditions.

- `NOAA_Atlas_Volume_Codes.json`: Metadata that maps the NOAA Atlas 14 volume number to the volume code. [Source](https://hdsc.nws.noaa.gov/hdsc/pfds/pfds_gis.html)

- `NOAA_Temporal_Areas_US.geojson`: geojson file containing the vector ploygons of the NOAA Atlas 14 temporal distribution areas. This file was constructed using the individual vector ploygons for each volume. [Source](https://hdsc.nws.noaa.gov/hdsc/pfds/pfds_temporal.html)

- `Temporal_Distribution_Data_Map.json`: Metadata used to extract the temporal distribution data from the csv files saved within the __Temporal_Distributions__ folder.

- `Temporal_Quartile_Ranks.xlsx`: Excel Workbook that contains the percentage of precipitation events whose temporal distributions are represented by those in each quartile of a specific volume/region/duration. [Source](https://www.nws.noaa.gov/oh/hdsc/currentpf.html)


*The ([CN Method](https://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb1044171.pdf)) is currently the only transform method in use for this project. Other transforms are available and can be adopted into the tool with minor modifications.

---

## Workflow

1. Run [PrecipTable](PrecipTable.ipynb) in order to calculate the area-averaged precipitation frequency table for the specified durations as well as to determine the NOAA Atlas 14 volume and region.
    ```
      Inputs:
        1. A vector polygon of the area of interest
        2. Optional/as needed: 
            - Precipitaiton event durations, the standard are 6, 12, 24, and 96 hour.
            - The polygon's projection if it can not be determined automatically
      Outputs:
        1. A spreadsheet with the area-averaged precipitation frequency table for each duration and a sheet with the NOAA Atlas 14 volume and region numbers.
    ```
    
    
2. Run [PM-EventsTable](PM-EventsTable.ipynb) Create runoff time-series events for hydraulic simulation

    ```
      Inputs:
        1. Precipitation spreadsheet (calculated in step 1)
        2. Curve number for the AOI
        3. NOAA Atlas 14 volume number
        4. Storm durations
        5. Filenames and paths for outputs
        6. EventsTable.ipynb

      Outputs:
        1. Precipitation statistics for each duration
        2. HTML copy of notebook
    ```



---

### Documentation

Complete project documentation can be found in [read the docs](https://dewberry.github.io/pfra-hydromet/about/).
