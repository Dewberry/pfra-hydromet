[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Dewberry/pfra-hydromet/master)

# Description

__pfra-hydromet__ is a collection of tools for developing excess precipitation scenarios
for input to hydraulic models using:
  1. Meteorological data
  2. Hydrologic transform
  3. Random sampling
  4. Convolution algorithm for  grouping

These tools ([jupyter notebooks](https://jupyter.org/) ) ingest data from the NOAA Hydrometeorological Design Studies Center ([HDSC](https://www.nws.noaa.gov/oh/hdsc/index.html)) and return unique, weighted runoff events suitable for use in 2D hydraulic *rain-on-grid* models. Executed notebooks should
be saved as documentation of the inputs, outputs, and results for a given project location.

__NOTE__: [EventsTable](EventsTable.ipynb) is currently the primary notebook for developing
excess precipitation scenarios, and is managed and called by [*Papermill*](https://pypi.org/project/papermill/).
Within this repo, papermill is designed to act as a *manager* to maintain consistency in computation,
and ensure cells are executed in order. *Manager* notebooks are designated with the `PM-` prefix.


## Contents

---

##### Notebooks

1. [__PrecipTable__](PrecipTable.ipynb): Retrieve NOAA Atlas 14 precipitation statisics
at an Area of Interest (AOI).

2. [__PM-EventsTable__](PM-EventsTable.ipynb): Manager notebook that executes `EventsTable` and executes.

3. [__EventsTable__](EventsTable.ipynb): Calculates excess rainfall using the NOAA Atlas 14 mean precipitation data, temporal distributions, and the curve number (CN)* transform. The output is a set of unique, weighted
excess precipitation time series.

##### DataRepository

 - __Temporal datasets__ (csv): Temporal distributions of observed rainfall patterns broken down by volume, region, duration, and quartile [NOAA Published](https://hdsc.nws.noaa.gov/hdsc/pfds/pfds_temporal.html) are saved as csv's here. Datasets here were compiled from source files into csv format for uniform formatting.

- `Temporal_Distribution_Data_Map.json` Mapping data used to extract the temporal distributions from
regional datasets to a uniform csv format.

- `Temporal_Quartile_Ranks.xlsx` contains the percentage of precipitation events whose temporal distributions are represented by those in each quartile. [Source](https://www.nws.noaa.gov/oh/hdsc/currentpf.html).

- `NEH630_Table_10_1.json` is a formatted copy of table 10-1 from the National Engineering
Handbook [Chapter 10](https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/NEHhydrology/ch10.pdf.) listing
CN values over differing antecedent moisture conditions.

*The ([CN Method](https://www.nrcs.usda.gov/Internet/FSE_DOCUMENTS/stelprdb1044171.pdf))
is currently the only transform method in use for this project. Other transforms are available
and can be adopted into the tool with minor modifications.


## Usage

---

1. Calculate area averaged precipitation frequencies for an area of interest:

  - Run [PrecipTable](PrecipTable.ipynb).

    ```
      Inputs:
        1. A vector polygon for the area of interest (AOI)
        2. NOAA Atlas 14 volume number
        3. Storm durations

      Outputs:
        1. Precipitation spreadsheet with statistics for each duration

    ```
2. Create runoff time-series events for hydraulic simulation

  - Run [PM-EventsTable](PM-EventsTable.ipynb).

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

Complete project documentation can be found in [read the docs](https://dewberry.github.io/pfra-hydromet/).
