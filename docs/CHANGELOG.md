Changelog
=========

0.11.0
------

**General**
 - Binder build tested and linked to readme. 

**Tests** 
 - Added unit test example for core.hydromet.


0.10.0
------

*Public Release*

Tools in this repository were developed by Dewberry in its role as part of [STARR II](http://www.starr-team.com/starr/Pages/default.aspx), a FEMA Production and Technical Services provider,
to facilitate the automation of hydrologic computations for probabilistic flood risk studies.

Computational approaches
were developed in coordination with engineers at [COMPASS](https://www.aecom.com/press-releases/aecom-announced-today-that-a-joint-venture-it-leads-has-been-awarded-a-contract-with-a-ceiling-of-us600-million-from-the-u-s-department-of-homeland-securitys-federal-emergency-management-ag/),
and guidance from engineers and scientists at
[FEMA](https://www.fema.gov/),
[USACE](https://www.usace.army.mil/), [USGS](https://www.usgs.gov/), and [NOAA](https://www.noaa.gov/).

**PrecipTable** (Notebook)

-  Download gridded data from [NOAA Atlas 14](https://hdsc.nws.noaa.gov/hdsc/pfds/) PFDS:
  - User provides a vector polygon for Area of Interest (AOI)
  - User provides selected storm durations
- Calculate average values within an AOI for:
  - Mean precipitation
  - Upper & lower confidence limits


**EventsTable** (Notebook)

-  Execute random event sampling using normal and lognormal distributions
-  Randomly pair precipitation totals (volume) with temporal distributions (shape) to create a suite of events
-  Performs convolution calculations to group like events based on volume and shape relationships
-  Develop probabilities for paired events
-  Write intermediate and metadata files to disk

**Core**
-  Codebase for operations performed in notebooks.
