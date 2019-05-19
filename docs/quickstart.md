
## Overview of tools:

Example input and output parameters/files.

## Example Inputs


`filenames and paths` for the precipitation frequency table, DataRepository, and location to save the outputs.

`Vector File` covering area of interest (hydraulic domain).

`NOAA data` Atlas 14 Volume number, region, and duration.

`CN` Curve Number calculated for the area of interest.

#### Convolution Parameters (*epsilons)

`tempEpsilon` The number of hours over which to resample the incremental excess rainfall during the first convolution.


`tempEpsilon2` The number of hours over which to resample the incremental excess rainfall during the final convolution.


`convEpsilon` The maximum allowable percent difference in incremental excess between two event at any given resampled time step.


`volEpsilon` The maximum allowable percent difference between the two events total runoff volume.

**Suggested epsilon values are given. Epsilon values should be adapted for each project area by analyzing the curve fitting
results. Plots can be shown for every group by setting `display_plots=True`*


#### Optional Features

Check functions for default options for using a known seed, testing, debugging, etc. (e.g. `seed=88`)



## Example Outputs


### PrecipTable
Precipitation Frequency estimates for AOI.


| Tr |  Lower 90%  | Expected Value |  Upper 90% |  
|:-:|:-:|:-:|:-:|
|2| 5.25    | 6.32 | 7.60 |
|5| 6.85  | 8.27| 9.96 |
|10|  8.20 | 9.95 | 12.05|
|25|  10.06  | 12.45 | 15.86 |
|50|  11.48 | 14.56 | 18.75 |
|100| 12.85  | 16.86 | 22.22 |
|200| 14.20  | 19.33 | 26.19 |
|500| 16.20  | 22.89 | 31.84 |
|1000| 17.71| 25.79 | 36.11 |   


### Random Events Table
Randomly chosen events & metadata before convolution.

|Return Period| Ann. Exc. Prob. |ARI |Log10_ARI |  Expected Value | Lower (90%) |Upper (90%) |Quartile  |  Sigma  | Fitted Lower (90%) Limit  |  Fitted Upper (90%) Limit  |  Random Precipitation|
|:-:|:-:|:-:|:-:|:-:|:-:|:-: |:-:|:-:| :-:|:-:| :-:|
|2.000240382|0.499939912| 1.442945202|0.366686304 |3.15751785  |2.876612547 |3.507262528 |3   |0.078317368| 2.855989829| 3.49088042 | 3.308578113|
|2.001189958|0.499702687| 1.44393337 |0.367370897 |3.157996455 |2.877045943 |3.507791811 |3   |0.078224163| 2.856763942| 3.490992541| 3.331452983|
|2.00267425 |0.49933233 | 1.445477879|0.368439979 |3.158744005 |2.877722879 |3.508618516 |3   |0.078224124| 2.857440325| 3.491818745| 3.453517084|
|...|...| ...|... |... |... |... |... |...| ...| ...| ...|
|923.2492014|0.001083131| 922.749111 |6.827357379 |13.55155136 |11.4917289  |14.67457848 |3   |0.087869372| 12.10830888| 15.16682025| 12.88399871|
|1119.173615|0.000893516| 1118.67354 |7.019898923 |14.11398151 |11.92280701 |15.27555124 |3   |0.088696172| 12.59748493| 15.81303532| 14.49259104|
|3000       |0.000333333| 2999.499972|8.006200878 |17.38236261 |14.3978362  |18.76211492 |3   |0.092819505| 15.43292272| 19.57804981| 18.76211492|


### Curve Number Table
Randomly chosen CN number and metadata.

|Random Sample| Lower  | Expected Value  |Upper|   Sigma |  Fitted Lower Limit| Fitted Upper Limit|  Random CN|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|1       |67  |83  |93  |0.116325785 |71.50469526 |96.34332368 |88|
|2       |67  |83  |93  |0.116325785 |71.50469526 |96.34332368 |89|
|3       |67  |83  |93  |0.116325785 |71.50469526 |96.34332368 |93|
|...|...|...|...|...|...|...|...|
|5107    |67  |83  |93  |0.116325785 |71.50469526 |96.34332368 |90|
|5108    |67  |83  |93  |0.116325785 |71.50469526 |96.34332368 |87|
|5109    |67  |83  |93  |0.116325785 |71.50469526 |96.34332368 |93|



### Final Events Table
Final events to be modeled.

|hours|E0001|E0002|E0003|...|E0395|E0396|E0397|E0398|E0399|  
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|0  |0  |0| 0|...|  0|  0|  0|  0|  0|    
|0.5|0.010| 0.0047|0.0030|...|0.01359|0.01038|0.00674|0.00234|0.00555|
|1  |0.013| 0.00948|    0.0043| ...|    0.0447| 0.0584| 0.0124| 0.00356|    0.008906|
|1.5    |0.10|  0.067680|0.0294|...|    0.1660| 0.15102|    0.1220| 0.0198| 0.0433|
|...|...|...|...|...|...|   ...|...|...|...|
|22 |0.13|  0.1534| 0.0807| ...| 0.1878| 0.1707|0.234|   0.07389|    0.149328|
|23.5   |0.103| 0.21576|    0.14100| ...|0.1373|  0.1408| 0.317|  0.15162|    0.221391|
|24|0.08732|    0.25004|    0.1836| ...|0.1132| 0.1293| 0.3403| 0.2053| 0.2430|


### Event Weight Table

| |Weight|
|:-:|:-:|
|E0001 |0.006207654|
|E0002 |0.005313752|
|...|...|
|E0964 |0.000019172|
|E0965 |0.000071544|


### Metadata

Additional data developed in intermediated calculations included for traceability/reproducibility.


```json
[{
    # Run Data (Event_Duration_Quartile_Decile_CurveNumber)
    'RunInfo': 

            {'E60001': 'E1_6Hr_Q1_D10_CN79',
             'E60002': 'E2_6Hr_Q1_D90_CN78',
             ...}

    # Cumulative (raw) precipitation time-series
    'precip':  
    
            {'E60001': {'0.0': 0.0,
                        '0.5': 0.8010830472988398,
                        '1.0': 1.4541398793359375,
                        '1.5': 1.8372665541310351,
                        '2.0': 2.057129020916858,
                        '2.5': 2.1550875457224223,
                        '3.0': 2.170325538469955,
                        '3.5': 2.176856106790326,
                        '4.0': 2.176856106790326,
                        '4.5': 2.176856106790326,
                        '5.0': 2.176856106790326,
                        '5.5': 2.176856106790326,
                        '6.0': 2.176856106790326},
             'E60002': { ...
             }


    # Cumulative excess-precipitation time-series
    'cum_excess':

            {'E60001': {'0.0': 0.0,
                        '0.5': 0.02479673948334961,
                        '1.0': 0.23766036959891948,
                        '1.5': 0.4300482251417235,
                        '2.0': 0.5562285716308201,
                        '2.5': 0.6155457935761868,
                        '3.0': 0.6249312656222388,
                        '3.5': 0.6289662984090281,
                        '4.0': 0.6289662984090281,
                        '4.5': 0.6289662984090281,
                        '5.0': 0.6289662984090281,
                        '5.5': 0.6289662984090281,
                        '6.0': 0.6289662984090281},
             'E60002': { ...
             }

    # Incremental excess-precipitation time-series
    'incr_excess'

            {'E60001': {'0.0': 0.0,
                        '0.5': 0.02479673948334961,
                        '1.0': 0.21286363011556986,
                        '1.5': 0.19238785554280402,
                        '2.0': 0.1261803464890966,
                        '2.5': 0.059317221945366705,
                        '3.0': 0.009385472046052001,
                        '3.5': 0.004035032786789294,
                        '4.0': 0.0,
                        '4.5': 0.0,
                        '5.0': 0.0,
                        '5.5': 0.0,
                        '6.0': 0.0},

             'E60002': { ...
             }

    # Curve groups (similar excess-precipitation time-series)
    'groups'

        {'E0001': ['E60002',
                   'E60930',
                   'E60551',
                      ...  ,
                    E64421',
                    E64423',
                   'E64483'],
         'E0002': [...
         ]
        }


    # Test Statistic 
    'test_stat'
      {'E0001': [0.235358,
                 0.251304,
                 0.214125,
                    ...
                 0.18715,
                 0.184935,
                 0.174473],
       'E0002': [...
     ]
}]
```
