# core.hydromet

## parse_filename
```python
parse_filename(zip_name:str, reg:str) -> dict
```
Builds a dictionary with the region, recurrance interval, duration,
and statistic type using the zip_name and region.

## get_masked_mean_atlas14
```python
get_masked_mean_atlas14(gdf:'GeoDataFrame', raster:str) -> float
```
Masks the Atlas 14 precipitation raster by the passed polygon and then
calculates the average precipitation for the masked polygon.

## get_input_data
```python
get_input_data(precip_table_dir:str, duration:int, lower_limit:int=2, display_print:bool=True) -> pandas.core.frame.DataFrame
```
Extracts the precipitation frequency data for the specified duration
from an Excel sheet and returns the dataframe with the data.

## get_temporal_map
```python
get_temporal_map(data_dir:str, filename:str, vol:int, reg:int, dur:int, display_print:bool=True) -> dict
```
Reads the json file containing the temporal distribution data metadata
and returns the data map and number of rows to skip for the specified
volume, region, and duration.

## get_temporals
```python
get_temporals(temporal_dir:str, vol:int, reg:int, dur:int, qmap:dict, display_print:bool=True) -> pandas.core.frame.DataFrame
```
Reads the csv file containing the temporal distributions for the
specified volume, region, and duration. Rows with NaNs for an index
are dropped. Data was downloaded from:
https://hdsc.nws.noaa.gov/hdsc/pfds/pfds_temporal.html

## get_quartile_rank
```python
get_quartile_rank(data_dir:str, filename:str, vol:int, reg:int, dur:int, display_print:bool=True) -> list
```
Extracts the quartile ranks for the specified volume, region, and
duration. The quartile rank corresponds to the percentage of
precipitation events whose temporal distributions are represented
by those in a specific quartile.

## get_duration_weight
```python
get_duration_weight(data_dir:str, filename:str, vol:int, reg:int, dur:int, display_print:bool=True) -> list
```
Extracts the duration weight for the specified volume, region, and
duration. The duration weight corresponds to the percentage of
precipitation events with the specified duration.

## get_CN_distribution
```python
get_CN_distribution(data_dir:str, filename:str, CN:int, display_print:bool=True) -> dict
```
Open the json file containing the curve number values for different
antecedent moisture conditions and return the values for the
specified curve number.

## extrap_add_ari
```python
extrap_add_ari(df:pandas.core.frame.DataFrame, display_print:bool=True) -> pandas.core.frame.DataFrame
```
Calls the add_ari function to update the dataframe and
then calls the extrapolate_extremes function in order to extrapolate
the confidence limits and expected value of the precipitation amount
for the 2000 and 3000 year return periods.

## add_ari
```python
add_ari(df:pandas.core.frame.DataFrame) -> pandas.core.frame.DataFrame
```
Calculates the annual exceedance probability (AEP),
average recurrance interval (ARI), and log of the ARI and adds the
results to the original dataframe.

## extrapolate_extremes
```python
extrapolate_extremes(df:pandas.core.frame.DataFrame, rp:int, ycol:str) -> float
```
Extrapolates the ycol for the specified return period.

## generate_random_samples
```python
generate_random_samples(samplesize:int, seed:int=None, display_print:bool=True) -> pandas.core.frame.DataFrame
```
Selects the specified number of random samples from a continuous
normal distribution, calculates the inverse of the sample, and saves
the results in a dataframe with column "Tr", where "Tr" is the
recurrance interval.

## Truncate_Random_Events
```python
Truncate_Random_Events(r_events:pandas.core.frame.DataFrame, lower_limit:int=2, upper_limit:int=3000) -> pandas.core.frame.DataFrame
```
Removes events with recurrance intervals less than the lower_limit
(typically 2 years) and sets recurrance intervals greater than the
upper limit (typically 3000 years) eqaul to the upper limit.

## events_table_random
```python
events_table_random(raw_precip:pandas.core.frame.DataFrame, events_table:pandas.core.frame.DataFrame) -> pandas.core.frame.DataFrame
```
Calls the add_ari function to update the dataframe and then calls the
scipy_interp function in order calculate the expected value, lower
(90%) confidence limits, and upper (90%) confidence limits for the
events_table given the raw_precip dataframe.

## scipy_interp
```python
scipy_interp(raw_precip:pandas.core.frame.DataFrame, df:pandas.core.frame.DataFrame, ynew:str='Expected Value') -> pandas.core.frame.DataFrame
```
Interpolates the ynew values for the passed df given the Log10_ARI
and ynew valuea contained within the raw_precip dataframe.

## find_optimal_curve_std
```python
find_optimal_curve_std(df:pandas.core.frame.DataFrame, lower:str='Lower (90%)', upper:str='Upper (90%)', sdev:float=0.15) -> pandas.core.frame.DataFrame
```
Calculates/optimizes the standard deviation of the lognormal
distribution using the expected value, lower confidence limit/value,
and the upper confidence limit/value. The sum of the squared residuals
of the lower and upper confidence limits/values is used as the test
statistic (this statistic is minimized). Note that the sdev is the
initial estimate of the standard deviation. The fitted values should
be compared to the lower and upper confidence limits/values to
validate the optimization. Note: additional code exists at the end of
the script containing this function which can be edited in order to
improve the fit of the standard devation for CN.

## RandomizeData
```python
RandomizeData(df:pandas.core.frame.DataFrame, number:int, results_dir:str, AOI:str, duration:int=24, quartile:int=None, seed:int=None, sampling_distro:str='Lognorm', variable:str='Precipitation', lower:str='Lower (90%)', upper:str='Upper (90%)', plot:bool=False, display_print:bool=True) -> pandas.core.frame.DataFrame
```
Randomly selects a value (precipitation or curve number) from the log-
normal distribution given the expected value and optimized standard
devation for each recurrance interval/event.

## join_rdata_tables
```python
join_rdata_tables(rdata_tables:list, type:str, display_print:bool=True) -> pandas.core.frame.DataFrame
```
Concatenates the dataframe elements of the passed list producing a
single dataframe. This resulting dataframe's index is set from 1 to
the length of the dataframe.

## get_quartiles
```python
get_quartiles(raw_temporals:pandas.core.frame.DataFrame, dur:int, qrank:list, qmap:dict, vol:int, reg:int, plot:bool=False) -> dict
```
For each quantile, extract the temporal data from the raw_temporals
dataframe, convert the data to numeric, store the data in a dictionary,
and plot the deciles.

## map_quartiles_deciles
```python
map_quartiles_deciles(n_samples:int=75, seed:int=None, plot:bool=False, display_print:bool=True) -> pandas.core.frame.DataFrame
```
Constructs a dataframe containing randomly selected deciles for the
specified number of samples (events).

## prep_cn_table
```python
prep_cn_table(CN:int, arc_data:dict) -> pandas.core.frame.DataFrame
```
Constructs a dataframe with the average/expected curve number (CN),
the dry/lower CN, and the wet/upper CN. The dry, average, and wet
curve numbers refer to different antecedent runoff conditions, which
were obtained from NEH Part 630, Chapter 10, Table 10-1
(https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/NEHhydrology/ch10.pdf)

## populate_event_precip_data
```python
populate_event_precip_data(random_cns:pandas.core.frame.DataFrame, temporals:pandas.core.frame.DataFrame, random_precip_table:pandas.core.frame.DataFrame, data_table:pandas.core.frame.DataFrame, curve_group:dict, dur:int=24, adjust_CN_less24:bool=False) -> (<class 'pandas.core.frame.DataFrame'>, <class 'pandas.core.frame.DataFrame'>, <class 'pandas.core.frame.DataFrame'>, <class 'pandas.core.frame.DataFrame'>)
```
Calculates cumulative and incremental runoff for each event using a
randomly selected precipitation amount, quartile specific temporal
distribution, and curve number.

## update_CN
```python
update_CN(CN:int, duration:int, grid_avg_precip:float) -> (<class 'int'>, <class 'float'>, <class 'float'>)
```
Adjusts the curve number (CN), potential maximum retention after
runoff begins (S), and intial abstraction (Ia) for durations less than
24 hours. Contact Kaveh Zomorodi: kzomorodi@Dewberry.com for
additional details regarding the adj_CN equation.

## S_24hr
```python
S_24hr(CN:int) -> float
```
Calculates the potential maximum retention after runoff begins (S), in
inches.

## IA_24hr
```python
IA_24hr(s24:float) -> float
```
Calculats the inital abstraction (Ia) as a function of the maximum
potentail rention (S). Lim et al. (2006) suggest that a 5% ratio of
Ia to S is more appropriate for urbanized areas instead of the more
commonly used 20% ratio
(https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1752-1688.2006.tb04481.x).

## QCN_24hr
```python
QCN_24hr(grid_avg_precip:float, s24:float) -> float
```
Calculates runoff using equation 10-11 of NEH Part 630, Chapter 10
(https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/NEHhydrology/ch10.pdf).

## infiltration_24hr
```python
infiltration_24hr(grid_avg_precip:float, s24:float, qcn_24:float) -> float
```
Calculates the actual retention (or infilitration) after runoff
begins, in inches using equation 10-7 of NEH Part 630, Chapter 10
(https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/NEHhydrology/ch10.pdf).

## calculate_excess
```python
calculate_excess(precip:float, ia:float, s:float) -> float
```
Calculates runoff using the curve number approach. See equation 10-9
of NEH 630, Chapter 10
(https://www.wcc.nrcs.usda.gov/ftpref/wntsc/H&H/NEHhydrology/ch10.pdf)

## adjust_incremental
```python
adjust_incremental(raw:pandas.core.series.Series, excess:pandas.core.series.Series) -> pandas.core.series.Series
```
Calculates the incremental runoff depth (depth/timestep) using the
cumulative_to_incremental function, and then redistributes the first
non-zero incremental runoff value over the prior timesteps using the
incremental precipitation as a weighting function.

## cumulative_to_incremental
```python
cumulative_to_incremental(vector:pandas.core.series.Series) -> pandas.core.series.Series
```
Converts the cumulative depth (precipitation or runoff) into the
incremental depth, i.e. the depth/timestep (rate).

## convert_tempEpsilon
```python
convert_tempEpsilon(tempEpsilon:float, incr_excess:pandas.core.frame.DataFrame) -> int
```
Converts the tempEpsilon from the number of hours to the number of
corresponding timesteps.

## bin_sorting_dev
```python
bin_sorting_dev(incr_excess:pandas.core.frame.DataFrame, nbins:int, display_print:bool=True) -> list
```
Computes the histogram of the series data with the specified number
of bins and returns the results as a list.

## get_bin_slice
```python
get_bin_slice(incr_excess:pandas.core.frame.DataFrame, binstart:float, binstop:float) -> pandas.core.frame.DataFrame
```
Slices the passed dataframe based on the events whose total runoff is
bound by binstart and binstop.

## prep_data_for_convolution
```python
prep_data_for_convolution(dataslice:pandas.core.frame.DataFrame, adj_tempEpsilon:int) -> pandas.core.frame.DataFrame
```
The runoff for each column (event) in the passed dataframe is
calculated from zero to 24 hours for the intervals of length
tempEpsilon*timstep (30 minutes).

## test_shapes
```python
test_shapes(dataslice:pandas.core.frame.DataFrame, col:str, adj_tempEpsilon:int) -> nptyping.Array
```
Calculates the total runoff for each interval, where the interval
width is equal to tempEpsilon times the timestep (30 minutes).

## conv_ts
```python
conv_ts(curve_test_df:pandas.core.frame.DataFrame, convEpsilon:float=150.0, volEpsilon:float=50.0) -> (<class 'dict'>, <class 'list'>)
```
For each event combination, a test statistic is calculated in order
to quantify the similarity between the two temporal distributions.
Note that in this function's code, "c" and "nc" refer to "column"
and "next column", respectively.

## test_stat
```python
test_stat(c_df:pandas.core.frame.DataFrame, nc_df:pandas.core.frame.DataFrame, c:str, nc:str, convEpsilon:float, volEpsilon:float) -> float
```
Calculates a test statistic that quantifies the similarity between
the two curves defined by "c" and "nc" within the passed dataframes.
Note that in this function's code, "c" and "nc" refer to "column"
and "next column", respectively.

## group_curves
```python
group_curves(test_dic:dict, test_values:list, events:list, test_stat_threshold:float=0.0) -> dict
```
If the test statistic for a particular pair of events is greater than
the threshold and neither of the events are already in a group, add
them to a new group. Add all curves that are not a part of a group,
to their own group.

## calc_mean_curves
```python
calc_mean_curves(curve_group:dict, dataslice:pandas.core.frame.DataFrame) -> pandas.core.frame.DataFrame
```
Calculate the mean of the temporal distributions within each group.

## check_upd_curv
```python
check_upd_curv(all_groups:dict, updated_curves:pandas.core.frame.DataFrame, df:pandas.core.frame.DataFrame, convEpsilon:float, volEpsilon:float, test_stat_threshold:float) -> (<class 'dict'>, <class 'pandas.core.frame.DataFrame'>)
```
The temporal distribution for each event within a group used to
calculate a mean temporal distribution is compared to that mean
temporal distribution using the same test statistic used to intially
combine the distributions into groups. If the test statistic for that
distribution is less than the test statistic threshold, the
distribution and its corresponding subgroup are removed from the
overall group used to calculate the mean curve.
The subgroup and remainder of the original group are assigned to new,
separate groups. Once all distributions have been checked against
their mean distributions, the new groups are used to calculated
updated mean distributions.

## extract_list
```python
extract_list(nested_list:list) -> list
```
Extract all of the elements from the sublists within the list and
return the elements as a list.

## map_curve_groups
```python
map_curve_groups(curve_group:dict, curve_group1:dict, ungroup:bool=False) -> dict
```
Map the temporary event keys back to the orignal event IDs to keep a
record of events within each group.

## renumber_dic_keys
```python
renumber_dic_keys(updated_group:dict, group_start_num:int) -> dict
```
Renumber the dictionary keys so that they are ascending.

## final_test_stat
```python
final_test_stat(updated_group:dict, updated_curves:pandas.core.frame.DataFrame, df:pandas.core.frame.DataFrame, convEpsilon:float, volEpsilon:float) -> dict
```
For each group of distributions, the test statistic for each temporal
distribution and corresponding mean temporal distribution (the group
average) is calculated.

## dic_to_list
```python
dic_to_list(dic:dict, get_set:bool=False) -> list
```
Extracts the values from each key within a dictionary and returns the
values as a single list.

## Calc_Group_Weight
```python
Calc_Group_Weight(final_groups:dict, duration_weight:float, display_print:bool=True) -> dict
```
Calculates the weight of each group of curves, such that the sum of
all the weights adds to the duration_weight.

## Rename_Final_Groups
```python
Rename_Final_Groups(curve_weight:dict, dur:int) -> dict
```
Sorts the groups by their weight and then renames the groups so that
the group with the largest weight is designed E0001 and the group with
the next largest weight is designated E0002 (for the 6 hour duration).
The thounsands place is set to 0, 1, 2, 3 for the 6, 12, 24, and 96
hour durations, respectively. A dictionary mapping the original group
names to the new group names is returned.

## dic_key_to_str
```python
dic_key_to_str(orig_dic:dict) -> dict
```
Converts the keys of the passed dictionary to strings.

