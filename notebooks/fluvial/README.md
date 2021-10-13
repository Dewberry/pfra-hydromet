# Description
These jupyter notebooks ingest [HEC-SSP](https://www.hec.usace.army.mil/software/hec-ssp/) .rpt files containing flow frequency data for a specific USGS Stream Gage calculated at a range of confidence limits and return a series of events, statified by the annual exceedance probability, with discharge based on the mean flow frequency curve. This approach relies on:
   1. [Bulletin 17C](https://pubs.usgs.gov/tm/04/b05/tm4b5.pdf) flow frequency analysis
   2. Mean flow frequency curve
   3. Stratified sampling
   


__NOTE__: [SSP_to_Mean_Curve](SSP_to_Mean_Curve.ipynb) is the primary notebook for calculating the mean flow frequency curve and can be called by [*Papermill*](https://pypi.org/project/papermill/). Papermill is designed to act as a *manager* to maintain consistency in computation, and ensure cells are executed in order. *Manager* notebooks are designated with the `PM-` prefix. Executed notebooks should be saved as documentation of the inputs, outputs, and results for a given project location.

---

## Contents
- [__PM_Sampler_Ops__](PM_Sampler_Ops.ipynb):  Manager notebook that executes `SSP_to_Mean_Curve`, `Stratified_Sampler`, and `Make_Production_Run_List`.

- [__SSP_to_Mean_Curve__](SSP_to_Mean_Curve.ipynb): Calculates the mean flow frequency curve using [Bulletin 17C](https://pubs.usgs.gov/tm/04/b05/tm4b5.pdf) confidence limits calculated in [HEC-SSP](https://www.hec.usace.army.mil/software/hec-ssp/).

- [__Stratified_Sampler__](Stratified_Sampler.ipynb): Calculates the weight of a specified number of annual exceedance probabilities/recurrence intervals uniformly selected between the minimum and maximum value within log space.

- [__Make_Production_Run_List__](Make_Production_Run_List.ipynb): Calculates the discharge for each annual exceedance probability (AEP) within the weights table using the mean flow frequency curve.

---

## Workflow

1. Run [HEC-SSP](https://www.hec.usace.army.mil/software/hec-ssp/) version 2.1 or 2.2 for a range of confidence limits and annual exceedance probabilities to calculate mean flow frequency curves for each confidence limit.

    ```
      Inputs:
        1. Annual peak flow imported using the "Data Importer" in HEC-SSP. Ensure "USGS Website", "Flow", and "Annual Peak Data" are selected; click "Get USGS Station ID's by State"; select a gage; and click "Import to Study DSS File".
        2. Create a new Bulletin 17 analysis, being sure to name the analysis using the following convention: "USGSStationID_UpperCL" where "UpperCL" is the upper confidence limit that will be specified in the "Confidence Limits" field of the "Bulletin 17 Editor". Enter the upper confidence limit and the lower limit is populated. 
        3. Specify the range of annual exceedance probabilities by clicking on the "Options" table and entering a list of annual exceedance probabilities (units in percent). For example, use values between 1E-10 and 99.0.
        4. Run the Bulletin 17 analysis for a range of confidence limits. For example, use upper confidence limits of 0.60, 0.70, 0.80, 0.90, 0.95, 0.99, 0.995, and 0.999.
        
      Outputs:
        1. .rpt files containing the flow frequency data for the selected USGS Stream Gage calculated at a range of confidence limits and annual exceedance probabilities.
    ```
    
    
2. Run [PM_Sampler_Ops](PM_Sampler_Ops.ipynb) which executes [SSP_to_Mean_Curve](SSP_to_Mean_Curve.ipynb) to calculate the mean flow frequency curve, [Stratified_Sampler](Stratified_Sampler.ipynb) to calculate the weight of a specified number of events, and [Make_Production_Run_List](Make_Production_Run_List.ipynb) to interpolate the discharge for each event from the mean flow frequency curve.

    ```
      Inputs:
      1. The .rpt files from step 1.
      2. Options such as the upper confidence limit to restrict the confidence limit range and the amount to adjust the discharge for each annual exceedance probability if the discharge does not increase with decreasing probability.
      
      Outputs:
      1. The mean flow frequency curve table (annual exceedance probability verses discharge).
      2. A table of annual exceedance probabilities and their corresponding weight.
      3. A table with the annual exceedance probability, discharge, and weight.
      4. Excel workbook with a sheet for the mean flow frequency curve table, the annual exceedance probability (AEP) vs weight table, and the discharge verses weight table.
      5. A copy of each of the notebooks executed by PM_Sampler_Ops.ipynb.
      
    ```