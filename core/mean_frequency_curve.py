import numpy as np
import pandas as pd
import pathlib as pl
from IPython.core.display import display
from meanffc import (list_ssp_files, make_ssp_table, monotonic_test, zvar, binq, interp_aep, zvar_inv, mean_aep, 
                     ffc_summary, interp_q, make_directories, plot_ssp_transform, plot_ssp_interp, 
                     plot_ssp_interptrans, plot_ssp_meanmed, plot_ssp_meanmedffc, plotly_ssp_meanffc)


def main(gage_id: str, inputs_path: str = '', outputs_dir: str = '', data_type: str = 'Q', version: str = '2_2', 
         max_cl: float = 0.99, adj_flows: float = 1.0, round_decimals: int = 1, verbose: bool = True, 
         display_plots: bool = True, extrapolate: bool = True, exclude_tails: bool = True) -> pd.DataFrame:
    """Calculates the mean frequency curve for a standard set of annual exceedance probabilities given the frequency 
       curves calculated at a range of confidence limits. The inputs_path is the path to either a directory containing 
       the .rpt files from [HEC-SSP](https://www.hec.usace.army.mil/software/hec-ssp/) or the path to a csv where the 
       index is the AEP, the columns are confidence limits, and the values are discharge or precipitation depth. By 
       default the mean curve table is returned but not saved, however, if the user specifies the outputs_dir, the table 
       will be saved to a csv. The data_type by default is "Q" for discharge, however, the user may specify a data_type 
       of "P" for precipitation which will update the plot and column labels.
    """
    standard_aep = [0.9, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 5E-04, 2E-04, 1E-04, 5E-05, 2E-05, 1E-05,
                    5E-06, 2E-06, 1E-06]    

    inputs_path = pl.Path(inputs_path)
    if inputs_path.suffix == '':
        if verbose:
            print('Specified inputs_path is to a directory, identifying .rpt files:')
        ssp_results = list_ssp_files(inputs_path, max_cl, verbose)     
        df = make_ssp_table(ssp_results, version) 
    elif inputs_path.suffix == '.csv':
        df = pd.read_csv(inputs_path, index_col = 0)
        assert df.index.name == 'AEP', ('The first column of the specified csv must be titled "AEP" and contain annual '
                                       'exceedance probabilites.')
        if verbose:
            print('Specified inputs_path is to a csv, loaded table.\n')
    else:
        raise NameError('The specified inputs_path is not a directory of .rpt files or a path to a csv file as '
                        'expected.')
    if data_type not in ['Q', 'P']:
        raise ValueError('The specified data_type may only be "Q" or "P" corresponding to discharge or precipitation, '
                         'respectively')               
    df = monotonic_test(df, adj_flows, verbose)
    cl = list(map(float, df.columns)) 
    clz = zvar(cl) 
    aep = df.index 
    aepz = zvar(aep) 
    data = np.log10(df) 
    q = binq(data)
    res = interp_aep(data, q, clz, aepz, extrapolate)    
    restrans = zvar_inv(res, cl)
    aepm = mean_aep(restrans, exclude_tails) 
    aepmz = zvar(aepm)  
    standard_aepz = zvar(standard_aep) 
    table = ffc_summary(standard_aep, standard_aepz) 
    table['Q_Mean_cfs'] = np.round(10**table['AEPz'].apply(interp_q(aepmz, q)), round_decimals)
    table['Q_Median_cfs'] = np.round(10**table['AEPz'].apply(interp_q(aepz, np.array(data['0.5']))), round_decimals)
    mean_curve_table = table.copy().drop(columns=['AEPz'])
    if display_plots:    
        plot_ssp_meanmedffc(table, gage_id, data_type) 
    if data_type == 'P':
        mean_curve_table.rename(columns = {'Q_Mean_cfs': 'P_Mean_in', 'Q_Median_cfs': 'P_Median_in'}, inplace = True)
    if verbose:
        display(mean_curve_table.head(2))    
    if outputs_dir!='':
        make_directories([outputs_dir], verbose)
        mean_curve_table.to_csv(pl.Path(outputs_dir)/f'MeanCurve_{gage_id}.csv')     
    return mean_curve_table

if __name__== "__main__":
    main()