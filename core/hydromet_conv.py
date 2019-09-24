from hydromet import*


#---------------------------------------------------------------------------#
def main(binData: list, incr_excess: pd.DataFrame, tempE: float, 
                                convE: float, volE: float, tsthresh: float, 
                                        display_print: bool=True) -> dict:
    '''Function for grouping incremental excess rainfall curves using a novel
       test statistic that quantifies the incremental and cumulative 
       volumentric differences between two curves. The mean of each group of
       curves is calculated and used in place of the original set of curves
       in order to improve modeling efficiency by reducing redundancy.

       Parameters
       ----------
       binData: a list of excess rainfall amounts which represent the bounds 
                of the binned excess rainfall events.
       incr_excess: a dataframe containing the incremental excess rainfall 
                    for a suite of randomly generated events.
       tempE: The number of hours over which to resample the excess rainfall 
              event timestep as a float. 
       convE:The maximum allowable percent difference in excess rainfall 
             between two curves at any timestep as a float.
       volE: The maximum allowable precent difference in the total excess 
             rainfall between two curves as a float
       tsthresh: The convolution test statistic threshold as a float. If a 
                 test statistic is above this threshold then two curves are
                 considered to be quantitatively similiar and are grouped.
       display_print: Bool specifying whether to display print statements.
   
       Returns
       -------
       results: A dictionary containing the final and mid-bin incremental
                excess events as dataframes, the resampled mid-bin incremental
                excess events as a dataframe, and the final and mid-bin group 
                IDs as dictionaries. 

    '''  
    penult_curves = pd.DataFrame()
    penult_groups = {}
    penult_tests = {} 
    group_start_num = 0
    adj_tempE = convert_tempEpsilon(tempE, incr_excess)
    for i, b in enumerate(binData): 
        start = time.time()
        if b[1] == binData[0][1]:   
            binstart, binstop = 0, binData[1][1]
        if b[1] == binData[-1][1]: 
            binstart, binstop = b[1], b[1]+100
        else:
            binstart, binstop = binData[i][1], binData[i+1][1]   
        # Get data given the bin threshold    
        dataslice = get_bin_slice(incr_excess, binstart, binstop)
        idx0 = list(dataslice.sum()[dataslice.sum()==0].index)
        idx = list(set(list(dataslice.columns))-set(idx0))
        n_zero =  len(idx0)
        n_nonzero = len(idx)
        if n_zero >= 1:
            curve = dataslice[idx0[0]].copy()
            curve.name = group_start_num
            penult_curves = pd.concat([penult_curves, curve], axis=1)
            penult_groups[group_start_num] = list(dataslice[idx0].columns)
            penult_tests[group_start_num] = [1.0]
            group_start_num += 1
        if n_nonzero > 1:
            # Prep data for convolution test 
            curv_df = prep_data_for_convolution(dataslice[idx], adj_tempE)
            # Perform the convolution to group curves:
            test_dic, test_values = conv_ts(curv_df, convE, volE) 
            events = list(curv_df.columns)                                                                         
            n_init = len(events) 
            all_groups = group_curves(test_dic, test_values, events, tsthresh)                                 
            upd_curv = calc_mean_curves(all_groups, curv_df)                                                   
            updated_group, upd_curv = check_upd_curv(all_groups, upd_curv, curv_df, convE, volE, tsthresh)
            n_fin =  len(upd_curv.columns)
            if display_print: print('Final Nonzero Groups in Count 0: {0};'
                                    ''.format(n_fin), 'Max Test Stat: {0}'
                                                ''.format(max(test_values)))   
            count = 1
            # Repeat the convolution with the grouped curves:
            while count <=10 and max(test_values) >= tsthresh and 1 < n_fin < n_init:
                test_dic, test_values = conv_ts(upd_curv, convE, volE)
                events = list(upd_curv.columns)
                n_init = len(events)
                curve_group = group_curves(test_dic, test_values, events, tsthresh)
                all_groups = map_curve_groups(updated_group, curve_group)
                upd_curv = calc_mean_curves(all_groups, curv_df)
                updated_group, upd_curv = check_upd_curv(all_groups, upd_curv, curv_df, convE, volE, tsthresh)
                n_fin = len(upd_curv.columns)
                if display_print: print('Final Nonzero Groups in Count {0}: '
                            '{1};'.format(count, n_fin), 'Max Test Stat: {0}'
                                                ''.format(max(test_values)))
                count+=1
            # Save the results:
            group_stop_num = group_start_num + upd_curv.columns.values.max()+1
            reordered_curves = upd_curv.reindex(sorted(upd_curv.columns), axis=1)
            reordered_curves.columns = np.arange(group_start_num, group_stop_num)
            reordered_group = renumber_dic_keys(updated_group, group_start_num)
            penult_curves = pd.concat([penult_curves, calc_mean_curves(reordered_group, dataslice[idx])], axis=1)
            penult_groups.update(reordered_group)
            penult_tests.update(final_test_stat(reordered_group, reordered_curves, curv_df, convE, volE))
            group_start_num = group_stop_num
            # Copy one of the bins for plotting/QC purposes
            if i==int(len(binData)/2.0):
                midbin_curve_df = curv_df
                midbin_group = reordered_group
                midbin_curves = reordered_curves
        if n_nonzero == 1:
            curve = dataslice[idx].copy()
            curve.columns = [group_start_num]
            penult_curves = pd.concat([penult_curves, curve], axis=1)
            penult_groups[group_start_num] = [dataslice[idx].columns[0]]
            penult_tests[group_start_num] = [1.0]
            group_start_num += 1
        if display_print: print('Processed Bin {} with {} nonzero curves in '
                            '{} Minutes'.format(i, dataslice[idx].shape[1], 
                                            round(time.time()-start)/60, 3))
    results = {'penult_curves': penult_curves,
               'penult_groups': penult_groups, 
               'penult_tests': penult_tests,
               'midbin_curve_df': midbin_curve_df,
               'midbin_group': midbin_group,
               'midbin_curves': midbin_curves}
    return results

if __name__== "__main__":
    main()


#---------------------------------------------------------------------------#    