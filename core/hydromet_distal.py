from hydromet import*


#---------------------------------------------------------------------------#
def main(md: dict, weights_dic: dict, durations: list, mainBCN: str, CN: int, 
    arc_data: dict, Project_Area: str,  Pluvial_Model: str, distalBCN: str, 
                outputs_dir: plib, time_idx_ordinate: str, run_dur_dic: dict, 
                        pluvial_BC_units: str, adjust_CN_less24: bool = False, 
                remove_intermediates: bool = True, display_print: bool = True, 
        plot: bool = True, pad_forcing: bool = True, uniform_pad: bool = True, 
                                                    pad_num: int = 2) -> None:
    '''Extracts data from the metadata dictionary, calculates random curve 
       numbers, performs the excess rainfall calculation, groups the events,
       saves the grouped incremental excess rainfall and metadata, and plots
       the results.

       Parameters
       ----------
       md: metadata dictionary which contains the metadata for the 
           precipitation events that were generated in EventsTable.ipynb.
       weights_dic: weights dictionary which contains the weight of each curve
                    group in each duration.
       durations: the event durations as a list, i.e ['H06', 'H12', ...]
       mainBCN: the main domain name as a string. This is the domain used to
                calculate the metadata contained with the metadata dictionary. 
       CN: the curve number for the distal domain as an integer.
       arc_data: a dictionary containing the AMCI and AMCIII values for the
                 specified curve number.
       Project_Area: project area name as a string.
       Pluvial_Model: name of the pluvial model as a string.
       distalBCN: the distal domain name as string. This is the domain whose
                  excess rainfall events are to be calculated using the 
                  provided curve number. 
       outputs_dir: The path for saving the outputs, including intermediate
                    and final results.
       time_idx_ordinate: The ordinate of the time index, i.e. minutes, 
                          hours, days, etc.
       run_dur_dic: dictionary containing the run duration for each event
                    duration.
       pluvial_BC_units: The units of the pluvial boundary condition, i.e. 
                         of the excess rainfall applied to this boundary.
       adjust_CN_less24: Bool specifying whether to adjust the curve number
                         when the storm duration is less than 24 hours.
       remove_intermediates: Bool specifying whether to remove the 
                             intermediate randomized data files once they 
                             have been added to the final metadata file.
       display_print: Bool specifying whether to display print statements.
       plot: Bool specifying whether to display plots.
       pad_forcing: Bool specifying whether to pad the forcing time series. 
       uniform_pad: Bool specifying whether to uniformly pad the forcing 
                    time series or pad it so that it is the same length as
                    the run duration.
       pad_num: The number of zeros to uniformly pad the forcing time series.
       
       Returns
       -------
       None

    '''
    outfiles = []
    for dur in durations:
        idur = int(dur.replace('H', ''))
        if display_print: 
              print('Calculating excess rainfall and grouping the {} hour '
                                    'duration for {}'.format(idur, distalBCN))
        scen = md[dur]['BCName'][mainBCN]
        groups = scen['groups']
        precip = scen['precip']
        metadata = scen['events_metadata']
        eventID = metadata['EventID']
        nevents = len(eventID.keys())
        params = scen['parameters']
        seed = params['seed']
        tempE = params['tempEpsilon']
        tempE2 = params['tempEpsilon2']
        convE = params['convEpsilon']
        volE = params['volEpsilon']
        df_CN = prep_cn_table(CN, arc_data) 
        fitted_cn = find_optimal_curve_beta_dist_S(df_CN)
        fnbase = 'Dur{0}_tempE{1}_convE{2}_volE{3}'.format(idur, tempE, 
                                                                convE, volE)
        fn_CN = "Rand_CN_{0}_{1}_Se{2}.csv".format(distalBCN, fnbase, seed)
        random_cns = RandomizeData(fitted_cn, nevents, outputs_dir, 
                fn_CN, seed = seed, variable = 'CN', display_print = False)
        cum_excess, final_precip, incr_excess = calc_excess_rainfall(eventID, 
                                precip, random_cns, idur, adjust_CN_less24)
        final_curves = calc_mean_curves(groups, incr_excess) 
        fn_Excess = 'Excess_Rainfall_{0}_{1}.csv'.format(distalBCN, fnbase) 
        final_curves.to_csv(outputs_dir/fn_Excess)
        outfiles.append(fn_Excess)
        dic_metadata = {}
        updated_metadata = {}
        dic_metadata['groups'] = dic_key_to_str(groups)
        dic_metadata['precip'] = final_precip.to_dict()
        dic_metadata['cum_excess'] = cum_excess.to_dict()
        dic_metadata['incr_excess'] = incr_excess.to_dict()
        dic_metadata['parameters'] = {'seed':seed, 'tempEpsilon': tempE,
                                      'tempEpsilon2': tempE2, 
                                      'convEpsilon': convE,
                                      'volEpsilon': volE}
        for k in metadata:
              if 'CN' not in k:
                  updated_metadata[k] = metadata[k]
        df = pd.read_csv(outputs_dir/fn_CN, index_col = 'E')
        if remove_intermediates:
              os.remove(outputs_dir/fn_CN)
        new_col = []
        for col in list(df.columns):
              if 'CN' not in col:
                  new_col.append(col+' CN')
              else:
                  new_col.append(col)
        df.columns = new_col       
        for col in df.columns:
              updated_metadata[col] = df[col].to_dict()  
        dic_metadata['events_metadata'] = updated_metadata
        fn_MD = 'Metadata_{0}_{1}.json'.format(distalBCN, fnbase)
        with open(outputs_dir/fn_MD, 'w') as f:
              json.dump(dic_metadata, f)
        outfiles.append(fn_MD)
        if plot:
              plot_rainfall_and_excess(final_precip, cum_excess, idur)
              y_max = final_curves.max().max()
              plot_grouped_curves(final_curves, y_max) 
    outfiles = extract_list(outfiles)
    excess_dic = combine_distal_results(outfiles, outputs_dir, 'Excess',
                distalBCN, time_idx_ordinate, pluvial_BC_units, run_dur_dic, 
                                                        remove_intermediates)
    if pad_forcing: 
        excess_dic = pad_pluvial_forcing(excess_dic, uniform_pad, pad_num, 
                                                            verbose = False)
    fn_final = '{0}_{1}_{2}'.format(Project_Area, Pluvial_Model, distalBCN)
    with open(outputs_dir/'{0}.json'.format(fn_final),'w') as f:
        json.dump(excess_dic, f)     
    metadata = combine_distal_results(outfiles, outputs_dir, 'Metadata', 
                            distalBCN, remove_ind_dur = remove_intermediates)
    with open(outputs_dir/'{0}_Metadata.json'.format(fn_final),'w') as f:
        json.dump(metadata, f)   
    if plot: 
        plot_amount_vs_weight(weights_dic, excess_dic, mainBCN, distalBCN)                  
    return 
if __name__== "__main__":
    main()


#---------------------------------------------------------------------------#    