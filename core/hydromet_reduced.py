from hydromet import*


#---------------------------------------------------------------------------#
def main(EventsTable: dict, durations: list, BCN: str, 
            rand_rate_cap: bool=False, rate: float=None, maxcap: float=None, 
                    minrate: float=None, maxrate: float=None, seed: int=None, 
                                          display_print: bool=True) -> list:
    '''Calculates the reduced excess rainfall for each event within the 
       EventsTable dictionary. 
        
       Parameters
       ----------
       Eventstable: A dictionary containing the incremental excess rainfall 
                    for a suite of randomly generated events.
       durations: A list of event durations, i.e ['H06', 'H12', ...]
       BCN: The boundary condition name as a string, i.e 'D01'
       rand_rate_cap: Bool indicating whether to randomly select the stormwater
                      removal rate and capacity or whether to use the 
                      user-specified values. 
       rate: Stormwater removal rate as a float.
       maxcap: Maximum stormwater system capacity as a float.
       minrate: The minimum stormwater removal rate as a float. This is used
                 when rand_rate_cap is set to 'True'.
       maxrate: The maximum stormwater removal rate as a float. This is used
                when rand_rate_cap is set to 'True'.
       seed: The random number generator seed as an integer.
       display_print: Bool specifying whether to display print statements..
   
       Returns
       -------
       results: A list of dictionaries, which included the incremental reduced
                excess rainfall for each event, the incremental stormwater 
                amount for each event, and the metadata.
    '''  
    RTab = {}
    STab = {}
    SW_variables = {}
    for d in durations:
        dic_dur = EventsTable[d]
        tord = dic_dur['time_idx_ordinate']
        run_dur = dic_dur['run_duration_days']
        tidx = dic_dur['time_idx']
        pluvial_BC_units = dic_dur['pluvial_BC_units']
        if display_print: print('Duration:', d)
        ts = determine_timestep(dic_dur, display_print)
        if rand_rate_cap:
            minrate30 = minrate*(ts*2.0)
            maxrate30 = maxrate*(ts*2.0)
            SW = storm_water_simulator(minrate30, maxrate30, ts, seed, display_print)
            adj_rate = SW[0]
            maxcap = SW[1]
            seed = SW[2]
        else:
            adj_rate = rate*(ts*2.0)
            if display_print: 
                print('Adjusted rate:', adj_rate, 'Adjusted maximum '
                                                        'capacity:', maxcap)
        dic_BCN = {}
        dic_BCN_SW = {}
        dic_reduced = {}
        dic_stormwater = {}
        dic_events = dic_dur['BCName'][BCN]
        for event in dic_events.keys():
            unred = dic_events[event]
            red = reduced_excess(unred, adj_rate, maxcap)
            dic_reduced[event] = red
            dic_stormwater[event] = list(np.array(unred)-np.array(red))
            dic_BCN[BCN] = dic_reduced
            dic_BCN_SW[BCN] = dic_stormwater
        RTab[d] = {'time_idx_ordinate': tord, 'run_duration_days': run_dur, 
                    'time_idx': tidx, 'pluvial_BC_units': pluvial_BC_units, 
                                                          'BCName': dic_BCN}
        STab[d] = {'time_idx_ordinate': tord, 'run_duration_days': run_dur, 
                    'time_idx': tidx, 'pluvial_BC_units': pluvial_BC_units,
                                                       'BCName': dic_BCN_SW}
        SW_variables[d] = {'Rate': adj_rate, 'Capacity': maxcap, 'Seed':seed}
    results = [RTab, STab, SW_variables]
    return results

if __name__== "__main__":
    main()


#---------------------------------------------------------------------------#    