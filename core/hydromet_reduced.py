from hydromet import*


#---------------------------------------------------------------------------#
def main(EventsTable: dict, durations: list, BCN: list, 
            rand_rate_cap: bool=False, rate: float=None, maxcap: float=None, 
                    minrate: float=None, maxrate: float=None, seed: int=None, 
                                        display_print: bool=True) -> list:
    '''Calculates the reduced excess rainfall for each event within the 
       EventsTable dictionary. 
    '''
    RTab = {}
    STab = {}
    SW_variables = {}
    for d in durations:
        dic_dur = EventsTable[d]
        tord = dic_dur['time_idx_ordinate']
        tidx = dic_dur['time_idx']
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
                print('Rate:', adj_rate, 'Maximum Capacity:', maxcap)
        dic_BCN = {}
        dic_BCN_SW = {}
        for name in BCN:
            dic_events = dic_dur['BCName'][name]
            dic_reduced = {}
            dic_stormwater = {}
            for event in dic_events.keys():
                unred = dic_events[event]
                red = reduced_excess(unred, adj_rate, maxcap)
                dic_reduced[event] = red
                dic_stormwater[event] = list(np.array(unred)-np.array(red))
            dic_BCN[name] = dic_reduced
            dic_BCN_SW[name] = dic_stormwater
        RTab[d] = {'time_idx_ordinate':tord,'time_idx':tidx,'BCName':dic_BCN}
        STab[d]={'time_idx_ordinate':tord,'time_idx':tidx,'BCName':dic_BCN_SW}
        SW_variables[d] = {'Rate': adj_rate, 'Capacity': maxcap, 'Seed':seed}
    results = [RTab, STab, SW_variables]
    return results

if __name__== "__main__":
    main()


#---------------------------------------------------------------------------#    