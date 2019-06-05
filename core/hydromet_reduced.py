from hydromet import*


#---------------------------------------------------------------------------#

def main(durations) -> :
    '''TBD
    '''
    ReducedTable = {}
    StormwaterTable = {}
    SW_variables = {}
    for dur in durations:
        dic_dur = EventsTable[dur]
        time_ord = dic_dur['time_idx_ordinate']
        time_idx = dic_dur['time_idx']
        print('Duration:', dur)
        ts = determine_timestep(dic_dur, display_print)
        adj_rate, maximum_capacity, seed = storm_water_simulator(minrate, maxrate, ts, seed, display_print)
        dic_BCN = {}
        dic_BCN_SW = {}
        for name in BCN:
            dic_events = dic_dur['BCName'][name]
            dic_reduced = {}
            dic_stormwater = {}
            for event in dic_events.keys():
                unreduced = dic_events[event]
                reduced = calculate_reduced_excess(unreduced, adj_rate, maximum_capacity)
                dic_reduced[event] = reduced
                dic_stormwater[event] = list(np.array(unreduced)-np.array(reduced))
            dic_BCN[name] = dic_reduced
            dic_BCN_SW[name] = dic_stormwater
        ReducedTable[dur] = {'time_idx_ordinate':time_ord,'time_idx':time_idx, 'BCName': dic_BCN}
        StormwaterTable[dur] = {'time_idx_ordinate':time_ord,'time_idx':time_idx, 'BCName': dic_BCN_SW}
        SW_variables[dur] = {'Rate': adj_rate, 'Capacity': maximum_capacity, 'Seed':seed}
    return results

if __name__== "__main__":
    main()


#---------------------------------------------------------------------------#    