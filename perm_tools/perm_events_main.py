import perm_event.perm_event as perm_ev_calc
import pandas as pd
filename = 'perm.nc'
chains = [0,1,2,3]
chains_letter = ['A','B','C','D']
(ref_z_1, ref_z_2, ref_xy_1, ref_xy_2) = (0,1,1,2)
list_events_total = []
for chain_id in chains:
    (list_events_chain, n_eventos) = perm_ev_calc.perm_event_calculator(filename, chain_id, ref_z_1, ref_z_2, ref_xy_1, ref_xy_2, chains_letter)
    print(f'Chain {chains_letter[chain_id]} perm events: {n_eventos}')
    list_events_total += list_events_chain

event_df = pd.DataFrame(list_events_total)
event_df.to_csv('event_df.csv')