from lib.classes import *
from lib.functions import *
import matplotlib.pyplot as plt


def histogram_calculator(filename, chain_id,ref_z_1_order, ref_z_2_order, ref_xy_1_order, ref_xy_2_order,bins, pore_radius, pedazos=5):
    n_ref_atoms = 3
    first_non_ref_atom = 12
    ref_z_1_atom = ref_z_1_order + chain_id*n_ref_atoms
    ref_z_2_atom = ref_z_2_order + chain_id*n_ref_atoms
    ref_xy_1_atom = ref_xy_1_order + chain_id*n_ref_atoms
    ref_xy_2_atom = ref_xy_2_order + chain_id*n_ref_atoms
    data = open_dataset(filename)
    atoms_coordinates = data['coordinates']
    top_atom, low_atom = get_z_top_low(atoms_coordinates, ref_z_1_atom, ref_z_2_atom)
    Pore = pore_traject(atoms_coordinates, top_atom, low_atom, ref_xy_1_atom, ref_xy_2_atom, pore_radius)
    fragments_marks = [i for i in range(len(atoms_coordinates)) if i%(len(atoms_coordinates)/pedazos) == 0]
    array_list = []
    for i in range(1,len(fragments_marks)):
        start = fragments_marks[i-1]
        stop = fragments_marks[i]
        filtered_array = filter_by_radius(atoms_coordinates[start:stop,first_non_ref_atom:,:],Pore[start:stop])
        array_list.append(filtered_array)
        start = fragments_marks[-1]
    filtered_array = filter_by_radius(atoms_coordinates[start:,first_non_ref_atom:,:],Pore[start:])
    array_list.append(filtered_array)
    filtered_array = np.concatenate(tuple(array_list))
    arr_hist,edges = compute_histogram(filtered_array,bins)
    return(arr_hist,edges,filtered_array.shape[0])

    
def pmf_calculator(arr_hist, edges, array_shape, unit='kJ/mol', temp=298):
    density,edges = compute_density(arr_hist,edges, array_shape)
    N_avog = 6.02e23
    kb = {'kcal':3.297623483e-27, 'kJ':1.380649e-26, 'kJ/mol':N_avog * 1.380649e-26, 'kcal/mol':N_avog * 3.297623483e-27}
    pmf = -kb[unit]*temp*np.log10(density)
    return (pmf,edges)
