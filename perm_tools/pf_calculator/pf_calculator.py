from lib.classes import *
from lib.functions import *
import matplotlib.pyplot as plt
def pf_calculator(filename, chain_id, ref_z_1_order, ref_z_2_order, ref_xy_1_order, ref_xy_2_order,timestep, pore_radius = 6):
    n_ref_atoms = 3
    first_non_ref_atom = 12
    drop_msd_points = 10
    vol_h2o = 2.99003322259e-23
    ref_z_1_atom = ref_z_1_order + chain_id*n_ref_atoms
    ref_z_2_atom = ref_z_2_order + chain_id*n_ref_atoms
    ref_xy_1_atom = ref_xy_1_order + chain_id*n_ref_atoms
    ref_xy_2_atom = ref_xy_2_order + chain_id*n_ref_atoms

    data = open_dataset(filename)
    atoms_coordinates = data['coordinates']
    total_atoms = len(data.dimensions['atom'])
    total_frames = len(data.dimensions['frame'])
    top_atom, low_atom = get_z_top_low(atoms_coordinates, ref_z_1_atom, ref_z_2_atom)
    
    Pore = pore_traject(atoms_coordinates, top_atom, low_atom, ref_xy_1_atom, ref_xy_2_atom, pore_radius)
    atom_list = list(range(first_non_ref_atom,total_atoms))
    compendio_atomos = atoms_inside_pore(atoms_coordinates, atom_list, Pore)
    coord_atoms_in_pore = atoms_coordinates[:,compendio_atomos,:]
    dz = coord_atoms_in_pore[1:,:,2] - coord_atoms_in_pore[:-1,:,2]
    dz = np.insert(dz, dz.shape[0],np.zeros(shape = dz.shape[1]), axis = 0)
    dz_exp = np.expand_dims(dz, axis=2)
    atoms_dz_array = np.concatenate(
        (coord_atoms_in_pore,dz_exp),
        axis=2)
    atoms_dz_array = drop_dz(atoms_dz_array, Pore)
    n_array = compute_n(atoms_dz_array, Pore)
    fragments = 200
    n_msd = compute_msd(n_array, fragments)
    
    
    time_axis = np.arange(n_msd[drop_msd_points:].size)
    pf = (regresion_lineal(time_axis, n_msd[drop_msd_points:]))*vol_h2o/(2*timestep)
    #print(f'{pf/1e-14:.4f}e-14')
    msd_plot = plt.plot(time_axis, n_msd[drop_msd_points:], 'bo')
    return (pf, msd_plot, compendio_atomos)


