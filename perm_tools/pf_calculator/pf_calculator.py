from lib.classes import *
from lib.functions import *
import matplotlib.pyplot as plt


filename = 'perm.nc'

ref_z_1_atom = 0
ref_z_2_atom = 1
ref_xy_1_atom = 0
ref_xy_2_atom = 2

data = open_dataset(filename)
atoms_coordinates = data['coordinates']
total_atoms = len(data.dimensions['atom'])
total_frames = len(data.dimensions['frame'])
top_atom, low_atom = get_z_top_low(atoms_coordinates, ref_z_1_atom, ref_z_2_atom)
first_non_ref_atom = 3
Pore = pore_traject(atoms_coordinates, top_atom, low_atom, ref_xy_1_atom, ref_xy_2_atom)
atom_list = list(range(first_non_ref_atom,total_atoms))

compendio_atomos = atoms_inside_pore(atoms_coordinates, atom_list, Pore)
coord_atoms_in_pore = atoms_coordinates[:,compendio_atomos,:]
dz = coord_atoms_in_pore[1:,:,2] - coord_atoms_in_pore[:-1,:,2]
dz_exp = np.expand_dims(dz, axis=2)
atoms_dz_array = np.concatenate(
    (coord_atoms_in_pore[1:,:,:],dz_exp),
     axis=2)
atoms_dz_array = drop_dz(coord_atoms_in_pore,atoms_dz_array, Pore)
n_array = np.insert(compute_n(atoms_dz_array, Pore),0,0)
n_fragments = 200
n_msd = compute_msd(n_array, n_fragments)
time_axis = np.arange(n_msd.size)
plt.plot(time_axis, n_msd, 'bo')
plt.show()


