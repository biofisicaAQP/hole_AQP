from ...hole_AQP.perm_event import *


filename = 'name.nc'

ref_z_1_atom = 0
ref_z_2_atom = 1
ref_xy_1_atom = 0
ref_xy_2_atom = 2

data = pore_fn.open_dataset(filename)
atoms_coordinates = data['coordinates']
total_atoms = len(data.dimensions['atom'])
total_frames = len(data.dimensions['frame'])
top_atom, low_atom = pore_fn.get_z_top_low(atoms_coordinates, ref_z_1_atom, ref_z_2_atom)
first_non_ref_atom = 3
Pore = pore_fn.pore_traject(atoms_coordinates, top_atom, low_atom, ref_xy_1_atom, ref_xy_2_atom)
atom_list = list(range(first_non_ref_atom,total_atoms))

compendio_atomos = pore_fn.atoms_inside_pore(atoms_coordinates)

