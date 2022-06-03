print('__file__={0:<35} | __name__={1:<20} | __package__={2:<20}'.format(__file__,__name__,str(__package__)))
from lib.classes import *
from lib.functions import *


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
print(compendio_atomos)

