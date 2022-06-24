from tkinter import N
import pandas as pd
from lib.functions import *

def perm_event_calculator(filename, chain_id, ref_z_1_order, ref_z_2_order, ref_xy_1_order, ref_xy_2_order, chains_letter, pore_radius = 6):
    n_ref_atoms = 3
    first_non_ref_atom = 12
    ref_z_1_atom = ref_z_1_order + chain_id*n_ref_atoms
    ref_z_2_atom = ref_z_2_order + chain_id*n_ref_atoms
    ref_xy_1_atom = ref_xy_1_order + chain_id*n_ref_atoms
    ref_xy_2_atom = ref_xy_2_order + chain_id*n_ref_atoms

    data = open_dataset(filename)
    atoms_coordinates = data['coordinates']
    total_atoms = len(data.dimensions['atom'])
    total_frames = len(data.dimensions['frame'])
    top_atom, low_atom = get_z_top_low(atoms_coordinates, ref_z_1_atom, ref_z_2_atom)
    first_non_ref_atom = 12
    Pore = pore_traject(atoms_coordinates, top_atom, low_atom, ref_xy_1_atom, ref_xy_2_atom, pore_radius)
    atom_list = list(range(first_non_ref_atom,total_atoms))
    n_eventos = 0

    compendio_atomos = atoms_inside_pore(atoms_coordinates, atom_list, Pore)
    list_events = []
    for atom in compendio_atomos:
        entro = False
        for frame in range(1,total_frames):
            if frame % 10000 == 0:
                print('Atom:', atom, ' Frame:', frame)
            is_in_cylinder_prev = is_in_cylinder(atoms_coordinates, frame-1,atom,Pore)
            is_in_cylinder_now = is_in_cylinder(atoms_coordinates, frame, atom, Pore)

            if is_in_cylinder_now and is_in_cylinder_prev:
                continue

            if is_in_cylinder_now and not is_in_cylinder_prev:
                ingreso = por_donde_pasa(atoms_coordinates, frame, top_atom, low_atom, atom)
                entro = True
                frame_entrada = frame
            
            if is_in_cylinder_prev and not is_in_cylinder_now and entro:
                salida = por_donde_pasa(atoms_coordinates, frame, top_atom, low_atom, atom)
                if ingreso != salida:
                    n_eventos += 1
                    frame_salida = frame
                    event_dict = {'Chain':chains_letter[chain_id],'Atom': atom, 'Start':frame_entrada, 'End': frame_salida}
                    list_events.append(event_dict)
                    print(n_eventos)
                entro = False
    print(n_eventos)
    return(list_events, n_eventos)
