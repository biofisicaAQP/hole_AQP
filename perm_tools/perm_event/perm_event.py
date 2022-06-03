import netCDF4 as nc
import numpy as np
from lib import classes
from lib.functions import *


filename = 'AOX_ref_atoms.nc'
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
n_eventos = 0

compendio_atomos = atoms_inside_pore(atoms_coordinates)

for atom in compendio_atomos:
    entro = False
    for frame in range(1,total_frames):
        if frame % 100 == 0:
            print('Atom:', atom, ' Frame:', frame)
        is_in_cylinder_prev = is_in_cylinder(atoms_coordinates, frame-1,atom,Pore)
        is_in_cylinder_now = is_in_cylinder(atoms_coordinates, frame, atom, Pore)

        if is_in_cylinder_now and is_in_cylinder_prev:
            continue

        if is_in_cylinder_now and not is_in_cylinder_prev:
            ingreso = por_donde_pasa(atoms_coordinates, frame, top_atom, low_atom, atom)
            entro = True
        
        if is_in_cylinder_prev and not is_in_cylinder_now and entro:
            salida = por_donde_pasa(atoms_coordinates, frame, top_atom, low_atom, atom)
            if ingreso != salida:
                n_eventos += 1
            print(n_eventos)
            entro = False
print(n_eventos)





#n_eventos = 0

#buscar átomos que pasan por el poro alguna vez: hacer lista atom_pore_list DONE

#iterar por la lista de átomos... DONE
    #iterar por cada frame, desde el 1... DONE
        #evaluar si en el frame n-1 está dentro del cilindro(booleano) DONE
        #evaluar si en el frame n está dentro del cilindro (booleano) DONE

            #si en n-1 no estaba y en el n está dentro del cilindro: DONE
            #evaluar si está más cerca del top o del bottom (ver por donde entró) DONE
            #almacenar tag del átomo de entrada DONE

            #si en n no está en el cilindro:
            #evaluar si está más cerca del top o del bottom (ver por donde salió)
            #almacenar tag del átomo de salida
                #si tag_atom_salida != tag_atom_entrada:
                #n_eventos += 1