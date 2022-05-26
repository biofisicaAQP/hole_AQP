import netCDF4 as nc
import numpy as np
import pore
import sys 


def open_dataset(filename):
    '''Lee el archivo a evaluar
        Pre: filename: archivo tipo netCDF4.
        Post: fn: data del tipo netCDF.
    '''
    fn = nc.Dataset(filename)
    
    return fn
def get_z_top_low(atoms_coordinates, atom_1, atom_2):
    avg_atom1 = atoms_coordinates[:,atom_1,2].mean()
    avg_atom2 = atoms_coordinates[:,atom_2,2].mean()
    if avg_atom1 > avg_atom2:
        top_atom = atom_1
        low_atom = atom_2        
    else:
        top_atom = atom_2
        low_atom = atom_1
    return(top_atom, low_atom)

def pore_traject(atoms_coordinates,top_atom, low_atom, ref_xy_1_atom, ref_xy_2_atom, pore_radius = 6):
    '''
    Calcula la trayectoria del poro  
       
       Pre: atoms_coordinates: arreglo 3D de la forma (frames, atoms, xyz_coordinates)
            pore_radius: un entero (radio del cilindro)
       Pos: Devuelve un objeto (poro) a partir de los átomos de referencia.
    '''

    # Coordenadas de los átomos de referencia (3)
    a_top = atoms_coordinates[:,top_atom,2]
    a_low = atoms_coordinates[:,low_atom,2]
    a1_xy = atoms_coordinates[:,ref_xy_1_atom,0:2]
    a2_xy = atoms_coordinates[:,ref_xy_2_atom,0:2]

    # Se llama a la clase pore y se determina la trayectoria del 
    # poro a partir de las coordenadas de los átomos de referencia y el radio.

    Pore = pore.Pore_cylinder_traj(a_top,a_low,a1_xy,a2_xy,pore_radius)      #cosa a mejorar: poder elegir 
                                                                             #los átomos de referencia.
                                                                        
    return Pore

def hacer_comparacion_un_solo_saque(atoms_coordinates, frame, atom_list,Pore):
    '''Devuelve una lista de átomos que se encuentran dentro del cilindro
    en un frame determinado
    Pre:
    atoms_coordinates: array 3D de la forma [frames, atoms, coordenadas_xyz]
    frame: un entero
    atom_list: una lista de enteros
    Pore: un objeto Pore_cylinder_traj'''
    pore_cylinder = Pore[frame]
    atom_coord = atoms_coordinates[frame,atom_list,:]
    distance_to_center = (((atom_coord[:,0:2] - pore_cylinder.xy_center)**2).sum(axis = 1))**(1/2)
    condition_1 = distance_to_center < pore_cylinder.radius
    condition_2 = (pore_cylinder.low < atom_coord[:,2]) & (atom_coord[:,2] < pore_cylinder.top)
    true_list = np.where((condition_1 & condition_2) == True)[0].tolist()
    filtered_atom_list = [atom_list[i] for i in true_list]
    return filtered_atom_list

def is_in_cylinder(atoms_coordinates, frame, atom, Pore):
    pore_cylinder = Pore[frame]
    atom_coord = atoms_coordinates[frame,atom,:]
    distance_to_center = (((atom_coord[0:2] - pore_cylinder.xy_center)**2).sum(axis = 0))**(1/2)
    condition_1 = distance_to_center < pore_cylinder.radius
    condition_2 = (pore_cylinder.low < atom_coord[2]) & (atom_coord[2] < pore_cylinder.top)
    return (condition_1 & condition_2)

def dist_interatom(atom_coordinates, frame, atom_1 ,atom_2):
    coord_atom_1 = atom_coordinates[frame, atom_1,:]
    coord_atom_2 = atom_coordinates[frame, atom_2,:]
    distancia = np.sqrt((coord_atom_1[0]-coord_atom_2[0])**2+
                            (coord_atom_1[1]-coord_atom_2[1])**2+
                            (coord_atom_1[2]-coord_atom_2[2])**2) 
    return distancia

def por_donde_pasa(atoms_coordinates, frame, top_atom, low_atom, atom):
    dist_to_low = dist_interatom(atoms_coordinates, frame, low_atom, atom)
    dist_to_top = dist_interatom(atoms_coordinates, frame, top_atom, atom)
    if dist_to_low < dist_to_top:
        return low_atom
    else:
        return top_atom


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

compendio_atomos = []
for n_frame, frame in enumerate(atoms_coordinates):
    lista_true_atomos = hacer_comparacion_un_solo_saque(atoms_coordinates, n_frame, atom_list, Pore)
    compendio_atomos.extend(lista_true_atomos)
compendio_atomos = list(set(compendio_atomos))

for atom in compendio_atomos:
    entro = False
    for frame in range(1,total_frames):
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