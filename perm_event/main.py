import netCDF4 as nc
import numpy as np
import pore
import sys 


def open_dataset(filename):
    '''Lee el archivo a evaluar
        Pre: filename: archivo tipo netCDF4.
        Post: atoms_coordinates: un arreglo 3D de la forma (frames, atoms, xyz) 
              de los átomos a evaluar.
    '''
    fn = nc.Dataset(filename)
    atoms_coordinates = fn['coordinates']
    return atoms_coordinates

def pore_traject(atoms_coordinates, pore_radius = 6):
    '''
    Calcula la trayectoria del poro  
       
       Pre: atoms_coordinates: arreglo 3D de la forma (frames, atoms, xyz_coordinates)
            pore_radius: un entero (radio del cilindro)
       Pos: Devuelve un objeto (poro) a partir de los átomos de referencia.
    '''

    # Coordenadas de los átomos de referencia (3)
    a_top = atoms_coordinates[:,0,2]
    a_low = atoms_coordinates[:,2,2]
    a1_xy = atoms_coordinates[:,1,0:2]
    a2_xy = atoms_coordinates[:,2,0:2]

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
    true_list = np.where((condition_1 & condition_2) == True)
    return atom_list[true_list]



filename = 'trajectory.nc'
with nc.Dataset(filename) as dataset:
    atoms_coordinates = dataset['coordinates']

Pore = pore_traject(atoms_coordinates)
atom_list = [] #TODO

compendio_atomos = []
for n_frame, frame in enumerate(atoms_coordinates):
    lista_true_atomos = hacer_comparacion_un_solo_saque(atoms_coordinates, n_frame, atom_list, Pore)
    compendio_atomos.append(lista_true_atomos)



n_eventos = 0

#buscar átomos que pasan por el poro alguna vez: hacer lista atom_pore_list

#iterar por la lista de átomos...
    #iterar por cada frame, desde el 1...
        #evaluar si en el frame n-1 está dentro del cilindro(booleano)
        #evaluar si en el frame n está dentro del cilindro (booleano)

            #si en n-1 no estaba y en el n está dentro del cilindro:
            #evaluar si está más cerca del top o del bottom (ver por donde entró)
            #almacenar tag del átomo de entrada

            #si en n no está en el cilindro:
            #evaluar si está más cerca del top o del bottom (ver por donde salió)
            #almacenar tag del átomo de salida
                #si tag_atom_salida != tag_atom_entrada:
                #n_eventos += 1