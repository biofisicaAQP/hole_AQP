from matplotlib.pyplot import axis
import netCDF4 as nc
import numpy as np
from lib import classes
from sklearn import linear_model



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
    Pore = classes.Pore_cylinder_traj(a_top,a_low,a1_xy,a2_xy,pore_radius)                                                                        
    return Pore


def regresion_lineal(columna_x, columna_y):
    ajuste = linear_model.LinearRegression()
    ajuste.fit(columna_x.reshape((-1,1)),columna_y)
    return ajuste.coef_[0]


def compute_msd(n_array, fragments):
    splitted_array = np.array(np.split(n_array, fragments))
    split_to_zero = (splitted_array.T-splitted_array[:,0]).T
    msd = np.mean(split_to_zero**2, axis = 0)
    return msd


def compute_n(atoms_dz_array, pore):
    dn_array = (atoms_dz_array[:,:,3].sum(axis = 1))/pore.length
    n_array = np.cumsum(dn_array)
    return n_array


def drop_dz(atoms_dz_array, pore):
    copy_dz_array = atoms_dz_array.copy()
    distance_to_center = (((atoms_dz_array[:,:,0:2] - np.expand_dims(pore.xy_center, axis=1))**2).sum(axis = 2))**(1/2)
    condition_1 = ((pore.low >= atoms_dz_array[:,:,2].T) | (atoms_dz_array[:,:,2].T >= pore.top)).T
    condition_2 = (distance_to_center >= pore.radius)
    copy_dz_array[1:,:,3][((condition_1[1:] | condition_2[1:]) & (condition_1[:-1] | condition_2[:-1]))] = 0
    return copy_dz_array


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
    condition_2 = (pore_cylinder.low.min() < atom_coord[:,2]) & (atom_coord[:,2] < pore_cylinder.top.max())
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



def atoms_inside_pore(atoms_coordinates, atom_list, Pore):
    compendio_atomos = []
    for n_frame, frame in enumerate(atoms_coordinates):
        if n_frame % 500 == 0:
            print(n_frame)
        lista_true_atomos = hacer_comparacion_un_solo_saque(atoms_coordinates, n_frame, atom_list, Pore)
        compendio_atomos.extend(lista_true_atomos)
    compendio_atomos = list(set(compendio_atomos))
    return compendio_atomos