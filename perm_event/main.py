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