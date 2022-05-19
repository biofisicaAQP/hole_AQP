'''Envoltorio para acceder a las coordenadas de un cilindro (poro) definido por
   3 átomos de referencia con su respectiva trayectoria.
'''
class Pore_cylinder:
    '''Se accede a las coordenadas del poro.
       
      Pre: top: (float) coordenada en z del átomo superior.
           low: (float) coordenada en z del átomo inferior.
           xy_1, xy_2: (lista o arreglos) puntos en xy
               diametralmente opuestos en el cilindro.
           radius: (float) radio del cilindro
      Pos: arreglo con las coordenadas del cilindro (poro)
    '''
    def __init__(self, top, low, xy_1, xy_2, radius):
        self.top =  top
        self.low = low
        self.xy_1 = xy_1
        self.xy_2 = xy_2
        self.xy_center = (self.xy_1+self.xy_2)/2
        self.radius = radius

class Pore_cylinder_traj(Pore_cylinder):
    '''Se accede a las coordenadas del poro con múltiples frames
       Pre: arrays de longitudes iguales (frame).
       Pos: se obtiene la trayectoria del poro.
      
    '''
    def __getitem__(self,frame):
        return Pore_cylinder(
            self.top[frame],
            self.low[frame],
            self.xy_1[frame],
            self.xy_2[frame],
            self.radius
        )

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