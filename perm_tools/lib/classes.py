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
        self.length = abs(self.top - self.low)

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