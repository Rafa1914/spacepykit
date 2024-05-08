import numpy as np

class Vector:
    '''
    Class that represents a vector in a 3D reference frame.

    Attributes:
        x:           Position in x-direction\n
        y:           Position in y-direction\n
        z:           Position in z-direction\n
        theta_x:     Director angle with respecto to x-axis [rad]\n
        theta_y:     Director angle with respecto to y-axis [rad]\n
        theta_z:     Director angle with respecto to z-axis [rad]\n
        magnitude:   Magnitude of the vector\n
    '''
    def __init__(self,x:float,y:float,z:float) -> None:
        self._x = x
        self._y = y
        self._z = z
        self._magnitude = np.sqrt(x**2+y**2+z**2)
        self._theta_x = np.arccos(x/self._magnitude)
        self._theta_y = np.arccos(y/self._magnitude)
        self._theta_z = np.arccos(z/self._magnitude)
    
    @property
    def x(self) -> float:
        return self._x
    
    @property
    def y(self) -> float:
        return self._y
    
    @property
    def z(self) -> float:
        return self._z
    
    @property
    def magnitude(self) -> float:
        return self._magnitude
    
    @property
    def theta_x(self) -> float:
        return self._theta_x
    
    @property
    def theta_y(self) -> float:
        return self._theta_y
    
    @property
    def theta_z(self) -> float:
        return self._theta_z
    
    def __str__(self) -> str:
        return f'x = {str(self.x)}\ny = {str(self.y)}\nz = {str(self.z)}\nMagnitude = {str(self.magnitude)}\nTheta_x = {np.rad2deg(self.theta_x)}\nTheta_y = {np.rad2deg(self.theta_y)}\nTheta_z = {np.rad2deg(self.theta_z)}'      
    
    def __add__(self,other):
        x = self.x+other.x
        y = self.y+other.y
        z = self.z+other.z
        return Vector(x,y,z)
    
    def __sub__(self,other):
        x = self.x-other.x
        y = self.y-other.y
        z = self.z-other.z
        return Vector(x,y,z)

def dot_product(a:Vector,b:Vector) -> float:
    '''
    Calculates de dot product between both vectors.
    '''
    return a.x*b.x + a.y*b.y + a.z*b.z

def cross_product(a:Vector,b:Vector) -> Vector:
    '''
    Calculates the cross product betwenn both vectors.
    '''
    x = (a.y*b.z-a.z*b.y)
    y = (a.z*b.x-a.x*b.z)
    z = (a.x*b.y-a.y*b.x)
    return Vector(x,y,z)

def find_projection(a:Vector,b:Vector) -> float:
    '''
    Calculate the projection of vector a over vector b.
    '''
    return dot_product(a,b)/b.magnitude

def find_angle_between(a:Vector,b:Vector) -> float:
    '''
    Calculates the angle between vectors a and b.

    Returns:
    theta: Angle between the vectors [rad].
    '''
    theta = np.arccos(dot_product(a,b)/(a.magnitude*b.magnitude))
    return theta