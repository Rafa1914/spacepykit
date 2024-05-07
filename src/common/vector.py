import numpy as np

class Vector:
    '''
    Class that represents a vector in a 3D reference frame.

    Attributes:
        x           Position in x-direction
        y           Position in y-direction
        z           Position in z-direction
        magnitude   Magnitude of the vector
    '''
    def __init__(self,x:float,y:float,z:float) -> None:
        self._x = x
        self._y = y
        self._z = z
        self._magnitude = np.sqrt(x**2+y**2+z**2)
    
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
    
    def __str__(self) -> str:
        return f'x = {str(self.x)}\ny = {str(self.y)}\nz = {str(self.z)}\nMagnitude = {str(self.magnitude)}'