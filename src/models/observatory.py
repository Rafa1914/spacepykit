from src.orbit.orbit_determination import get_geodetic_latitude, sidereal_time
import numpy as np

class Observatory():
    def __init__(self,longitude:float,latitude:float,altitude:float) -> None:
        '''
        Constructor of observation's class.

        Attributes:
        longitude: Site location's longitude [deg]
        latitude: Site location's latitude [deg]
        altitude: Site location's altitude [km]
        '''
        self._longitude = longitude
        self._latitude = latitude
        self._altitude = altitude
        self._latitude_geodetic = np.rad2deg(get_geodetic_latitude(latitude))

    
    @property
    def longitude(self) -> float:
        return self._longitude
    @property
    def latitude(self) -> float:
        return self._latitude
    @property
    def altitude(self) -> float:
        return self._altitude
    @property
    def latitude_geodetic(self) -> float:
        return self._latitude_geodetic