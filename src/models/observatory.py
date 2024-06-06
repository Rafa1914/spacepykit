from src.models.observation import Observation
from src.orbit.orbit_determination import get_geodetic_latitude, sidereal_time
import numpy as np

class Observatory():
    _observations = [Observation]
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

    def add_observation(self,observation:Observation):
        observation.sidereal_time = sidereal_time(self._longitude,date=observation.time),
        self._observations.append(observation)
    
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