from src.models.observation import Observation

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

    def add_observation(self,observation:Observation):
        self._observations.append(observation)