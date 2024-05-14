from datetime import datetime

class Observation():
    def __init__(self,ra:float,dec:float,time:datetime) -> None:
        '''
        Constructor of observatory's class.

        Attributes:
        ra: Right ascension angle of observation [deg]
        dec: Declination angle of observation [deg]
        time: Time of observation
        '''
        self._ra = ra
        self._dec = dec
        self._time = time