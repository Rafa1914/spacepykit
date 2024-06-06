from datetime import datetime


class Observation:
    def __init__(self, ra: float, dec: float, time: datetime) -> None:
        """
        Constructor of observation's class.

        Attributes:
        ra: Right ascension angle of observation [deg]
        dec: Declination angle of observation [deg]
        time: Time of observation [UTC]
        """
        self._ra = ra
        self._dec = dec
        self._time = time
        self._sidereal_time = None

    @property
    def ra(self):
        return self._ra

    @property
    def dec(self):
        return self._dec

    @property
    def time(self):
        return self._time

    @property
    def sidereal_time(self):
        return self._sidereal_time

    @sidereal_time.setter
    def sidereal_time(self, value):
        self._sidereal_time = value
