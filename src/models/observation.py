from datetime import datetime
from src.models.observatory import Observatory
from src.orbit.orbit_determination import sidereal_time


class Observation:
    def __init__(self, ra: float, dec: float, time_utc: datetime) -> None:
        """
        Constructor of observation's class.

        Attributes:
        ra: Right ascension angle of observation [deg]
        dec: Declination angle of observation [deg]
        time: Time of observation [UTC]
        """
        self._ra = ra
        self._dec = dec
        self._time_utc = time_utc
        self._time_seconds = time_utc.timestamp()
        self._sidereal_time = None
        self._observatory = None

    @property
    def ra(self):
        return self._ra

    @property
    def dec(self):
        return self._dec

    @property
    def time_utc(self):
        return self._time_utc

    @property
    def time_seconds(self):
        return self._time_seconds

    @property
    def sidereal_time(self):
        return self._sidereal_time

    @property
    def observatory(self):
        return self._observatory

    @observatory.setter
    def observatory(self, value: Observatory):
        self._sidereal_time = sidereal_time(value.longitude, date=self.time_utc)
        self._observatory = value
