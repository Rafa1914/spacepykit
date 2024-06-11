class TLE:
    """
    Class to encapsulate the representation of a TLE.
    """

    def __init__(
        self, a: float, i: float, ra: float, e: float, w: float, theta: float
    ) -> None:
        self._a = a
        self._i = i
        self._ra = ra
        self._e = e
        self._w = w
        self._theta = theta

    @property
    def a(self):
        """Semi-major axis [km]"""
        return self._a

    @property
    def i(self):
        """Orbit inclination [°]"""
        return self._i

    @property
    def ra(self):
        """Right ascension of the ascending node [°]"""
        return self._ra

    @property
    def e(self):
        """Eccentricity [-]"""
        return self._e

    @property
    def w(self):
        """Argumento of perigee [°]"""
        return self._w

    @property
    def theta(self):
        """True anomaly [°]"""
        return self._theta

    @property
    def to_list(self):
        """List containing [a,i,ra,e,w,theta]"""
        return [self._a, self._i, self._ra, self._e, self._w, self._theta]
