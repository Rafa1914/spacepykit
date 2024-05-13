from src.common.vector import *
from src.common.constants import EARTH_OBLATENESS,EARTH_RADIUS
from src.common.time import sidereal_time
from datetime import datetime

def radec_from_r(r : Vector) -> tuple[float,float]:
    '''
    Calculates de Right Ascension and Declination for the vector r in its reference frame.

    Reference: Curtis - Orbital Mechanics for Engineering Students

    Atrributes:
    r: Position vector

    Returns:
    alpha: Right ascension angle [rad]
    delta: Declination angle [rad]
    '''
    cossine_direction_x = r.x/r.magnitude
    cossine_direction_y = r.y/r.magnitude
    cossine_direction_z = r.z/r.magnitude

    delta = np.arcsin(cossine_direction_z)
    alpha = np.arccos(cossine_direction_x/np.cos(delta)) if cossine_direction_y > 0 else 2*np.pi - np.arccos(cossine_direction_x/np.cos(delta))
    return (alpha,delta)

def get_geodetic_latitude(latitude:float) -> float:
    '''
    Calculates the geodetic latitude for a given geocentric latitude.

    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
    latitude: Latitude location [deg]

    Returns:
    Returns the geodetic latitude [rad]
    '''
    phi = np.deg2rad(latitude)
    f = EARTH_OBLATENESS
    geodetic_latitude = np.arctan(np.tan(phi)*(1-f)**2)
    return geodetic_latitude

def get_topocentric_origin(longitude:float,latitude:float,altitude:float,date:datetime = None,theta_g: float = None) -> Vector:
    '''
    Get the topocentric coordinate system origin in the equatorial geocentric reference frame.

    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
    longitude: Longitude location [deg]
    latitude: Latitude location [deg]
    altitude: Altitude location [km]
    date: Date and time [UT]
    theta_g: Sidereal time at Greenwich [deg] 

    Returns:
    The relative position vector of the origin in equatorial geocentric reference frame.
    '''
    # Validation
    if date is None and theta_g is None:
        raise Exception('Either date or theta_g parameter must have a value.')
    # Calculation
    f = EARTH_OBLATENESS
    r_e = EARTH_RADIUS
    phi = get_geodetic_latitude(latitude)
    r_phi = r_e/(np.sqrt(1-(2-f)*f*(np.sin(phi))**2)) # Seidelmann, P.K., 1992. Explanatory Supplement to the Astronomical Almanac
    r_c = r_phi + altitude
    r_s = np.power(1-f,2)*r_phi + altitude
    theta = np.deg2rad(sidereal_time(longitude=longitude,theta_g=theta_g,date=date))
    r = Vector(x = r_c*np.cos(phi)*np.cos(theta), y = r_c*np.cos(phi)*np.sin(theta), z = r_s*np.sin(phi))
    return r