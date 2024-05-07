from src.common.vector import *
from datetime import datetime

def radec_from_r(r : Vector) -> tuple[float,float]:
    '''
    Calculates de Right Ascension and Declination for the vector r in its reference frame.

    Atrributes:
    r: Position vector

    Returns:
    alpha: Right ascension angle [rad]
    ddelta: Declination angle [rad]
    '''
    cossine_direction_x = r.x/r.magnitude
    cossine_direction_y = r.y/r.magnitude
    cossine_direction_z = r.z/r.magnitude

    delta = np.arcsin(cossine_direction_z)
    alpha = np.arccos(cossine_direction_x/np.cos(delta)) if cossine_direction_y > 0 else 2*np.pi - np.arccos(cossine_direction_x/np.cos(delta))
    return (alpha,delta)

def julian_day(date:datetime) -> float:
    '''
    Calculate the julian day considering Boulet (1991) formulation.

    Attributes:
    date: Date and time at UTC

    Returns:
    jd: Julian day.
    '''
    year = date.year
    month = date.month
    day = date.day
    j0 = 367*year - int((7*(year + int((month + 9)/12)))/4) + int((275*month)/9) + day + 1721013.5 # Boulet (1991) - Methods of Orbit Detemination for the Microcomputer
    ut = date.hour + date.minute/60 + date.second/3600
    jd = j0 + ut/24
    return jd