from src.common.vector import *

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

