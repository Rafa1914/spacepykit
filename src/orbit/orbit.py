from src.common.constants import *
from scipy.optimize import fsolve
import numpy as np

def stumpff_S(z:float) -> float:
    '''
    Calculates de Stumpff function S(z).
    Reference: Reference: Curtis - Orbital Mechanics for Engineering Students
    '''
    if z > 0:
        z_sqrt = np.sqrt(z)
        return (z_sqrt-np.sin(z_sqrt))/(pow(z_sqrt,3))
    elif z < 0:
        z_sqrt = np.sqrt(-z)
        return (np.sinh(z_sqrt)-z_sqrt)/(pow(z_sqrt,3))
    else:
        return 1/6

def stumpff_C(z:float) -> float:
    '''
    Calculates de Stumpff function C(z).
    Reference: Reference: Curtis - Orbital Mechanics for Engineering Students
    '''
    if z > 0:
        return (1-np.cos(z))/z
    elif z < 0:
        return (np.cosh(np.sqrt(-z))-1)/(-z)
    else:
        return 0.5


def find_universal_anomaly(dt:float,r:float,vr:float,a:float) -> float:
    '''
    Calculates the universal anomaly for a given time, state vector and semimajor axis.
    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
    dt: Time since the universal variable is 0 (Could assume as time since perigee) [s]
    r: Radius at time t0 (Could assume radius of perigee) [km]
    v: Velocity at time t0 (Could assume velocity of perigee) [km/h]
    a: Semimajor axis [km]

    Returns:
    x: Universal anomaly (km^1/2)
    '''
    mu_sqrt = np.sqrt(EARTH_GRAVITATIONAL_PARAMETER)
    alpha = 1/a
    f = lambda x : (r*vr*pow(x,2)*stumpff_C(alpha*pow(x,2)))/mu_sqrt + (1-alpha*r)*pow(x,3)*stumpff_S(alpha*pow(x,2)) + r*x - dt*mu_sqrt
    x0 = mu_sqrt*np.abs(alpha)*dt
    x = fsolve(f,x0)
    return x
