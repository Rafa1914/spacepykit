from src.common.constants import *
from src.common.vector import *
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
        return (1-np.cos(np.sqrt(z)))/z
    elif z < 0:
        return (np.cosh(np.sqrt(-z))-1)/(-z)
    else:
        return 0.5


def find_universal_anomaly(dt:float,r:float,vr:float,alpha:float) -> float:
    '''
    Calculates the universal anomaly for a given time, state vector and semimajor axis.
    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
    dt: Time since the universal variable is 0 (Could assume as time since perigee) [s]
    r: Radius at time t0 (Could assume radius of perigee) [km]
    v: Velocity at time t0 (Could assume velocity of perigee) [km/h]
    alpha: Reciprocal of Semimajor axis [km^-1]

    Returns:
    x: Universal anomaly (km^1/2)
    '''
    mu_sqrt = np.sqrt(EARTH_GRAVITATIONAL_PARAMETER)
    f = lambda x : (r*vr*pow(x,2)*stumpff_C(alpha*pow(x,2)))/mu_sqrt + (1-alpha*r)*pow(x,3)*stumpff_S(alpha*pow(x,2)) + r*x - dt*mu_sqrt
    x0 = mu_sqrt*np.abs(alpha)*dt
    x = fsolve(f,x0)
    return x[0]

def lagrange_f_and_g_from_x(x:float,r:float,dt:float,alpha:float) -> tuple[float,float]:
    '''
    Calculates de Lagrange's coefficients from the universal anomaly.
    The coefficients are the f and g of the expression: 
    r(t-t0) = f*r0 + g*v0
    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
    x: Universal anomaly [km^1/2]
    r: Magnitude of r0 [km]
    dt: Time interval since t0
    alpha: Reciprocal of the semimajor axis [km^-1]

    Returns:
    f: Coefficient [-]
    g: Coefficient [s]
    '''
    z = alpha*pow(x,2)
    f = 1 - pow(x,2)*stumpff_C(z)/r
    g = dt - pow(x,3)*stumpff_S(z)/np.sqrt(EARTH_GRAVITATIONAL_PARAMETER)
    return (f,g)

def lagrange_fdot_and_gdot_from_x(x:float,r:float,r0:float,dt:float,alpha:float) -> tuple[float,float]:
    '''
    Calculates de Lagrange's coefficients from the universal anomaly.
    The coefficients are the f and g of the expression: 
    v(t-t0) = f_dot*r0 + g_dot*v0
    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
    x: Universal anomaly [km^1/2]
    r: Magnitude of r [km]
    r0: Magnitude of r0 [km]
    dt: Time interval since t0
    alpha: Reciprocal of the semimajor axis [km^-1]

    Returns:
    f: Coefficient [-]
    g: Coefficient [s]
    '''
    z = alpha*pow(x,2)
    f_dot = np.sqrt(EARTH_GRAVITATIONAL_PARAMETER)*(alpha*pow(x,3)*stumpff_S(z)-x)/(r*r0)
    g_dot = 1 - pow(x,2)*stumpff_C(z)/r
    return (f_dot,g_dot)

def rv_from_r0v0(r0:Vector,v0:Vector,dt:float) -> tuple[Vector,Vector]:
    """This function computes the state vector (R,V) from the initial state vector (R0,V0) and the elapsed time.

    Args:
        r0 (Vector): Initial position vector
        v0 (Vector): Initial velocity vector
        dt (float): Elapsed time [s]

    Returns:
        tuple[Vector,Vector]: State vector resultant
    """
    vr = find_projection(v0,r0)
    alpha = 2/r0.magnitude - dot_product(v0,v0)/EARTH_GRAVITATIONAL_PARAMETER
    x = find_universal_anomaly(dt,r0.magnitude,vr,alpha)
    f,g = lagrange_f_and_g_from_x(x,r0.magnitude,dt,alpha)
    r = r0*f + v0*g
    f_dot,g_dot = lagrange_fdot_and_gdot_from_x(x,r.magnitude,r0.magnitude,dt,alpha)
    v = r0*f_dot + v0*g_dot
    return(r,v)


