from datetime import datetime
from src.common.utils import *
from src.common.constants import *

def julian_day(date:datetime) -> float:
    '''
    Calculate the julian day considering Boulet (1991) formulation.

    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
    date: Date and time at UTC

    Returns:
    jd: Julian day.
    '''
    j0 = julian_day_0h_boulet(date.year,date.month,date.day)
    ut = date.hour + date.minute/60 + date.second/3600
    jd = j0 + ut/24
    return jd

def julian_day_0h_boulet(year,month,day) -> float:
    '''
    Calculate the julian day at 0h UT considering Boulet (1991) formulation.

    Reference: Boulet (1991) - Methods of Orbit Detemination for the Microcomputer

    Attributes:
    year: Year between 1901 and 2099
    month: Month between 1 and 12
    day: Day between 1 and 31

    Returns:
    j0: Julian day at 0h UT.
    '''
    j0 = 367*year - int((7*(year + int((month + 9)/12)))/4) + int((275*month)/9) + day + 1721013.5 # Boulet (1991) - Methods of Orbit Detemination for the Microcomputer
    return j0

def sidereal_time_greenwich_0h_seidelmann(date:datetime) -> float:
    '''
    Calculate the sidereal time on Greenwich Meriridian for a date at 0h UT using Seidelmann (1992) formulation.
    
    Reference: Seidelmann, 1992 - Explanatory Supplement to the Astronomical Almanac
    '''
    j0 = julian_day_0h_boulet(date.year,date.month,date.day)
    t0 = (j0-J2000)/36525 #Centuries passed since J2000
    theta_g0 = 100.4606184 + 36000.77004*t0 + 0.000387933*t0**2 - 2.583e-8*t0**3 # Seidelmann, 1992 - Explanatory Supplement to the Astronomical Almanac
    return convert_angle_0_to_360_range(theta_g0)

def sidereal_time_greenwich(date:datetime) -> float:
    '''
    Calculate the sidereal time on Greenwich Meriridian for a date and time.

    Reference: Curtis - Orbital Mechanics for Engineering Students
    '''
    theta_g0 = sidereal_time_greenwich_0h_seidelmann(date)
    ut = date.hour + date.minute/60 + date.second/3600
    theta_g = theta_g0 + 360.98564724*ut/24
    return convert_angle_0_to_360_range(theta_g)

def sidereal_time(longitude:float, date:datetime = None, theta_g = None) -> float:
    '''
    Calculate the sidereal time for a longitude in a specific date and time.
    
    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
    longitude: Site's longitude [deg]
    date: Date and time at UTC
    theta_g: Siedereal time at Greenwich [deg]

    Returns:
    theta: Sidereal time [deg].
    '''
    # Validation
    if date is None and theta_g is None:
        raise Exception('Either date or theta_g parameter must have a value.')
    #Calculation
    if theta_g is None:
        theta_g = sidereal_time_greenwich(date)
    theta = theta_g + longitude
    return convert_angle_0_to_360_range(theta)