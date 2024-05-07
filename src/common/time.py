from datetime import datetime

def julian_day(date:datetime) -> float:
    '''
    Calculate the julian day considering Boulet (1991) formulation.

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

    Attributes:
    year: Year between 1901 and 2099
    month: Month between 1 and 12
    day: Day between 1 and 31

    Returns:
    j0: Julian day at 0h UT.
    '''
    j0 = 367*year - int((7*(year + int((month + 9)/12)))/4) + int((275*month)/9) + day + 1721013.5 # Boulet (1991) - Methods of Orbit Detemination for the Microcomputer
    return j0

def sidereal_time(date:datetime,longitude:float) -> float:
    '''
    Calculate the sidereal time for a longitude in a specific date and time.

    Attributes:
    date: Date and time at UTC
    longitude: Site's longitude [deg]

    Returns:
    theta: Sidereal time [deg].
    '''
    j0 = julian_day_0h_boulet(date.year,date.month,date.day)
    j2000 = julian_day_0h_boulet(year=2000,month=1,day=1)
    t0 = (j0-j2000)/36525 #Centuries passed since J2000
    theta_g0 = 100.4606184 + 36000.77004*t0 + 0.000387933*t0**2 - 2.583e-8*t0**3 # Seidelmann, 1992 - Explanatory Supplement to the Astronomical Almanac
    ut = date.hour + date.minute/60 + date.second/3600
    theta_g = theta_g0 + 360.98564724*ut/24
    theta = theta_g + longitude
    if 0 <= theta <= 360:
        return theta
    else:
        if theta > 360:
            theta = theta - 360
            return theta
        else:
            theta = 360 + theta
            return theta