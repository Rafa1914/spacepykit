from src.common.vector import *
from src.common.reference_frame import *
from src.common.constants import *
import numpy.polynomial.polynomial as poly


def r2v2_from_angles_only_gauss(
        latitude_geodetic:float,
        altitude:float,
        local_sidereal_times:list[float],
        RAs:list[float],
        DECs:list[float],
        times:list[float]) -> tuple[Vector,Vector]:
    '''
    Perform the Gauss angles-only method to provide the position and velocity vectors at second observation
    in a three set observation.

    Attributes:
        latitude_geodetic: Site's geodetic latitude [deg]
        altitude: Site's altitude [km]
        local_sidereal_times: Site's sidereal local times of observation [deg]
        RAs: Right ascension angles from observation [deg]
        DECs: Declination angles from observation [deg]
        times: Times of observation
    
    Returns:
        r2: Object's position vector on Geocentrical equatorial reference frame at second observation [km]
        v2: Object's velocity vector on Geocentrical equatorial reference frame at second observation [km/s]
    '''
    R = []           #Topocentric reference frame origin 
    versor_rho = []  #Relative position versor

    for i in range(0,3):
        R.append(get_topocentric_origin(latitude=latitude_geodetic,altitude=altitude,theta=local_sidereal_times[i],is_geocentric_latitude=False))
        versor_rho.append(versor_from_radec(ra = RAs[i], dec = DECs[i]))

    #Time's interval:
    tau_1 = times[0]-times[1]
    tau_3 = times[2]-times[1]
    tau = times[2]-times[0]

    p = [] #Auxiliary variable
    p.append(cross_product(versor_rho[1],versor_rho[2]))
    p.append(cross_product(versor_rho[0],versor_rho[2]))
    p.append(cross_product(versor_rho[0],versor_rho[1]))

    D0 = dot_product(versor_rho[0],p[0])

    #Constants D's:
    D = []
    for i in range(0,3):
        row = []
        for j in range(0,3):
            row.append(dot_product(R[i],p[j]))
        D.append(row)

    #Constants A and B:
    A = (-D[0][1]*(tau_3/tau) + D[1][1] + D[2][1]*(tau_1/tau))/D0
    B = (1/(6*D0))*(D[0][1]*(tau_3**2 - tau**2)*tau_3/tau + D[2][1]*(tau**2-tau_1**2)*tau_1/tau)

    #Constant E:
    E = dot_product(R[1],versor_rho[1])

    #Polynomial's coeficients -> (r2)^8 + a(r2)^6 + b(r2)^3 + c:
    mu = EARTH_GRAVITATIONAL_PARAMETER #Earth graviational constant [km^3/s^2]
    a = -(A**2 + 2*A*E + R[1].magnitude**2)
    b = -2*mu*B*(A+E)
    c = -(mu*B)**2

    #Polynomial's roots:
    coef = [c,0,0,b,0,0,a,0,1]
    roots = poly.polyroots(coef)

    #Initial estimation of r2 magnitude:
    condition = (np.isreal(roots)) & (roots > 0)
    roots_filtered = roots[condition]
    r2_star = np.real(roots_filtered[0])

    #Relative position's:
    rho1_mod = ((6*(D[2][0]*tau_1/tau_3 + D[1][0]*tau/tau_3)*r2_star**3 + mu*D[2][0]*(tau**2-tau_1**2)*tau_1/tau_3)/(6*r2_star**3 + mu*(tau**2-tau_3**2))-D[0][0])/D0
    rho2_mod = A + mu*B/(r2_star**3)
    rho3_mod = ((6*(D[0][2]*tau_3/tau_1 - D[1][2]*tau/tau_1)*r2_star**3 + mu*D[0][2]*(tau**2-tau_3**2)*tau_3/tau_1)/(6*r2_star**3 + mu*(tau**2-tau_1**2))-D[2][2])/D0
    rho_mod = [rho1_mod,rho2_mod,rho3_mod]
    rho = []
    for i in range(0,3):
        rho.append(Vector(x = rho_mod[i]*versor_rho[i].x, y = rho_mod[i]*versor_rho[i].y, z = rho_mod[i]*versor_rho[i].z))


    #r's calculus:
    r = [] #Position vector on geocentrical equatorial reference frame
    for i in range(0,3):
        r.append(R[i]+rho[i])

    #Lagrange's coefficients:
    f1 = 1-0.5*mu*(tau_1**2)/(r2_star**3)
    f3 = 1-0.5*mu*(tau_3**2)/(r2_star**3)
    g1 = tau_1 - (1/6)*mu*(tau_1**3)/(r2_star**3)
    g3 = tau_3 - (1/6)*mu*(tau_3**3)/(r2_star**3)

    #Velocity:
    v2 = (r[2]*f1 - r[0]*f3)*((f1*g3-f3*g1)**(-1))
    return (r[1],v2)