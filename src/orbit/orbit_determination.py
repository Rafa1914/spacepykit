from src.models.tle import TLE
from src.common.vector import *
from src.common.reference_frame import *
from src.common.constants import *
from src.orbit.orbit import *
import numpy.polynomial.polynomial as poly
from src.models.observatory import Observatory
from src.models.observation import Observation
import pandas as pd


def r2v2_from_angles_only_gauss(
    latitude_geodetic: float,
    altitude: float,
    local_sidereal_times: list[float],
    RAs: list[float],
    DECs: list[float],
    times: list[float],
    perform_improvement: bool = True,
) -> tuple[Vector, Vector]:
    """
    Perform the Gauss angles-only method to provide the position and velocity vectors at second observation
    in a three set observation.
    Reference: Curtis - Orbital Mechanics for Engineering Students

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
    """
    r2, v2, R, rho, [f1, g1, f3, g3], D0, D = r2v2_from_angles_only_gauss_core(
        latitude_geodetic, altitude, local_sidereal_times, RAs, DECs, times
    )

    if perform_improvement:
        r2, v2 = improve_r2v2_from_gauss_angles_only(
            r2, v2, R, rho, f1, g1, f3, g3, D0, D, times
        )

    return (r2, v2)


def r2v2_from_angles_only_gauss_core(
    latitude_geodetic: float,
    altitude: float,
    local_sidereal_times: list[float],
    RAs: list[float],
    DECs: list[float],
    times: list[float],
) -> tuple[Vector, Vector, list[Vector], list[Vector], list[float], float, list[list]]:
    """
    Perform the Gauss angles-only core method to provide an initial estimation of the position and velocity vectors at second observation
    in a three set observation.
    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
        latitude_geodetic: Site's geodetic latitude [deg]
        altitude: Site's altitude [km]
        local_sidereal_times: Site's sidereal local times of observation [deg]
        RAs: Right ascension angles from observation [deg]
        DECs: Declination angles from observation [deg]
        times: Times of observation

    Returns:
        r2: Object's position vector on Geocentrical Equatorial Reference Frame at second observation [km]
        v2: Object's velocity vector on Geocentrical Equatorial Reference Frame at second observation [km/s]
        R: Origin of Topocentric Equatorial Reference Frame at each observation [km]
        rho: Object's position vector on Topocentric Equatorial Reference Frame at each observation [km]
        [f1,g1,f3,g3]: Lagrange's estimated coefficients
        D0: Auxiliary Constant
        D: Auxiliary Matrix
    """
    R = []  # Topocentric reference frame origin
    versor_rho = []  # Relative position versor

    for i in range(0, 3):
        R.append(
            get_topocentric_origin(
                latitude=latitude_geodetic,
                altitude=altitude,
                theta=local_sidereal_times[i],
                is_geocentric_latitude=False,
            )
        )
        versor_rho.append(versor_from_radec(ra=RAs[i], dec=DECs[i]))

    # Time's interval:
    tau_1 = times[0] - times[1]
    tau_3 = times[2] - times[1]
    tau = times[2] - times[0]

    p = []  # Auxiliary variable
    p.append(cross_product(versor_rho[1], versor_rho[2]))
    p.append(cross_product(versor_rho[0], versor_rho[2]))
    p.append(cross_product(versor_rho[0], versor_rho[1]))

    D0 = dot_product(versor_rho[0], p[0])

    # Constants D's:
    D = []
    for i in range(0, 3):
        row = []
        for j in range(0, 3):
            row.append(dot_product(R[i], p[j]))
        D.append(row)

    # Constants A and B:
    A = (-D[0][1] * (tau_3 / tau) + D[1][1] + D[2][1] * (tau_1 / tau)) / D0
    B = (1 / (6 * D0)) * (
        D[0][1] * (tau_3**2 - tau**2) * tau_3 / tau
        + D[2][1] * (tau**2 - tau_1**2) * tau_1 / tau
    )

    # Constant E:
    E = dot_product(R[1], versor_rho[1])

    # Polynomial's coeficients -> (r2)^8 + a(r2)^6 + b(r2)^3 + c:
    mu = EARTH_GRAVITATIONAL_PARAMETER  # Earth graviational constant [km^3/s^2]
    a = -(A**2 + 2 * A * E + R[1].magnitude ** 2)
    b = -2 * mu * B * (A + E)
    c = -((mu * B) ** 2)

    # Polynomial's roots:
    coef = [c, 0, 0, b, 0, 0, a, 0, 1]
    roots = poly.polyroots(coef)

    # Initial estimation of r2 magnitude:
    condition = (np.isreal(roots)) & (roots > 0)
    roots_filtered = roots[condition]
    r2_star = np.real(roots_filtered[0])

    # Relative position's:
    rho1_mod = (
        (
            6 * (D[2][0] * tau_1 / tau_3 + D[1][0] * tau / tau_3) * r2_star**3
            + mu * D[2][0] * (tau**2 - tau_1**2) * tau_1 / tau_3
        )
        / (6 * r2_star**3 + mu * (tau**2 - tau_3**2))
        - D[0][0]
    ) / D0
    rho2_mod = A + mu * B / (r2_star**3)
    rho3_mod = (
        (
            6 * (D[0][2] * tau_3 / tau_1 - D[1][2] * tau / tau_1) * r2_star**3
            + mu * D[0][2] * (tau**2 - tau_3**2) * tau_3 / tau_1
        )
        / (6 * r2_star**3 + mu * (tau**2 - tau_1**2))
        - D[2][2]
    ) / D0
    rho_mod = [rho1_mod, rho2_mod, rho3_mod]
    rho = []
    for i in range(0, 3):
        rho.append(
            Vector(
                x=rho_mod[i] * versor_rho[i].x,
                y=rho_mod[i] * versor_rho[i].y,
                z=rho_mod[i] * versor_rho[i].z,
            )
        )

    # r's calculus:
    r = []  # Position vector on geocentrical equatorial reference frame
    for i in range(0, 3):
        r.append(R[i] + rho[i])

    # Lagrange's coefficients:
    f1 = 1 - 0.5 * mu * (tau_1**2) / (r2_star**3)
    f3 = 1 - 0.5 * mu * (tau_3**2) / (r2_star**3)
    g1 = tau_1 - (1 / 6) * mu * (tau_1**3) / (r2_star**3)
    g3 = tau_3 - (1 / 6) * mu * (tau_3**3) / (r2_star**3)

    # Velocity:
    v2 = (r[2] * f1 - r[0] * f3) * ((f1 * g3 - f3 * g1) ** (-1))

    return (r[1], v2, R, rho, [f1, g1, f3, g3], D0, D)


def classic_orbital_elements_from_rv(
    r: Vector, v: Vector
) -> tuple[float, float, float, float, float, float]:
    """
    From the state vector, get the classical orbital elements for the object's orbit.
    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
    r: Position vector [km]
    v: Velocity vector [km/s]

    Returns:
    h: Specific angular moment [km2/s]
    i: Inclination [deg]
    ra_ascending_node: Rigth ascension of the ascending node [deg]
    e: Eccentricity [-]
    w: Argument of perigee [deg]
    theta: True anomaly [deg]
    """
    v_r = dot_product(r, v) / r.magnitude  # Radial Velocity [km/s]
    h_vector = cross_product(r, v)  # Specific angular momentum vector
    i = np.arccos(h_vector.z / h_vector.magnitude)  # Inclination
    N = cross_product(Vector(x=0, y=0, z=1), h_vector)  # Node line vector
    e_vector = (
        r * (v.magnitude**2 - EARTH_GRAVITATIONAL_PARAMETER / r.magnitude)
        - v * r.magnitude * v_r
    ) * (
        1 / EARTH_GRAVITATIONAL_PARAMETER
    )  # Eccentricity vector
    # Rigth ascension of the ascending node
    if N.y >= 0:
        ra_ascending_node = np.arccos(
            N.x / N.magnitude,
        )
    else:
        ra_ascending_node = 2 * np.pi - np.arccos(N.x / N.magnitude)
    # Argument of perigee
    if e_vector.z >= 0:
        w = np.arccos(dot_product(N, e_vector) / (N.magnitude * e_vector.magnitude))
    else:
        w = 2 * np.pi - np.arccos(
            dot_product(N, e_vector) / (N.magnitude * e_vector.magnitude)
        )
    # True anomaly
    if v_r >= 0:
        theta = np.arccos(dot_product(e_vector, r) / (e_vector.magnitude * r.magnitude))
    else:
        theta = 2 * np.pi - np.arccos(
            dot_product(e_vector, r) / (e_vector.magnitude * r.magnitude)
        )
    return (
        h_vector.magnitude,
        np.rad2deg(i),
        np.rad2deg(ra_ascending_node),
        e_vector.magnitude,
        np.rad2deg(w),
        np.rad2deg(theta),
    )


def semimajor_axis_from_he(h: float, e: float) -> float:
    """
    Get the semimajor axis for an elliptical orbit from its specific angular moment and eccentricity.
    Reference: Curtis - Orbital Mechanics for Engineering Students

    Attributes:
    h: Specific angular momentum [km^2/s]
    e: Eccentricity

    Returns:
    a: Semimajor axis [km]
    """
    mu = EARTH_GRAVITATIONAL_PARAMETER
    rp = (pow(h, 2) / mu) * (1 / (1 + e))
    ra = (pow(h, 2) / mu) * (1 / (1 - e))
    a = 0.5 * (rp + ra)
    return a


def improve_r2v2_from_gauss_angles_only(
    r2: Vector,
    v2: Vector,
    R: list[Vector],
    rho: list[Vector],
    f1: float,
    g1: float,
    f3: float,
    g3: float,
    D0: float,
    D: list[list],
    t: list[float],
) -> tuple[Vector, Vector]:
    # Time's interval:
    tau_1 = t[0] - t[1]
    tau_3 = t[2] - t[1]

    # Improvement of r2,v2:
    diff = [1, 1, 1]
    n = 0
    n_max = 1000
    tol = 1e-8

    while (max(diff) > tol) and (n < n_max):
        # Iteration's number control:
        n = n + 1

        # Old parameters:
        rho_old = list(rho)
        f1_old = f1
        f3_old = f3
        g1_old = g1
        g3_old = g3

        # Initial Constants
        alpha = (2 / r2.magnitude) - pow(
            v2.magnitude, 2
        ) / EARTH_GRAVITATIONAL_PARAMETER
        vr = dot_product(v2, r2) / r2.magnitude
        x1 = find_universal_anomaly(tau_1, r2.magnitude, vr, alpha)
        x3 = find_universal_anomaly(tau_3, r2.magnitude, vr, alpha)
        f1, g1 = lagrange_f_and_g_from_x(x1, r2.magnitude, tau_1, alpha)
        f3, g3 = lagrange_f_and_g_from_x(x3, r2.magnitude, tau_3, alpha)
        f1 = (f1 + f1_old) / 2
        f3 = (f3 + f3_old) / 2
        g1 = (g1 + g1_old) / 2
        g3 = (g3 + g3_old) / 2
        c1 = g3 / (f1 * g3 - f3 * g1)
        c3 = -g1 / (f1 * g3 - f3 * g1)

        # Obtaining the new rhos:
        rho1_mod = (1 / D0) * (-D[0][0] + (1 / c1) * D[1][0] - (c3 / c1) * D[2][0])
        rho2_mod = (1 / D0) * (-c1 * D[0][1] + D[1][1] - c3 * D[2][1])
        rho3_mod = (1 / D0) * ((-c1 / c3) * D[0][2] + (1 / c3) * D[1][2] - D[2][2])
        rho_mod = [rho1_mod, rho2_mod, rho3_mod]
        rho = []
        for i in range(0, 3):
            rho.append(
                Vector(
                    x=rho_mod[i] * rho_old[i].x / rho_old[i].magnitude,
                    y=rho_mod[i] * rho_old[i].y / rho_old[i].magnitude,
                    z=rho_mod[i] * rho_old[i].z / rho_old[i].magnitude,
                )
            )

        # Obatinig the new position's vector on geocentrical equatorial reference frame:
        r = []
        for i in range(0, 3):
            r.append(R[i] + rho[i])

        # Obtaining the new state vector:
        r2 = r[1]
        v2 = (r[2] * f1 - r[0] * f3) * ((f1 * g3 - f3 * g1) ** (-1))

        # Verifying if the loop is over:
        diff = []
        for i in range(0, 3):
            diff.append(abs(rho_old[i].magnitude - rho[i].magnitude))
    return r2, v2


def get_observations(
    tdm: pd.DataFrame, observatory: Observatory, last_position: int
) -> list[Observation]:
    """Generate the observations from RA and DEC informations

    Args:
        tdm (pd.DataFrame): _description_
        observatory (Observatory): _description_
        last_position (int): _description_

    Returns:
        list[Observation]: _description_
    """
    # Posições Analisadas:
    initial = 0
    middle = int(last_position / 2)
    end = last_position - 1

    observations = []

    for i in [initial, middle, end]:
        ra = tdm["RA[deg]"].iloc[i]
        dec = tdm["DEC[deg]"].iloc[i]
        time_utc = tdm["Tempo (UTC)"].iloc[i]
        obs = Observation(ra, dec, time_utc)
        obs.observatory = observatory
        observations.append(obs)
    return observations


def get_state_vector(
    observations: list[Observation], observatory: Observatory
) -> tuple[Vector, Vector]:
    # Posições Analisadas:
    RAs = [observations[0].ra, observations[1].ra, observations[2].ra]
    DECs = [observations[0].dec, observations[1].dec, observations[2].dec]
    local_sidereal_time = [
        observations[0].sidereal_time,
        observations[1].sidereal_time,
        observations[2].sidereal_time,
    ]
    t = [
        observations[0].time_seconds,
        observations[1].time_seconds,
        observations[2].time_seconds,
    ]
    # Estimativa inicial do vetor de estado do objeto (r2,v2):
    r2, v2 = r2v2_from_angles_only_gauss(
        latitude_geodetic=observatory.latitude_geodetic,
        altitude=observatory.altitude,
        local_sidereal_times=local_sidereal_time,
        RAs=RAs,
        DECs=DECs,
        times=t,
    )
    return (r2, v2)


def get_orbital_parameters(r: Vector, v: Vector) -> TLE:

    # Obtenção dos elementos orbitais:
    h, i, ra, e, w, theta = classic_orbital_elements_from_rv(r, v)

    # Obtenção do semi-eixo maior:
    a = semimajor_axis_from_he(h, e)

    return TLE(a, i, ra, e, w, theta)

def get_iod(tdm: pd.DataFrame, observatory: Observatory, index: int = None):
    index = len(tdm["Tempo[s]"]) if index is None else index
    obs = get_observations(tdm, observatory, index)
    r, v = get_state_vector(obs, observatory)
    tle = get_orbital_parameters(r, v)
    return (obs, r, v, tle)

def get_iod_multiple_radec(tdm: pd.DataFrame, observatory: Observatory) -> pd.DataFrame:
    dt = []
    a = []
    i = []
    omega = []
    w = []
    e = []
    for index in range(1, len(tdm["Tempo[s]"])):
        try:
            obs, _, _, tle = get_iod(tdm, observatory, index+1)
            a.append(tle.a)
            e.append(tle.e)
            omega.append(tle.ra)
            w.append(tle.w)
            i.append(tle.i)
            dt.append(obs[2].time_seconds - obs[0].time_seconds)
        except:
            continue
    data = {"a": a, "e": e, "omega": omega, "w": w, "i": i, "dt": dt}
    df = pd.DataFrame(data=data)
    # Filtro de Outliers:
    df.dropna(inplace=True)
    df.where((df["e"] < 0.8), inplace=True)
    return df