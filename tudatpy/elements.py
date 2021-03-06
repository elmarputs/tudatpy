# from . import _orbital_element_conversions as _oec
from tudatpy import _orbital_element_conversions as _oec
import numpy as np

from astropy import units as u

m_s = u.m / u.s
m = u.m
rad = u.rad


@u.quantity_input
def spherical2cartesian(r: m,
                        lat: rad,
                        lon: rad,
                        speed: m_s,
                        fpa: rad,
                        heading: rad) -> np.ndarray:
    """
    Function to convert spherical state to cartesian.

    Parameters
    ----------
    r : float
        Position vector magnitude (m).
    lat : float
        Latitude (m).
    lon : float
        Longitude (m).
    speed : float
        Magnitude of velocity (m).
    fpa : float
        Flight path angle (rad).
    heading : float
        Heading angle (rad).

    Returns
    -------
    cartesian_state : ndarray
        Cartesian state represented as [Rx, Ry, Rz, Vx, Vz, Vy] with distance in (m) and speed in (m/s).

    """
    spherical_idx = _oec.SphericalOrbitalStateElementIndices
    spherical_state = np.zeros(6)
    spherical_state[int(spherical_idx.radius_index)] = r
    spherical_state[int(spherical_idx.latitude_index)] = lat
    spherical_state[int(spherical_idx.longitude_index)] = lon
    spherical_state[int(spherical_idx.speed_index)] = speed
    spherical_state[int(spherical_idx.flight_path_index)] = fpa
    spherical_state[int(spherical_idx.heading_angle_index)] = heading
    return _oec.convert_spherical_orbital_to_cartesian_state(spherical_state)


def keplerian2cartesian(mu, a, ecc, inc, raan, argp, nu):
    """
    Function to convert Keplerian state to cartesian.

    Parameters
    ----------
    mu : float
        Standard gravitational parameter (m^3 / s^2).
    a : float
        Semi-major axis (m).
    ecc : float
        Eccentricity (-).
    inc : float
        Inclination (rad).
    raan : float
        Right Ascension of the Ascending Node (rad).
    argp : float
        Argument of Perigee (rad).
    nu : float
        True anomaly (rad).

    Returns
    -------
    cartesian_state : ndarray
        Cartesian state represented as [Rx, Ry, Rz, Vx, Vz, Vy] with distance in (m) and speed in (m/s).

    """
    keplerian_idx = _oec.KeplerianElementIndices
    keplerian_state = np.zeros(6)
    keplerian_state[int(keplerian_idx.semi_major_axis_index)] = a
    keplerian_state[int(keplerian_idx.eccentricity_index)] = ecc
    keplerian_state[int(keplerian_idx.inclination_index)] = inc
    keplerian_state[int(keplerian_idx.longitude_of_ascending_node_index)] = raan
    keplerian_state[int(keplerian_idx.argument_of_periapsis_index)] = argp
    keplerian_state[int(keplerian_idx.true_anomaly_index)] = nu
    return _oec.convert_keplerian_to_cartesian_elements(keplerian_state, mu)


if __name__ == "__main__":
    res = spherical2cartesian(100, 2, 2, 1000, 1, 1)
    print(res)