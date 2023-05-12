"""
position_calculations.py

This script contains functions for performing calculations related to the
positions of celestial bodies in a star system.

The script includes three main functions:

1. km_to_au(distance_km) - Converts a distance in kilometers to astronomical units (AU).
2. calculate_position(body, time) - Calculates the x, y, z position of a celestial body
at a given time.
3. compute_max_apogee(celestial_bodies) - Computes the maximum apogee of all celestial bodies.

These functions are used to calculate and update the positions of celestial bodies in the star
system over time, taking into account the bodies' orbits, inclinations, and other factors.

The script is meant to be used in conjunction with the cosmology_initialisation.py script,
which initializes the celestial bodies and their attributes. After initialization, the
position_calculations.py script can be used to update the bodies' positions and perform other
calculations as needed.
"""

import numpy as np

def km_to_au(distance_km):
    """
    Convert a distance in kilometers to astronomical units (AU).

    Parameters:
    ----------
    distance_km : float
        The distance in kilometers to be converted to astronomical units.

    Returns:
    -------
    float
        The distance in astronomical units.
    """
    # Divide the distance in km by the number of kilometers in one AU (approx. 149,597,870.7 km)
    return distance_km / 149_597_870.7


def calculate_position(body, time):
    """
    Calculate the x, y, z position of a celestial body at a given time.

    This function recursively calculates the position of a celestial body in 3D space
    at a given time, taking into account the body's parent (if any), orbit radius,
    orbit period, and inclination. It returns the x, y, z coordinates of the body's
    position as a numpy array.

    Parameters:
    ----------
    body (CelestialBody): A celestial body object.
    time (float): The time at which to calculate the body's position, in arbitrary units.

    Returns:
    -------
    np.ndarray
        The x, y, z coordinates of the celestial body's position.
    """
    # Calculate the position of the parent body, if any
    if body.parent:
        parent_x, parent_y, parent_z = calculate_position(body.parent, time)
    else:
        parent_x, parent_y, parent_z = 0, 0, 0

    # Calculate the semi-major axis and eccentricity of the body's orbit
    semi_major_axis = (body.perigee + body.apogee) / 2
    if body.apogee != 0 and body.perigee != 0:
        eccentricity = (body.apogee - body.perigee) / (body.perigee + body.apogee)
    else:
        return np.array([0, 0, 0])

    # Calculate the angle of the body in its orbit at the given time
    if body.orbit_period != 0:
        angle = 2 * np.pi * (time % body.orbit_period) / body.orbit_period
    else:
        angle = 0

    # Calculate the distance of the body from its parent at the given time
    true_anomaly_dist = semi_major_axis * (1 - eccentricity**2) / (1 + eccentricity * np.cos(angle))

    # Calculate the x, y, z position of the body relative to its parent
    x = parent_x + true_anomaly_dist * np.cos(angle)
    y = parent_y + true_anomaly_dist * np.sin(angle)
    z = parent_z + true_anomaly_dist * np.sin(np.radians(body.inclination))

    return np.array([x, y, z])


def compute_max_apogee(celestial_bodies):
    """
    Calculate the maximum apogee across all celestial bodies in the system.

    Parameters:
    ----------
    celestial_bodies : dict
        A dictionary mapping names of celestial bodies to CelestialBody objects.

    Returns:
    -------
    float
        The maximum apogee value.
    """
    return max(body.apogee for body in celestial_bodies.values())

# def compute_max_apogee(celestial_bodies):
#     """
#     Compute the maximum apogee of all celestial bodies.

#     Parameters:
#     ----------
#     celestial_bodies (dict): a dictionary of celestial bodies.

#     Returns:
#     -------
#     float
#         The maximum apogee.
#     """
#     max_apogee = max(km_to_au(body.orbit_radius) for body in celestial_bodies.values())
#     return max_apogee
