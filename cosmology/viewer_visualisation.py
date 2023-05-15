"""
viewer_visualisation.py
This script contains several functions for visualising celestial bodies in a star system.

It includes the following functions:

- `calculate_moon_phase()`: This function calculates the phase of a moon based on its position
    relative to its parent planet and the star.

- `update_positions()`: This function calculates the positions of all celestial bodies in the
    star system at a given time.

- `planet_viewer()`: This function calculates and displays the visible celestial bodies in the
    sky from a specific location on a given planet, including their altitude, azimuth, and
    angular size.

To use these functions, you need to provide a star system (a dictionary where keys are the names
of celestial bodies and values are instances of the `CelestialBody` class), a current time, and
the positions of celestial bodies. For the `planet_viewer()` function, you also need to specify
the name of the planet where the observer is located, and the latitude and longitude of the
observer's location.

Please note that this script assumes that the positions of celestial bodies are provided in
astronomical units (AU), and time is provided in hours since epoch (0-23.9999). The positions
are expected to be calculated using the `calculate_position()` function from the
`position_calculations.py` script.

This script uses the numpy and math libraries and requires importing the `CelestialBody` class from
the `cosmology_initialisation.py` script.
"""

import math
import numpy as np
from position_calculations import calculate_position, km_to_au
from cosmology_initialisation import Skybox

# Minimum angular size for an object to be visible (1m3 object at a distance
# of about DISTANCE)
DISTANCE = 5000
MIN_ANGULAR_SIZE = math.radians(1 / DISTANCE)
OBSERVER_HEIGHT = 2


def calculate_moon_phase(moon_pos, planet_pos, star_pos):
    """
    Calculate the moon phase based on the positions of the moon, planet, and star.

    This function determines the moon phase by calculating the angle between the
    moon-planet vector and the planet-star vector. It then returns the corresponding
    moon phase as a string.

    Parameters:
    moon_pos (list or array-like): [x, y, z] coordinates of the moon in 3D space.
    planet_pos (list or array-like): [x, y, z] coordinates of the planet in 3D space.
    star_pos (list or array-like): [x, y, z] coordinates of the star in 3D space.

    Returns:
    str: The moon phase as a string, one of the following:
         "New Moon", "Waxing Crescent", "First Quarter", "Waxing Gibbous",
         "Full Moon", "Waning Gibbous", "Last Quarter", "Waning Crescent".
    """
    # Calculate the moon-planet and planet-star vectors
    moon_planet_vector = moon_pos - planet_pos
    planet_star_vector = planet_pos - star_pos

    # Calculate angle between vectors
    dot_product = np.dot(moon_planet_vector, planet_star_vector)
    moon_planet_distance = np.linalg.norm(moon_planet_vector)
    planet_star_distance = np.linalg.norm(planet_star_vector)
    angle_cos = dot_product / (moon_planet_distance * planet_star_distance)
    angle_rad = np.arccos(angle_cos)
    angle_deg = np.degrees(angle_rad)

    # Determine moon phase based on angle
    if angle_deg < 22.5 or angle_deg > 337.5:
        return "New Moon"
    elif 22.5 <= angle_deg < 67.5:
        return "Waxing Crescent"
    elif 67.5 <= angle_deg < 112.5:
        return "First Quarter"
    elif 112.5 <= angle_deg < 157.5:
        return "Waxing Gibbous"
    elif 157.5 <= angle_deg < 202.5:
        return "Full Moon"
    elif 202.5 <= angle_deg < 247.5:
        return "Waning Gibbous"
    elif 247.5 <= angle_deg < 292.5:
        return "Last Quarter"
    else:
        return "Waning Crescent"


def update_positions(current_time, star_system):
    """
    Update the positions of all celestial bodies and store them centrally.

    Parameters:
    current_time (float): The current time at which the positions should be calculated.
    star_system (dict): A dictionary containing CelestialBody objects representing a star system.

    Returns:
    dict: A dictionary containing the positions of all celestial bodies, with keys as body names
    and values as NumPy arrays of x, y, z coordinates.
    """
    positions = {}

    for body in star_system.values():
        position = calculate_position(body, current_time)
        positions[body.name] = position

    return positions



def planet_viewer(star_system ,current_time, positions, planet_name, lat, lon, skybox):
    """
    This function calculates and displays the visible celestial bodies in the sky, their altitude,
    azimuth, and angular size as seen by an observer on a specified planet in the solar system, at
    a given location and time.

    Parameters:
    current_time (float): The current time in hours since epoch (0-23.9999).
    positions (dict): A dictionary containing the positions of celestial bodies in the solar system.
    planet_name (str): The name of the planet on which the observer is located.
    lat (float): The latitude of the observer's location in degrees (-90 to 90).
    lon (float): The longitude of the observer's location in degrees (-180 to 180).

    Returns:
    numpy.ndarray: A structured array containing the visible celestial bodies, sorted by altitude in
    descending order. The array has the following fields:
    - 'name': The name of the celestial body (str).
    - 'altitude': The altitude of the celestial body in radians (float).
    - 'azimuth': The azimuth of the celestial body in radians (float).
    - 'angular_size': The angular size of the celestial body in radians (float).
    """
    planet_position = positions[planet_name]
    planet = star_system[planet_name]
    axial_tilt_rad = math.radians(planet.axial_tilt)

    observer_radius = km_to_au(planet.radius + OBSERVER_HEIGHT / 1000)
    observer_lat_rad = math.radians(lat)

    planet_rotation_angle = 2 * math.pi * (current_time % planet.rotation_period) / planet.rotation_period

    visible_bodies = []

    for body_name, body_position in positions.items():
        if body_name == planet_name:
            continue

        dx, dy, dz = body_position[0] - planet_position[0], body_position[1] - planet_position[1], body_position[2] - planet_position[2]
        distance_to_center = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

        distance = math.sqrt(observer_radius ** 2 + distance_to_center ** 2 - 2 * observer_radius * distance_to_center * math.cos(observer_lat_rad))

        declination = math.asin(dz / distance)
        right_ascension = math.atan2(dy, dx)

        lst = (planet_rotation_angle + lon) % (2 * math.pi)

        hour_angle = lst - right_ascension

        sin_latitude_tilted = math.sin(observer_lat_rad) * math.cos(axial_tilt_rad) - math.cos(observer_lat_rad) * math.sin(axial_tilt_rad) * math.cos(hour_angle)
        cos_latitude_tilted = math.cos(observer_lat_rad) * math.cos(axial_tilt_rad) + math.sin(observer_lat_rad) * math.sin(axial_tilt_rad) * math.cos(hour_angle)
        altitude = math.asin(sin_latitude_tilted * math.sin(declination) + cos_latitude_tilted * math.cos(declination) * math.cos(hour_angle))
        azimuth = math.atan2(-math.sin(hour_angle), cos_latitude_tilted * math.sin(declination) - sin_latitude_tilted * math.cos(declination) * math.cos(hour_angle))

        body = star_system[body_name]
        angular_size = math.atan2(body.radius, distance)

        if -math.pi / 2 <= altitude <= math.pi / 2 and angular_size >= MIN_ANGULAR_SIZE:
            visible_bodies.append((body_name, altitude, azimuth, angular_size))

    visible_bodies_np = np.array(visible_bodies, dtype=[('name', 'U20'), ('altitude', float), ('azimuth', float), ('angular_size', float)])
    visible_bodies_np.sort(order='altitude')
    visible_bodies_np = visible_bodies_np[::-1]

    for body_info in visible_bodies_np:
        print(f"{body_info['name']}: Altitude: {math.degrees(body_info['altitude']):.2f}°, Azimuth: {math.degrees(body_info['azimuth']):.2f}°, Angular Size: {math.degrees(body_info['angular_size']):.2f}°")

        body = star_system[body_info['name']]
        if body.body_type == "moon":
            moon_phase = calculate_moon_phase(positions[body_info['name']], positions[planet_name], positions["Sun"])
            print(f"{body_info['name']} is in the {moon_phase} phase")

    visible_bodies_np = np.array(visible_bodies, dtype=[('name', 'U20'), ('altitude', float), ('azimuth', float), ('angular_size', float)])
    visible_bodies_np.sort(order='altitude')
    visible_bodies_np = visible_bodies_np[::-1]

    visible_constellations = []

    # Loop through constellations in the skybox and calculate their altitude and azimuth
    for constellation_name, constellation in skybox.constellations.items():
        dx, dy, dz = Skybox.sphere_to_cartesian(skybox.radius, constellation.lat, constellation.long)

        distance_to_center = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        distance = math.sqrt(observer_radius ** 2 + distance_to_center ** 2 - 2 * observer_radius * distance_to_center * math.cos(observer_lat_rad))
        declination = math.asin(dz / distance)
        right_ascension = math.atan2(dy, dx)
        hour_angle = lst - right_ascension

        sin_latitude_tilted = math.sin(observer_lat_rad) * math.cos(axial_tilt_rad) - math.cos(observer_lat_rad) * math.sin(axial_tilt_rad) * math.cos(hour_angle)
        cos_latitude_tilted = math.cos(observer_lat_rad) * math.cos(axial_tilt_rad) + math.sin(observer_lat_rad) * math.sin(axial_tilt_rad) * math.cos(hour_angle)
        altitude = math.asin(sin_latitude_tilted * math.sin(declination) + cos_latitude_tilted * math.cos(declination) * math.cos(hour_angle))
        azimuth = math.atan2(-math.sin(hour_angle), cos_latitude_tilted * math.sin(declination) - sin_latitude_tilted * math.cos(declination) * math.cos(hour_angle))

        if -math.pi / 2 <= altitude <= math.pi / 2:
            visible_constellations.append((constellation_name, altitude, azimuth))

    visible_constellations_np = np.array(visible_constellations, dtype=[('name', 'U20'), ('altitude', float), ('azimuth', float)])
    visible_constellations_np.sort(order='altitude')
    visible_constellations_np = visible_constellations_np[::-1]

    for constellation_info in visible_constellations_np:
        print(f"{constellation_info['name']}: Altitude: {math.degrees(constellation_info['altitude']):.2f}°, Azimuth: {math.degrees(constellation_info['azimuth']):.2f}°")

    return visible_bodies_np, visible_constellations_np
