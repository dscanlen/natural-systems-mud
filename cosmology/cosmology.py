import os
import math
# import time
import random
import numpy as np
from evennia.utils import gametime
# run the processes multithreaded? probably not needed as it will only be used one calculation at a time. 
# potentially moon phase will be useful to run it their own threads

import matplotlib.pyplot as plt

# Constants
SECONDS_PER_DAY = 86400
SECONDS_PER_YEAR = 31536000
OBSERVER_HEIGHT = 2

# Minimum angular size for an object to be visible (1m3 object at a distance of about DISTANCE)
DISTANCE = 5000
MIN_ANGULAR_SIZE = math.radians(1 / DISTANCE)


class CelestialBody:
    """
    A class representing a celestial body in a star system.

    Attributes:
    -----------
    name : str
        The name of the celestial body.
    body_type : str
        The type of celestial body (e.g., 'star', 'planet', 'moon', 'comet').
    orbit_radius : float
        The average distance from the celestial body to its parent body in astronomical units (AU).
    orbit_period : float
        The time it takes for the celestial body to complete one orbit around its parent body in seconds.
    rotation_period : float
        The time it takes for the celestial body to complete one rotation around its axis in seconds.
    inclination : float
        The angle between the celestial body's orbital plane and the reference plane in degrees.
    parent : CelestialBody or None
        The parent celestial body that this body orbits, or None if this body is the central star.
    radius : float
        The radius of the celestial body in kilometers.
    axial_tilt : float
        The angle between the celestial body's rotation axis and the perpendicular to its orbital plane in degrees.
    """

    def __init__(self, name, body_type, orbit_radius, orbit_period, rotation_period, inclination, parent, radius, axial_tilt):
        self.name = name
        self.body_type = body_type
        self.orbit_radius = orbit_radius
        self.orbit_period = orbit_period
        self.rotation_period = rotation_period
        self.inclination = inclination
        self.parent = parent
        self.radius = radius
        self.axial_tilt = axial_tilt


def initialize_star_system():    
    """
    Initialize a star system with celestial bodies.

    This version currently supports only single star-systems.

    Returns:
    --------
    list
        A list of CelestialBody instances, including the central star, planets, moons, and comets in the star system.
    """

    SUN = CelestialBody("Sun", "star", 0, 0, 0, 0, None, 696340, 0)

    PLANETS = [
        CelestialBody("Mercury", "planet", 0.39, 87.97 * SECONDS_PER_DAY, 58.646 * SECONDS_PER_DAY, 7.00, SUN, 2440, 0.034),
        CelestialBody("Venus", "planet", 0.72, 224.70 * SECONDS_PER_DAY, -243.018 * SECONDS_PER_DAY, 3.39, SUN, 6052, 2.64),
        CelestialBody("Earth", "planet", 1.00, 365.26 * SECONDS_PER_DAY, 1 * SECONDS_PER_DAY, 0, SUN, 6371, 23.44),
        CelestialBody("Jupiter", "planet", 5.20, 11.86 * SECONDS_PER_YEAR, 0.41354 * SECONDS_PER_DAY, 1.31, SUN, 69911, 3.13),
        CelestialBody("Saturn", "planet", 9.58, 29.46 * SECONDS_PER_YEAR, 0.44401 * SECONDS_PER_DAY, 2.49, SUN, 58232, 26.73),
    ]

    MOONS = [
        CelestialBody("Moon1", "moon", 0.00257, 27 * SECONDS_PER_DAY, 27 * SECONDS_PER_DAY, 5.14, PLANETS[2], 1737, 0),
        CelestialBody("Moon2", "moon", 0.005, 50 * SECONDS_PER_DAY, 50 * SECONDS_PER_DAY, 7.37, PLANETS[2], 1000, 0),
    ]

    COMETS = [
        CelestialBody("Comet1", "comet", 35, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), SUN, 6, 0),
        CelestialBody("Comet2", "comet", 50, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), SUN, 10, 0),
        CelestialBody("Comet3", "comet", 80, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), SUN, 25, 0),
        CelestialBody("Comet4", "comet", 100, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), SUN, 92, 0),
    ]

    return [SUN] + PLANETS + MOONS + COMETS

STAR_SYSTEM = initialize_star_system()
STAR_SYSTEM_DICT = {body.name: body for body in STAR_SYSTEM}


def km_to_au(km):
    """
    Convert a distance in kilometers to astronomical units (AU).

    Parameters:
    km : float
        The distance in kilometers to be converted to astronomical units.

    Returns:
    float
        The distance in astronomical units.
    """
    # Divide the distance in km by the number of kilometers in one AU (approx. 149,597,870.7 km)
    return km / 149_597_870.7


def calculate_position(body, time):
    """
    Calculate the x, y, z position of a celestial body at a given time.
    
    This function recursively calculates the position of a celestial body in 3D space
    at a given time, taking into account the body's parent (if any), orbit radius,
    orbit period, and inclination. It returns the x, y, z coordinates of the body's
    position as a numpy array.
    
    Parameters:
    body (object): A celestial body object with the following attributes:
                   - parent (object or None): The parent body, if any, or None.
                   - orbit_radius (float): The radius of the body's orbit.
                   - orbit_period (float): The period of the body's orbit.
                   - inclination (float): The inclination of the body's orbit.
    time (float): The time at which to calculate the body's position, in arbitrary units.

    Returns:
    np.ndarray: The x, y, z coordinates of the celestial body's position as a numpy array.
    """
    # Calculate the position of the parent body, if any
    if body.parent:
        parent_x, parent_y, parent_z = calculate_position(body.parent, time)
    else:
        parent_x, parent_y, parent_z = 0, 0, 0

    # Calculate the angle of the body in its orbit at the given time
    if body.orbit_period != 0:
        angle = 2 * np.pi * (time % body.orbit_period) / body.orbit_period
    else:
        angle = 0

    # Calculate the x, y, z position of the body relative to its parent
    x = parent_x + body.orbit_radius * np.cos(angle)
    y = parent_y + body.orbit_radius * np.sin(angle)
    z = parent_z + body.orbit_radius * np.sin(body.inclination)

    return np.array([x, y, z])


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
    

def update_positions(current_time):
    """
    Update the positions of all celestial bodies and store them centrally.
    
    Parameters:
    current_time (float): The current time at which the positions should be calculated.
    
    Returns:
    dict: A dictionary containing the positions of all celestial bodies, with keys as body names and values as NumPy arrays of x, y, z coordinates.
    """
    positions = {}

    for body in STAR_SYSTEM:
        position = calculate_position(body, current_time)
        positions[body.name] = position
    
    return positions


def planet_viewer(current_time, positions, planet_name, lat, lon):
    """
    This function calculates and displays the visible celestial bodies in the sky, their altitude, azimuth, and angular size as seen by an observer on a specified planet in the solar system, at a given location and time.

    Parameters:
    current_time (float): The current time in hours since epoch (0-23.9999).
    positions (dict): A dictionary containing the positions of celestial bodies in the solar system.
    planet_name (str): The name of the planet on which the observer is located.
    lat (float): The latitude of the observer's location in degrees (-90 to 90).
    lon (float): The longitude of the observer's location in degrees (-180 to 180).

    Returns:
    numpy.ndarray: A structured array containing the visible celestial bodies, sorted by altitude in descending order. The array has the following fields:
    - 'name': The name of the celestial body (str).
    - 'altitude': The altitude of the celestial body in radians (float).
    - 'azimuth': The azimuth of the celestial body in radians (float).
    - 'angular_size': The angular size of the celestial body in radians (float).
    """
    planet_position = positions[planet_name]
    planet = STAR_SYSTEM_DICT[planet_name]
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

        body = STAR_SYSTEM_DICT[body_name]
        angular_size = math.atan2(body.radius, distance)

        if -math.pi / 2 <= altitude <= math.pi / 2 and angular_size >= MIN_ANGULAR_SIZE:
            visible_bodies.append((body_name, altitude, azimuth, angular_size))

    visible_bodies_np = np.array(visible_bodies, dtype=[('name', 'U20'), ('altitude', float), ('azimuth', float), ('angular_size', float)])
    visible_bodies_np.sort(order='altitude')
    visible_bodies_np = visible_bodies_np[::-1]

    for body_info in visible_bodies_np:
        print(f"{body_info['name']}: Altitude: {math.degrees(body_info['altitude']):.2f}°, Azimuth: {math.degrees(body_info['azimuth']):.2f}°, Angular Size: {math.degrees(body_info['angular_size']):.2f}°")

        body = STAR_SYSTEM_DICT[body_info['name']]
        if body.body_type == "moon":
            moon_phase = calculate_moon_phase(positions[body_info['name']], positions[planet_name], positions["Sun"])
            print(f"{body_info['name']} is in the {moon_phase} phase")

    return visible_bodies_np


def plot_visible_bodies(visible_bodies, file_name):
    """
    Plot the visible celestial objects in the context of the viewer and save the plot as an image file.

    Parameters
    ----------
    visible_bodies : numpy.ndarray
        A structured NumPy array containing information about the visible celestial bodies.
        The fields are 'name', 'altitude', 'azimuth', and 'angular_size'.
    file_name : str
        The name of the file to save the plot as (without the file extension).

    Returns
    -------
    None
    """
    # Create a polar plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    # Set the plot's properties
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rlim(0, 90)
    ax.set_rlabel_position(45)
    ax.set_yticklabels([])

    # Find the maximum angular size
    max_angular_size = np.amax(visible_bodies['angular_size'])

    # Plot each visible object as a circle and add a label
    for body in visible_bodies:
        radius = 90 - np.degrees(body['altitude'])
        angle = np.degrees(body['azimuth'])
        size = ((body['angular_size'] / max_angular_size) * 5000) / 10  # Scale factor for circle size
        ax.scatter(np.radians(angle), radius, s=size)
        ax.text(np.radians(angle), radius, body['name'], fontsize=10)

    # Save the plot
    plt.savefig(f"{file_name}.png")
    plt.close()


def main():
    """
    Main function to run the astrological system.
    """
    current_time = gametime.time()
    # current_time = time.time()
    # current_time = 45527321.453
    filename = 1
    os.chdir("/home/dan/develop/natural-systems-mud/cosmology/plot/")
    for i in range(1, 400):
        positions = update_positions(current_time)
        visible_bodies = planet_viewer(current_time, positions, 'Earth', 0, 0)
        plot_visible_bodies(visible_bodies, str(filename).zfill(3))  # add leading zeros to filename
        filename += 1
        current_time += 1800


if __name__ == "__main__":
    main()
