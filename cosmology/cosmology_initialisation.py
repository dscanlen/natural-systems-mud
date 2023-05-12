"""
cosmology_initialisation.py

This module contains the classes and functions necessary to initialize a star system. It includes
the creation of celestial bodies (such as stars, planets, moons, and comets), the initialization of
star systems, and the definition of constants related to time and space.

This module is structured as follows:

1. Constants: Definitions of various constants, like the number of seconds per day and per year,
    and observer height.

2. Classes: Definitions of classes to represent celestial bodies (Star, Constellation,
    Skybox, CelestialBody).

3. Initialization Function: A function to initialize a star system with celestial bodies.

"""

import random

# Constants
SECONDS_PER_DAY = 86400
SECONDS_PER_YEAR = 31536000

class Star:
    """
    A class to represent a star.

    Attributes:
    name (str): the name of the star.
    coordinates (tuple): a tuple representing the star's right ascension, declination, and distance.
    """

    def __init__(self, name, coordinates):
        """
        Constructs all the necessary attributes for the Star object.

        Parameters:
        name (str): the name of the star.
        coordinates (tuple): a tuple representing the star's right ascension, declination,
            and distance.
        """
        self.name = name
        self.coordinates = coordinates


class Constellation:
    """
    A class to represent a constellation.

    Attributes:
    name (str): the name of the constellation.
    stars (list): a list of Star objects that make up the constellation.
    """

    def __init__(self, name, stars):
        """
        Constructs all the necessary attributes for the Constellation object.

        Parameters:
        name (str): the name of the constellation.
        stars (list): a list of Star objects that make up the constellation.
        """
        self.name = name
        self.stars = stars


class Skybox:
    """
    A class to represent the skybox as a sphere surrounding the star system.

    Attributes:
    radius (float): the radius of the sphere.
    stars (dict): a dictionary mapping coordinates to Star objects.
    constellations (dict): a dictionary mapping coordinates to Constellation objects.
    """

    def __init__(self, radius, stars=None, constellations=None):
        """
        Constructs all the necessary attributes for the Skybox object.

        Parameters:
        radius (float): the radius of the sphere.
        stars (dict): a dictionary mapping coordinates to Star objects. Default is an empty dict.
        constellations (dict): a dictionary mapping coordinates to Constellation objects.
                Default is an empty dict.
        """
        self.radius = radius
        self.stars = stars if stars is not None else {}
        self.constellations = constellations if constellations is not None else {}


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
        The time it takes for the celestial body to complete one orbit around its parent
        body in seconds.
    rotation_period : float
        The time it takes for the celestial body to complete one rotation around its axis
        in seconds.
    inclination : float
        The angle between the celestial body's orbital plane and the reference plane in degrees.
    parent : CelestialBody or None
        The parent celestial body that this body orbits, or None if this body is the central star.
    radius : float
        The radius of the celestial body in kilometers.
    axial_tilt : float
        The angle between the celestial body's rotation axis and the perpendicular to its orbital
        plane in degrees.
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
        A list of CelestialBody instances, including the central star, planets, moons, and comets in
        the star system.
    """

    sun = CelestialBody("Sun", "star", 0, 0, 0, 0, None, 696340, 0)

    planets = [
        CelestialBody("Mercury", "planet", 0.39, 87.97 * SECONDS_PER_DAY, 58.646 * SECONDS_PER_DAY, 7.00, sun, 2440, 0.034),
        CelestialBody("Venus", "planet", 0.72, 224.70 * SECONDS_PER_DAY, -243.018 * SECONDS_PER_DAY, 3.39, sun, 6052, 2.64),
        CelestialBody("Earth", "planet", 1.00, 365.26 * SECONDS_PER_DAY, 1 * SECONDS_PER_DAY, 0, sun, 6371, 23.44),
        CelestialBody("Jupiter", "planet", 5.20, 11.86 * SECONDS_PER_YEAR, 0.41354 * SECONDS_PER_DAY, 1.31, sun, 69911, 3.13),
        CelestialBody("Saturn", "planet", 9.58, 29.46 * SECONDS_PER_YEAR, 0.44401 * SECONDS_PER_DAY, 2.49, sun, 58232, 26.73),
    ]

    moons = [
        CelestialBody("Moon1", "moon", 0.00257, 27 * SECONDS_PER_DAY, 27 * SECONDS_PER_DAY, 5.14, planets[2], 1737, 0),
        CelestialBody("Moon2", "moon", 0.005, 50 * SECONDS_PER_DAY, 50 * SECONDS_PER_DAY, 7.37, planets[2], 1000, 0),
    ]

    comets = [
        CelestialBody("Comet1", "comet", 35, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), sun, 6, 0),
        CelestialBody("Comet2", "comet", 50, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), sun, 10, 0),
        CelestialBody("Comet3", "comet", 80, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), sun, 25, 0),
        CelestialBody("Comet4", "comet", 100, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), sun, 92, 0),
    ]

    star_system = [sun] + planets + moons + comets

    return {body.name: body for body in star_system}


def initialize_skybox(radius):
    """
    Initialize a skybox for the star system.

    This function creates a Skybox object with a specified radius. The stars and constellations
    are populated with example stars and constellations.

    Parameters:
    -----------
    radius : float
        The radius of the skybox. This is the distance from the center of the star system to the
        edge of the skybox.

    Returns:
    --------
    Skybox
        A Skybox object with the specified radius and populated stars and constellations
        dictionaries.
    """

    # Initialize some stars with names and coordinates
    stars = {
        "Alpha Centauri": Star("Alpha Centauri", (14.63, -60.83, 4.367)),
        "Sirius": Star("Sirius", (6.75, -16.716, 8.611)),
        "Betelgeuse": Star("Betelgeuse", (5.919, 7.407, 643)),
    }

    # Initialize a constellation with the name and list of stars
    constellations = {
        "Canis Major": Constellation("Canis Major", [stars["Sirius"]]),
        "Orion": Constellation("Orion", [stars["Betelgeuse"]]),
    }

    return Skybox(radius, stars, constellations)