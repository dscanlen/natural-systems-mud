import math
import time
import random
# from evennia.utils import gametime
import matplotlib.pyplot as plt
from collections import namedtuple

# Constants
SECONDS_PER_DAY = 86400
SECONDS_PER_YEAR = 31536000
OBSERVER_HEIGHT = 2

# Minimum angular size for an object to be visible (1m3 object at a distance of about DISTANCE)
DISTANCE = 3000
MIN_ANGULAR_SIZE = math.radians(1 / DISTANCE)

# Celestial Body class
class CelestialBody:
    def __init__(self, name, orbit_radius, orbit_period, rotation_period, inclination, parent, radius):
        self.name = name
        self.orbit_radius = orbit_radius
        self.orbit_period = orbit_period
        self.rotation_period = rotation_period
        self.inclination = inclination
        self.parent = parent
        self.radius = radius


# Initialize solar system
def initialize_solar_system():
    SUN = CelestialBody("Sun", 0, 0, 0, 0, None, 696340)

    PLANETS = [
        CelestialBody("Mercury", 0.39, 87.97 * SECONDS_PER_DAY, 58.646 * SECONDS_PER_DAY, 7.00, SUN, 2440),
        CelestialBody("Venus", 0.72, 224.70 * SECONDS_PER_DAY, -243.018 * SECONDS_PER_DAY, 3.39, SUN, 6052),
        CelestialBody("Earth", 1.00, 365.26 * SECONDS_PER_DAY, 1 * SECONDS_PER_DAY, 0, SUN, 6371),
        CelestialBody("Jupiter", 5.20, 11.86 * SECONDS_PER_YEAR, 0.41354 * SECONDS_PER_DAY, 1.31, SUN, 69911),
        CelestialBody("Saturn", 9.58, 29.46 * SECONDS_PER_YEAR, 0.44401 * SECONDS_PER_DAY, 2.49, SUN, 58232),
    ]

    EARTH_MOONS = [
        CelestialBody("Moon1", 0.00257, 27 * SECONDS_PER_DAY, 27 * SECONDS_PER_DAY, 5.14, PLANETS[2], 1737),
        CelestialBody("Moon2", 0.005, 50 * SECONDS_PER_DAY, 50 * SECONDS_PER_DAY, 7.37, PLANETS[2], 1000),
    ]

    COMETS = [
        CelestialBody("Comet1", 35, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), SUN, 2),
        CelestialBody("Comet2", 50, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), SUN, 2),
        CelestialBody("Comet3", 80, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), SUN, 2),
        CelestialBody("Comet4", 100, 1 * SECONDS_PER_YEAR, 0, random.uniform(0, 180), SUN, 2),
    ]

    return [SUN] + PLANETS + EARTH_MOONS + COMETS

ALL_BODIES = initialize_solar_system()

# Helper functions
def calculate_position(body, time):
    """
    Calculate the x, y, z position of a celestial body at a given time.
    """
    if body.parent:
        parent_x, parent_y, parent_z = calculate_position(body.parent, time)
    else:
        parent_x, parent_y, parent_z = 0, 0, 0

    if body.orbit_period != 0:
        angle = 2 * math.pi * (time % body.orbit_period) / body.orbit_period
    else:
        angle = 0

    x = parent_x + body.orbit_radius * math.cos(angle)
    y = parent_y + body.orbit_radius * math.sin(angle)
    z = parent_z + body.orbit_radius * math.sin(body.inclination)

    return x, y, z


def update_positions():
    """
    Update the positions of all celestial bodies and store them centrally.
    """
    positions = {}
    # current_time = gametime.time()
    current_time = time.time()
    
    for body in ALL_BODIES:
        x, y, z = calculate_position(body, current_time)
        positions[body.name] = (x, y, z)
        print(f"{body.name} located at: {x}/{y}/{z}")
    
    return positions


def planet_viewer(positions, planet_name, lat, lon):
    """
    Display visible celestial bodies from the perspective of a viewer at a given latitude and longitude on a planet.
    """
    # Get the planet's position
    planet_position = positions[planet_name]
    planet_radius = [b.radius for b in ALL_BODIES if b.name == planet_name][0]

    visible_bodies = []

    for body_name, body_position in positions.items():
        if body_name == planet_name:
            continue

        # Calculate the relative x, y, and z coordinates of the celestial body
        dx, dy, dz = body_position[0] - planet_position[0], body_position[1] - planet_position[1], body_position[2] - planet_position[2]
        distance = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

        # Convert relative Cartesian coordinates to spherical coordinates (distance, inclination, azimuth)
        inclination = math.acos(dz / distance)
        azimuth = math.atan2(dy, dx)

        # Convert azimuth to a positive value in the range (0, 2π) or (0°, 360°)
        if azimuth < 0:
            azimuth += 2 * math.pi

        # Convert spherical coordinates to horizontal coordinates (altitude, azimuth)
        lat_rad, lon_rad = math.radians(lat), math.radians(lon)
        sin_altitude = (math.sin(inclination) * math.sin(lat_rad) + math.cos(inclination) * math.cos(lat_rad) * math.cos(azimuth - lon_rad)) * (planet_radius / (planet_radius + OBSERVER_HEIGHT))
        altitude = math.asin(sin_altitude)

        # Calculate angular size and check if it's greater than the minimum angular size
        body = [b for b in ALL_BODIES if b.name == body_name][0]
        angular_size = math.atan2(body.radius, distance)
        if -math.pi / 2 <= altitude <= math.pi / 2 and angular_size >= MIN_ANGULAR_SIZE:
            visible_bodies.append((body_name, altitude, azimuth, angular_size))

    # Sort visible bodies by altitude, and print them
    visible_bodies.sort(key=lambda x: x[1], reverse=True)
    for body_name, altitude, azimuth, angular_size in visible_bodies:
        print(f"{body_name}: Altitude: {math.degrees(altitude):.2f}°, Azimuth: {math.degrees(azimuth):.2f}°, Angular Size: {math.degrees(angular_size):.2f}°")

    # Sort visible bodies by altitude
    visible_bodies.sort(key=lambda x: x[1], reverse=True)
    return visible_bodies


def plot_visible_bodies(visible_bodies):
    """
    Plot the visible celestial objects in the context of the viewer.
    """
    # Create a polar plot
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

    # Set the plot's properties
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_rlim(0, 90)
    ax.set_rlabel_position(45)
    ax.set_yticklabels([])

    # Plot each visible object as a circle and add a label
    for body_name, altitude, azimuth, angular_size in visible_bodies:
        radius = 90 - math.degrees(altitude)
        angle = math.degrees(azimuth)
        size = math.degrees(angular_size)
        ax.scatter(math.radians(angle), radius, s=size)
        ax.text(math.radians(angle), radius, body_name, fontsize=10)

    # Display the plot
    plt.show()


def main():
    """
    Main function to run the astrological system.
    """
    positions = update_positions()
    visible_bodies = planet_viewer(positions, 'Earth', 0, 0)
    plot_visible_bodies(visible_bodies)

if __name__ == "__main__":
    main()
