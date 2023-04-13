import math
from evennia.utils import gametime
from collections import namedtuple

# Define celestial body
CelestialBody = namedtuple("CelestialBody", "name orbit_radius orbit_period rotation_period inclination parent radius")

# Constants
SECONDS_PER_DAY = 86400
SECONDS_PER_YEAR = 31536000

# Minimum angular size for an object to be visible (1m3 object at a distance of about DISTANCE)
DISTANCE = 3000
MIN_ANGULAR_SIZE = math.radians(1 / DISTANCE)

# Solar system data
SUN = CelestialBody("Sun", 0, 0, 0, 0, None, 696340)

PLANETS = [
    CelestialBody("Mercury", 0.39, 87.97 * SECONDS_PER_DAY, 58.646 * SECONDS_PER_DAY, 0, SUN, 2440),
    CelestialBody("Venus", 0.72, 224.70 * SECONDS_PER_DAY, -243.018 * SECONDS_PER_DAY, 0, SUN, 6052),
    CelestialBody("Earth", 1.00, 365.26 * SECONDS_PER_DAY, 1 * SECONDS_PER_DAY, 0, SUN, 6371),
    CelestialBody("Jupiter", 5.20, 11.86 * SECONDS_PER_YEAR, 0.41354 * SECONDS_PER_DAY, 0, SUN, 69911),
    CelestialBody("Saturn", 9.58, 29.46 * SECONDS_PER_YEAR, 0.44401 * SECONDS_PER_DAY, 0, SUN, 58232),
]

EARTH_MOONS = [
    CelestialBody("Moon1", 0.00257, 27 * SECONDS_PER_DAY, 27 * SECONDS_PER_DAY, 0, PLANETS[2], 1737),
    CelestialBody("Moon2", 0.005, 50 * SECONDS_PER_DAY, 50 * SECONDS_PER_DAY, 0, PLANETS[2], 1000),
]

COMETS = [
    CelestialBody("Comet1", 35, 1 * SECONDS_PER_YEAR, 0, 0, SUN, 2),
    CelestialBody("Comet2", 50, 1 * SECONDS_PER_YEAR, 0, 0, SUN, 2),
    CelestialBody("Comet3", 80, 1 * SECONDS_PER_YEAR, 0, 0, SUN, 2),
    CelestialBody("Comet4", 100, 1 * SECONDS_PER_YEAR, 0, 0, SUN, 2),
]

ALL_BODIES = [SUN] + PLANETS + EARTH_MOONS + COMETS

# Helper functions
def calculate_position(body, time):
    """
    Calculate the x, y, z position of a celestial body at a given time.
    """
    if body.parent:
        parent_x, parent_y, parent_z = calculate_position(body.parent, time)
    else:
        parent_x, parent_y, parent_z = 0, 0, 0

    angle = 2 * math.pi * (time % body.orbit_period) / body.orbit_period
    x = parent_x + body.orbit_radius * math.cos(angle)
    y = parent_y + body.orbit_radius * math.sin(angle)
    z = parent_z + body.orbit_radius * math.sin(body.inclination)

    return x, y, z

def update_positions():
    """
    Update the positions of all celestial bodies and store them centrally.
    """
    positions = {}
    current_time = gametime.time()
    
    for body in ALL_BODIES:
        x, y, z = calculate_position(body, current_time)
        positions[body.name] = (x, y, z)
    
    return positions

def planet_viewer(positions, planet_name, lat, lon):
    """
    Display celestial objects visible from a specified planet at a given latitude and longitude.
    """
    planet_position = positions[planet_name]
    visible_bodies = []

    for body_name, body_position in positions.items():
        body = next(body for body in ALL_BODIES if body.name == body_name)
        if body_name == planet_name:
            continue

        # Calculate the distance and position difference between the specified planet and the celestial body
        dx, dy, dz = body_position[0] - planet_position[0], body_position[1] - planet_position[1], body_position[2] - planet_position[2]
        distance = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

        # Convert latitude, longitude, and distance to Cartesian coordinates
        lat_rad, lon_rad = math.radians(lat), math.radians(lon)
        x_rel = distance * math.cos(lat_rad) * math.cos(lon_rad)
        y_rel = distance * math.cos(lat_rad) * math.sin(lon_rad)
        z_rel = distance * math.sin(lat_rad)

        # Calculate altitude (angle above horizon) and azimuth (angle along horizon)
        altitude = math.asin(z_rel / distance)
        azimuth = math.atan2(y_rel, x_rel)

        # Calculate angular size and check if it's greater than the minimum angular size
        angular_size = math.atan2(body.radius, distance)
        if -math.pi / 2 <= altitude <= math.pi / 2 and angular_size >= MIN_ANGULAR_SIZE:
            visible_bodies.append((body_name, altitude, azimuth, angular_size))

    # Sort visible bodies by altitude, and print them
    visible_bodies.sort(key=lambda x: x[1], reverse=True)
    for body_name, altitude, azimuth, angular_size in visible_bodies:
        print(f"{body_name}: Altitude: {math.degrees(altitude):.2f}°, Azimuth: {math.degrees(azimuth):.2f}°, Angular Size: {math.degrees(angular_size):.2f}°")

def main():
    """
    Main function to run the astrological system.
    """
    positions = update_positions()
    planet_viewer(positions, "Earth", 0, 0)

if __name__ == "__main__":
    main()
