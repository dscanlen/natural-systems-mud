import math
# import time
import random
# from evennia.utils import gametime
import matplotlib.pyplot as plt

# Constants
SECONDS_PER_DAY = 86400
SECONDS_PER_YEAR = 31536000
OBSERVER_HEIGHT = 2

# Minimum angular size for an object to be visible (1m3 object at a distance of about DISTANCE)
DISTANCE = 5000
MIN_ANGULAR_SIZE = math.radians(1 / DISTANCE)


# Celestial Body class
class CelestialBody:
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


# Initialize star system
def initialize_star_system():
    # Currently only supports single star-systems
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

# Helper functions
def km_to_au(km):
    return km / 149_597_870.7  # Divide the distance in km by the number of kilometers in one AU


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


def calculate_moon_phase(moon_pos, planet_pos, star_pos):
    moon_planet_vector = [moon_pos[0] - planet_pos[0], moon_pos[1] - planet_pos[1], moon_pos[2] - planet_pos[2]]
    planet_star_vector = [planet_pos[0] - star_pos[0], planet_pos[1] - star_pos[1], planet_pos[2] - star_pos[2]]

    # Calculate angle between vectors
    dot_product = sum(a * b for a, b in zip(moon_planet_vector, planet_star_vector))
    moon_planet_distance = math.sqrt(sum(a * a for a in moon_planet_vector))
    planet_star_distance = math.sqrt(sum(a * a for a in planet_star_vector))
    angle_cos = dot_product / (moon_planet_distance * planet_star_distance)
    angle_rad = math.acos(angle_cos)

    angle_deg = math.degrees(angle_rad)

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
    """
    positions = {}

    for body in STAR_SYSTEM:
        x, y, z = calculate_position(body, current_time)
        positions[body.name] = (x, y, z)
    
    return positions


def planet_viewer(current_time, positions, planet_name, lat, lon):
    planet_position = positions[planet_name]
    planet = [b for b in STAR_SYSTEM if b.name == planet_name][0]
    axial_tilt_rad = math.radians(planet.axial_tilt)

    observer_radius = km_to_au(planet.radius + OBSERVER_HEIGHT / 1000)  # incorporate observer height, converting it to km first
    observer_lat_rad = math.radians(lat)

    planet_rotation_angle = 2 * math.pi * (current_time % planet.rotation_period) / planet.rotation_period

    visible_bodies = []

    for body_name, body_position in positions.items():
        if body_name == planet_name:
            continue

        dx, dy, dz = body_position[0] - planet_position[0], body_position[1] - planet_position[1], body_position[2] - planet_position[2]
        distance_to_center = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

        # Adjust the distance by taking into account the observer's height
        distance = math.sqrt(observer_radius ** 2 + distance_to_center ** 2 - 2 * observer_radius * distance_to_center * math.cos(observer_lat_rad))

        # Calculate Right Ascension and Declination
        declination = math.asin(dz / distance)
        right_ascension = math.atan2(dy, dx)

        # Calculate Local Sidereal Time
        lst = (planet_rotation_angle + lon) % (2 * math.pi)

        # Calculate Hour Angle
        hour_angle = lst - right_ascension

        # Calculate Altitude and Azimuth (adjusted for axial tilt)
        sin_latitude_tilted = math.sin(observer_lat_rad) * math.cos(axial_tilt_rad) - math.cos(observer_lat_rad) * math.sin(axial_tilt_rad) * math.cos(hour_angle)
        cos_latitude_tilted = math.cos(observer_lat_rad) * math.cos(axial_tilt_rad) + math.sin(observer_lat_rad) * math.sin(axial_tilt_rad) * math.cos(hour_angle)
        altitude = math.asin(sin_latitude_tilted * math.sin(declination) + cos_latitude_tilted * math.cos(declination) * math.cos(hour_angle))
        azimuth = math.atan2(-math.sin(hour_angle), cos_latitude_tilted * math.sin(declination) - sin_latitude_tilted * math.cos(declination) * math.cos(hour_angle))

        body = [b for b in STAR_SYSTEM if b.name == body_name][0]
        angular_size = math.atan2(body.radius, distance)

        if -math.pi / 2 <= altitude <= math.pi / 2 and angular_size >= MIN_ANGULAR_SIZE:
            visible_bodies.append((body_name, altitude, azimuth, angular_size))

    visible_bodies.sort(key=lambda x: x[1], reverse=True)

    for body_name, altitude, azimuth, angular_size in visible_bodies:
        print(f"{body_name}: Altitude: {math.degrees(altitude):.2f}°, Azimuth: {math.degrees(azimuth):.2f}°, Angular Size: {math.degrees(angular_size):.2f}°")

        body = [b for b in STAR_SYSTEM if b.name == body_name][0]
        if body.body_type == "moon":
            moon_phase = calculate_moon_phase(positions[body_name], positions[planet_name], positions["Sun"])
            print(f"{body_name} is in the {moon_phase} phase")

    return visible_bodies
    


def plot_visible_bodies(visible_bodies, file_name):
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
    max_angular_size = max([angular_size for _, _, _, angular_size in visible_bodies])

    # Plot each visible object as a circle and add a label
    for body_name, altitude, azimuth, angular_size in visible_bodies:
        radius = 90 - math.degrees(altitude)
        angle = math.degrees(azimuth)
        size = ((angular_size / max_angular_size) * 5000) / 10 # Scale factor for circle size
        ax.scatter(math.radians(angle), radius, s=size)
        ax.text(math.radians(angle), radius, body_name, fontsize=10)

    # Display the plot
    # plt.show()
    plt.savefig(f"{file_name}.png")
    plt.close()



def main():
    """
    Main function to run the astrological system.
    """
    # current_time = gametime.time()
    # current_time = time.time()
    current_time = 45527321.453
    filename = 1
    for i in range(1, 400):
        positions = update_positions(current_time)
        visible_bodies = planet_viewer(current_time, positions, 'Earth', 0, 0)
        plot_visible_bodies(visible_bodies, str(filename).zfill(3))  # add leading zeros to filename
        filename += 1
        current_time += 1800

if __name__ == "__main__":
    main()
