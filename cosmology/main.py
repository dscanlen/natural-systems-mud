"""
main.py
This script serves as the main entry point for running the astronomical system. It initializes the
star system, calculates the positions of celestial bodies at a given time, and produces a graphical
representation of the visible bodies from a viewer's perspective on Earth.

Functions:
----------
main(): Initialize the star system, update positions of celestial bodies, create graphical
representations, and save the images to a specified directory.
"""

import os
import time
# from evennia import gametime
from cosmology_initialisation import initialize_star_system, initialize_skybox
from position_calculations import compute_max_apogee
from viewer_visualisation import update_positions, planet_viewer
from graphical_visualisation import plot_visible_bodies

def main():
    """
    Main function to run the astrological system.
    """
    # initialise everything:
    star_system = initialize_star_system()
    max_apogee = compute_max_apogee(star_system)
    sky_box = initialize_skybox(max_apogee + 10)

    # get the time:
    # current_time = gametime.time()
    current_time = time.time()
    # current_time = 45527321.453

    # setup sky plotting
    filename = 1
    os.chdir(os.getcwd() + "/cosmology/plot/")

    # test run
    for i in range(1, 400):
        positions = update_positions(current_time, star_system)
        visible_bodies, visible_constellations = planet_viewer(star_system, current_time, positions, 'Earth', 40, -74, sky_box)
        plot_visible_bodies(visible_bodies, visible_constellations, str(filename).zfill(3))  # add leading zeros to filename
        filename += 1
        current_time += 1800


if __name__ == "__main__":
    main()
