"""
graphical_visualisation.py

This script provides a function for creating a graphical representation of the visible
celestial bodies in the sky.
It takes as input a structured array of visible celestial bodies and a file name, and
creates a polar plot showing the location of each body in the sky. The plot is then saved
as an image file with the given file name.

Functions:
----------
plot_visible_bodies(visible_bodies, file_name): Generate a polar plot of the visible
celestial bodies and save the plot as an image file.
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_visible_bodies(visible_bodies, file_name):
    """
    Plot the visible celestial objects in the context of the viewer and save the plot as an
    image file.

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
        # Scale factor for circle size (this doesn't work as expected)
        size = ((body['angular_size'] / max_angular_size) * 5000) / 10
        ax.scatter(np.radians(angle), radius, s=size)
        ax.text(np.radians(angle), radius, body['name'], fontsize=10)

    # Save the plot
    plt.savefig(f"{file_name}.png")
    plt.close()
