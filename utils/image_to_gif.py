import os
import imageio

# Set the folder containing the images
img_folder = '/home/dan/develop/natural-systems-mud/cosmology/plot/'

# Set the output GIF file name and frame rate
gif_name = 'planetary_movmenets.gif'
fps = 10

# Get a list of all the image file names in the folder
img_names = os.listdir(img_folder)
img_names.sort()

# Create a list of image paths
img_paths = [os.path.join(img_folder, img_name) for img_name in img_names]

# Create the GIF using imageio
with imageio.get_writer(gif_name, mode='I', fps=fps) as writer:
    for img_path in img_paths:
        img = imageio.imread(img_path)
        writer.append_data(img)