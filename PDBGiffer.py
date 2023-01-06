import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import ArtistAnimation
plt.style.use('fivethirtyeight')

def generate_gif(input_file,output_dir=None):
    # read in the PDB file and extract the x, y, and z coordinates of the atoms
    with open(input_file, 'r') as f:
        lines = f.readlines()
    x, y, z = [], [], []
    for line in lines:
        if line[0:4] == 'ATOM':
            x.append(float(line[30:38]))
            y.append(float(line[38:46]))
            z.append(float(line[46:54]))

    # set up the figure and 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # create a list of angles for the gif
    angles = np.linspace(0, 360, 36)

    # create an empty list to store the plots
    plots = []

    # loop through the angles and create a plot for each angle
    for angle in angles:
        # convert the angle to radians
        angle_rad = np.deg2rad(angle)
        # rotate the x and y coordinates by the angle
        x_rotated = [x_i * np.cos(angle_rad) - y_i * np.sin(angle_rad) for x_i, y_i in zip(x, y)]
        y_rotated = [x_i * np.sin(angle_rad) + y_i * np.cos(angle_rad) for x_i, y_i in zip(x, y)]
        # create a scatter plot of the rotated coordinates
        plot = ax.scatter(x_rotated, y_rotated, z)
        # add the plot to the list of plots
        plots.append([plot])

    # create the gif using the list of plots
    anim = ArtistAnimation(fig, plots, interval=50, blit=True, repeat_delay=1000)

    if output_dir is not None:
        # If an output directory is specified, create the gif in the specified directory
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = os.path.join(output_dir, base_name + '.gif')
        anim.save(output_file, writer='imagemagick')
    else:
        # If no output directory is specified, create the gif in the current working directory
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        anim.save(base_name+'.gif', writer='imagemagick')
    
        
if __name__ == '__main__':
    if len(sys.argv) > 3:
        print("Too many arguments. PDBGiffer.py INPUT_PATH [OUTPUT_PATH]")
        exit()
        
    if len(sys.argv) < 2:
        print("Too few arguments. PDBGiffer.py INPUT_PATH [OUTPUT_PATH]")
        exit()
    if sys.argv[1].lower() == "help":
        print("Usage: generate_gif.py INPUT_FILE [OUTPUT_DIRECTORY]\n")
        print("Arguments:\n  INPUT_FILE       The input PDB file containing the protein structure.\n  [OUTPUT_DIRECTORY]  Optional directory to save the output gif. If not specified, the gif will be saved in the current working directory.\n")
        print("Description:\n  This script generates a gif of a spinning 3D animation of a protein from a PDB file. The gif is created by extracting the coordinates of the atoms in the protein from the PDB file and rendering them using the py3Dmol library. The gif will have 360 frames and will show the protein rotating 360 degrees around the y-axis. The atoms are displayed as sticks and are colored orange. The size of the atoms can be adjusted by setting the atom radius scale in the code.\n")
        print("Examples:\n  generate_gif.py protein.pdb\n  generate_gif.py protein.pdb /home/user/gifs")
        exit()    
    if sys.argv[2] and not os.path.isdir(sys.argv[2]):
        os.mkdir(sys.argv[2])
        print("Creating the directory: "+sys.argv[2])
    input_path = sys.argv[1]
    if sys.argv[2]:
        output_path = sys.argv[2]
    else:
        output_path = None
    if os.path.isfile(input_path):
        generate_gif(input_path,output_path)
    elif os.path.isdir(input_path):
        for root, dirs, files in os.walk(input_path):
            for file in files:
                if file.endswith('.pdb'):
                    pdb_file = os.path.join(root, file)
                    generate_gif(pdb_file,output_path)