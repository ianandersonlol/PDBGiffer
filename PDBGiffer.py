import os
import sys
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBExceptions import PDBConstructionException
from Bio.PDB.PDBIO import PDBIO
from PIL import Image, ImageDraw
import math
import time
import warnings
warnings.filterwarnings("ignore") # This is really bad programming. Don't do it. Ian is a bad programmer.

def generate_gif(pdb_path, gif_path=None):
    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure('structure', pdb_path)
    # Build a list of polypeptides (amino acid chains)
    ppb = PPBuilder()
    polypeptides = ppb.build_peptides(structure)

    # Set up the image size and background color
    size = (400, 400)
    background = (255, 255, 255)

    # Create an empty image
    image = Image.new('RGB', size, background)
    draw = ImageDraw.Draw(image)

    # Set the number of frames in the GIF
    num_frames = 36

    # Create the GIF frames
    frames = []
    for i in range(num_frames):
        # Rotate the structure
        print(rotation_matrix(i * 360 / num_frames))
        print(i * 360 / num_frames)
        structure = structure.rotate_about(axis=(0, 0, 1), point=(0, 0, 0), angle=i * 360 / num_frames)

        # Create a PDBIO object and draw the structure
        io = PDBIO()
        io.set_structure(structure)
        io.save(draw, polypeptides=polypeptides)

        # Save the frame to the list of frames
        frame = image.copy()
        frames.append(frame)


    # Save the frames as a GIF animation
        if gif_path is not None:
            # If an output directory is specified, create the gif in the specified directory
            base_name = os.path.splitext(os.path.basename(pdb_path))[0]
            output_file = os.path.join(gif_path, base_name + '.gif')
        else:
            # If no output directory is specified, create the gif in the current working directory
            base_name = os.path.splitext(os.path.basename(pdb_path))[0]
            frames[0].save(base_name+".gif", format='GIF', append_images=frames[1:], save_all=True, duration=100, loop=0)

    
def rotation_matrix(angle):
    """Returns a rotation matrix for rotating a protein structure by a given angle."""
    angle = math.radians(angle)
    c = math.cos(angle)
    s = math.sin(angle)
    return [[c, -s, 0], [s, c, 0], [0, 0, 1]]

    
        
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
                    start_time = time.time()
                    generate_gif(pdb_file,output_path)
                    end_time = time.time()
                    print("Took "+str(end_time - start_time)+"s to generate a gif for "+pdb_file)
                    
                    