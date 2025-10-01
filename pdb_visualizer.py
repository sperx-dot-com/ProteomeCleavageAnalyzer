import os
import sys

import pymol

from pymol import cmd
from pathlib import Path

def visualize_pdb_region(pdb_file, start_residue, end_residue, gene_name=None, matched_sequence=None, output_dir="./pymol_pictures", output_name=None):
    """
    Visualize a region of a protein structure from a PDB file using PyMOL.
    
    Parameters:
    -----------
    pdb_file : str
        Path to the PDB file
    start_residue : int
        Starting residue number for the highlighted region
    end_residue : int
        Ending residue number for the highlighted region
    output_dir : str, optional
        Directory to save the image (default: "pymol_pictures")
    output_name : str, optional
        Name of the output PNG file (default: same name as PDB file with .png extension)
    
    Returns:
    --------
    str
        Path to the generated PNG file
    """
   # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize PyMOL in true headless mode
    os.environ['PYMOL_QUIET'] = '1'
    pymol.pymol_argv = ['pymol', '-cQ']
    pymol.finish_launching()
    
    pdb_name = os.path.basename(pdb_file)
    uniprot_id = os.path.splitext(pdb_name)[0]

    # Determine output filename
    if output_name is None:
        # Use the same name as the PDB file but with .png extension
        pdb_name = os.path.basename(pdb_file)
        output_name_png = gene_name + "_" + os.path.splitext(pdb_name)[0] +  ".png"
        output_name_pse = gene_name + "_" + os.path.splitext(pdb_name)[0] + ".pse"

    output_path_png = os.path.join(output_dir, output_name_png)
    output_path_pse = os.path.join(output_dir, output_name_pse)

    
    # Load the PDB file
    cmd.load(pdb_file, "protein")
    
    # Show the whole protein as cartoon
    cmd.show_as("cartoon", "protein")
    
    # Define the region of interest
    region_selection = f"protein and resi {start_residue}-{end_residue}"
    cmd.select("region", region_selection)
    
    # Show the region in licorice representation
    cmd.show("licorice", "region")
    
    # Add labels with one-letter amino acid codes and numbers
    for i in range(start_residue, end_residue + 1):
        # cmd.label('''byca(region)''', 'oneletter+resi')
        cmd.label('''byca(region)''', 'oneletter')


    pymol.cmd.set('bg_rgb', '1.0, 1.0, 1.0')
    cmd.do('set label_bg_color, white')    
    cmd.do('set seq_view, 1')
    
    pymol.cmd.set('ray_opaque_background', 0)
    # Set some nice visualization parameters
    cmd.set("ray_shadows", "0")
    cmd.set("cartoon_fancy_helices", "1")
    cmd.set("cartoon_transparency", "0.1")
    cmd.set("label_size", "12")
    cmd.set("label_font_id", "11")  # Sans bold
    cmd.set("label_position", "(0,0,10)")

    #alphafold coroing
    cmd.do("run https://raw.githubusercontent.com/cbalbin-bio/pymol-color-alphafold/master/coloraf.py")
    cmd.do("coloraf protein")
    
    # Center the view on the region of interest
    #cmd.zoom("protein")
    cmd.orient("region", state=-1)
    cmd.zoom("protein", complete=1)

    #cmd.center("region")
    
    # Color the protein
    #cmd.color("marine", "protein")
    cmd.color_deep("red", 'region', 0)

    cmd.do('util.cnc("region",_self=cmd)')

    cmd.do('rotate y, 180, protein')

    # Ray trace and save the image
    cmd.ray(1200, 900)
    cmd.png(output_path_png, dpi=300)

    pymol_name= str(gene_name + "_" + uniprot_id)
    cmd.do(f"set_name protein, {pymol_name}")
    
    cmd.do(f"set_name region, {matched_sequence}")


    #Save the current session to a file
    cmd.save(output_path_pse)

    # Clean up
    cmd.delete("all")
    
    return output_path_png

class PDBVisualizer:
    def __init__(self, pdb_file=None, start_residue=None, end_residue=None):
        """
        Initialize the PDB visualizer.
        
        If all parameters are provided, automatically visualize the PDB file.
        
        Parameters:
        -----------
        pdb_file : str, optional
            Path to the PDB file
        start_residue : int, optional
            Starting residue number for the highlighted region
        end_residue : int, optional
            Ending residue number for the highlighted region
        """
        self.pdb_file = pdb_file
        self.start_residue = start_residue
        self.end_residue = end_residue
        
        # If all parameters are provided, run visualization
        if pdb_file and start_residue is not None and end_residue is not None:
            self.output_path = self.visualize()
    
    def visualize(self, pdb_file=None, start_residue=None, end_residue=None, output_dir="pymol_pictures"):
        """
        Visualize a region of a protein structure.
        
        Parameters:
        -----------
        pdb_file : str, optional
            Path to the PDB file (uses instance variable if not provided)
        start_residue : int, optional
            Starting residue number (uses instance variable if not provided)
        end_residue : int, optional
            Ending residue number (uses instance variable if not provided)
        output_dir : str, optional
            Directory to save the image (default: "pymol_pictures")
        
        Returns:
        --------
        str
            Path to the generated PNG file
        """
        # Use provided parameters or instance variables
        pdb_file = pdb_file or self.pdb_file
        start_residue = start_residue if start_residue is not None else self.start_residue
        end_residue = end_residue if end_residue is not None else self.end_residue
        
        # Validate parameters
        if not pdb_file:
            raise ValueError("PDB file path must be provided")
        if start_residue is None or end_residue is None:
            raise ValueError("Start and end residue numbers must be provided")
        
        # Update instance variables
        self.pdb_file = pdb_file
        self.start_residue = start_residue
        self.end_residue = end_residue
        
        # Call the visualization function
        return visualize_pdb_region(pdb_file, start_residue, end_residue, output_dir)

def main():
    """
    Main function for running the script as a standalone program.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Visualize a region of a protein structure from a PDB file.")
    parser.add_argument("pdb_file", help="Path to the PDB file")
    parser.add_argument("start_residue", type=int, help="Starting residue number for the highlighted region")
    parser.add_argument("end_residue", type=int, help="Ending residue number for the highlighted region")
    parser.add_argument("--output-dir", default="pymol_pictures", help="Directory to save the image (default: pymol_pictures)")
    parser.add_argument("--output-name", help="Name of the output PNG file (default: same name as PDB file with .png extension)")
    
    args = parser.parse_args()
    
    try:
        output_path = visualize_pdb_region(
            args.pdb_file, 
            args.start_residue, 
            args.end_residue, 
            args.output_dir,
            args.output_name
        )
        print(f"Image saved to: {output_path}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()