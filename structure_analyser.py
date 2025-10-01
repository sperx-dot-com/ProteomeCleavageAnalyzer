import os
import requests
import freesasa
from Bio import PDB
from Bio.PDB.DSSP import DSSP
import os
 
import pymol
from pymol import cmd

def fetch_alphafold_structure(uniprot_id, save_dir="structures"):
    """Fetch the latest AlphaFold structure for a given UniProt ID."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    os.makedirs(save_dir, exist_ok=True)
    pdb_path = os.path.join(save_dir, f"{uniprot_id}.pdb")
    
    # Check if the file already exists
    if os.path.exists(pdb_path):
        print(f"File for {uniprot_id} already exists at {pdb_path}. Skipping download.")
        return pdb_path
    
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(pdb_path, "wb") as f:
            f.write(response.content)
        return pdb_path
    else:
        raise ValueError(f"No Alphafold structure available for {uniprot_id}.")
    

def calculate_sasa(pdb_path, start, end):
    """Calculate Solvent Accessible Surface Area (SASA) for a defined region in a specific chain."""
    
    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)
    
    region = f"s1, resi {start}-{end}"
    sasa_of_region = freesasa.selectArea([region], structure, result)
    
    return sasa_of_region["s1"]

def get_secondary_structure(pdb_path, start, end, chain="A"):
    """Get secondary structure annotations for a given region in a specific chain using DSSP."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = next(structure.get_models())
    
    dssp = DSSP(model, pdb_path)
    sec_struct = {}
    for key in dssp.keys():
        if key[0] != chain:
            continue
        # Key format in DSSP is typically (chain, (het, resseq, icode)) or similar
        res_num = key[1][1] if isinstance(key[1], tuple) else key[1]
        if res_num in range(start, end + 1):
            sec_struct[res_num] = dssp[key][2]
    return sec_struct

def find_accession(fasta_file_path, gn_value):
    """Find UniProt accession number based on gene name in FASTA file."""
    if not os.path.exists(fasta_file_path):
        raise ValueError(f"FASTA file not found: {fasta_file_path}")
        
    with open(fasta_file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                # Split the header line into parts
                parts = line.strip().split('|')
                # The accession number is the second part
                accession = parts[1]
                
                # Check if the GN= value is present and matches the input
                gn_field = next((part for part in parts if 'GN=' in part), None)
                if gn_field:
                    # Extract the actual gene name from the GN field
                    gene_name = gn_field.split('GN=')[1].split()[0]
                    if gene_name == gn_value:
                        return accession

    return None  # Return None if no match is found


def count_residues_within_distance(pdb_file, chain, start, end, threshold):
    """
    Counts the number of residues (outside the target range)
    that have at least one atom within the given threshold (in Å)
    of any atom in the target residue range.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]  # Use the first model

    # Gather atoms from the target residue range
    target_atoms = []
    for res_id in range(start, end + 1):
        try:
            residue = model[chain][res_id]
            target_atoms.extend(list(residue.get_atoms()))
        except KeyError:
            continue

    interacting_residues = set()
    
    # Iterate over all residues in the same chain
    for residue in model[chain]:
        # Skip heteroatoms (e.g., water or ligands)
        if residue.id[0] != ' ':
            continue
        # Skip the target residues themselves
        if residue.id[1] in range(start, end + 1):
            continue
        
        # Check if any atom in this residue is within the threshold of any target atom
        for atom in residue.get_atoms():
            for target_atom in target_atoms:
                if (atom - target_atom) <= threshold:
                    interacting_residues.add(residue.id[1])
                    break  # No need to check more atoms for this residue
            else:
                continue
            break

    return len(interacting_residues)

def get_fasta_path_for_organism(organism_type):
    """Return the appropriate FASTA file path based on organism type."""
    if organism_type.lower() == "human":
        return 'uniprot_human.fa'
    elif organism_type.lower() == "ecoli":
        return 'uniprot_ecoli.fa'
    else:
        # Default to E. coli if unknown organism
        return 'uniprot_ecoli.fa'

def analyze_protein_region(gn_value, start, end, organism="ecoli", chain="A"):
    """
    Main function to analyze a protein region. Takes gene name, start and end positions
    and returns a dictionary with SASA, secondary structure, and residue counts within thresholds.
    
    Args:
        gn_value (str): Gene name
        start (int): Start residue position
        end (int): End residue position
        organism (str, optional): Organism type ("ecoli" or "human"). Defaults to "ecoli".
        chain (str, optional): Chain ID. Defaults to "A".
        
    Returns:
        dict: Dictionary containing analysis results
    """
    # Get the appropriate FASTA file path based on organism
    fasta_file_path = get_fasta_path_for_organism(organism)
    
    # Find UniProt ID from gene name
    uniprot_id = find_accession(fasta_file_path, gn_value)
    if not uniprot_id:
        raise ValueError(f"Could not find UniProt ID for gene name {gn_value} in {fasta_file_path}")
    
    result = {}
    uniprot = {"uniprot_id": uniprot_id}
    result.update(uniprot)
    

    # Fetch structure
    pdb_file = fetch_alphafold_structure(uniprot_id)
    
    # Calculate SASA
    sasa = calculate_sasa(pdb_file, start, end)
    
    # Get secondary structure
    sec_struct = get_secondary_structure(pdb_file, start, end, chain)
    
    # Calculate residues within thresholds
    thresholds = [5.0, 10.0, 15.0]
    residue_counts = {}
    for threshold in thresholds:
        count = count_residues_within_distance(pdb_file, chain, start, end, threshold)
        residue_counts[f"residues_within_{threshold}A"] = count
    
    # Create result dictionary
    structure_results = {
       # "uniprot_id": uniprot_id,
        "total_sasa": sasa,
        "secondary_structure": sec_struct
    }
    
    result.update(structure_results)

    # Add residue counts to result
    result.update(residue_counts)

    # Make pymol picture
    try:
        from pdb_visualizer import visualize_pdb_region
        visualize_pdb_region(pdb_file, start, end)
        print(f"PDB visualization generated for {gn_value} at {pdb_file}")
    except Exception as e:
        print(f"Could not generate PDB visualization: {e}")
    
    return result


# Example usage when running this file directly
if __name__ == "__main__":
    gn_value = "rpoD"
    start = 91
    end = 95
    
    result = analyze_protein_region(gn_value, start, end, organism="ecoli")
    print(result)

    """
    use in other 
    from protein_analysis import analyze_protein_region

    # Analyze a protein region
    result = analyze_protein_region("rpoD", 91, 95)
    print(result)

    # You can also specify a different organism if needed
    # result = analyze_protein_region("rpoD", 91, 95, organism="human")
    """