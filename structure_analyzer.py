import os
from structure_analyser import (
    fetch_alphafold_structure,
    calculate_sasa,
    get_secondary_structure,
    find_accession,
    count_residues_within_distance,
    get_fasta_path_for_organism
)

class StructureAnalyzer:
    """Class to handle structure analysis business logic"""
    
    def __init__(self):
        self.structures_dir = "structures"
        os.makedirs(self.structures_dir, exist_ok=True)
    
    def analyze_proteins(self, hits, organism="ecoli",  progress_callback=None):
        """
        Analyze structure for a list of protein hits
        
        Args:
            hits (list): List of hit dictionaries
            organism (str): Organism type ("ecoli" or "human")
            
        Returns:
            list: List of hits with structure analysis added
        """
        analyzed_hits = []
        
        for i, hit in enumerate(hits):
            try:
                length = len(hit['sequence'][0])
                start = hit['occurrences'][0][0] + 1
                end = start + length
                gene_name = hit['gene_id']
                matched_sequence = str(hit['sequence'][0]) + str(hit['p1prime'][0])
                
                # Get the appropriate FASTA file path based on organism
                fasta_file_path = get_fasta_path_for_organism(organism)
                
                # Find UniProt ID from gene name
                uniprot_id = find_accession(fasta_file_path, hit["gene_id"])
                if not uniprot_id:
                    # Add hit with "Not found" status but still include uniprot_id
                    hit.update({
                        "uniprot_id": "Not found",
                        "total_sasa": "N/A",
                        "secondary_structure": {},
                        "residues_within_5.0A": "N/A",
                        "residues_within_10.0A": "N/A",
                        "residues_within_15.0A": "N/A"
                    })
                    analyzed_hits.append(hit)
                    continue
                
                # Try to fetch and analyze structure
                try:
                    pdb_file = fetch_alphafold_structure(uniprot_id, self.structures_dir)
                    
                    # Calculate SASA
                    sasa = calculate_sasa(pdb_file, start, end)
                    
                    # Get secondary structure
                    #sec_struct = get_secondary_structure(pdb_file, start, end)
                    sec_struct = []
                    # Calculate residues within thresholds
                    thresholds = [5.0, 10.0, 15.0]
                    residue_counts = {}
                    for threshold in thresholds:
                        count = count_residues_within_distance(pdb_file, "A", start, end, threshold)
                        residue_counts[f"residues_within_{threshold}A"] = count
                    
                    # Update hit with structure analysis
                    hit.update({
                        "uniprot_id": uniprot_id,
                        "total_sasa": sasa,
                        "secondary_structure": sec_struct,
                        **residue_counts
                    })
                        # Make pymol picture
                    try:
                        from pdb_visualizer import visualize_pdb_region
                        visualize_pdb_region(pdb_file, start, end, gene_name, matched_sequence)
                        print(f"PDB visualization generated for gene: {hit['gene_id']} (uniprot_id: {uniprot_id}) at {pdb_file}")
                    except Exception as e:
                        print(f"Could not generate PDB visualization: {e}")



                except Exception as e:
                    print(f"Error analyzing structure for {hit['gene_id']}: {str(e)}")
                    # If structure analysis fails, still include uniprot_id but mark other fields as N/A
                    hit.update({
                        "uniprot_id": uniprot_id,
                        "total_sasa": "N/A",
                        "secondary_structure": {},
                        "residues_within_5.0A": "N/A",
                        "residues_within_10.0A": "N/A",
                        "residues_within_15.0A": "N/A"
                    })
                

                analyzed_hits.append(hit)
                
            except Exception as e:
                print(f"Error processing hit {hit['gene_id']}: {str(e)}")
                # If any other error occurs, add hit with "Not found" status
                hit.update({
                    "uniprot_id": "Not found",
                    "total_sasa": "N/A",
                    "secondary_structure": {},
                    "residues_within_5.0A": "N/A",
                    "residues_within_10.0A": "N/A",
                    "residues_within_15.0A": "N/A"
                })
                analyzed_hits.append(hit)
        
                    # Call the progress callback if provided
        if progress_callback:
            progress_callback(i)
        return analyzed_hits

    def format_secondary_structure(self, sec_struct_dict):
        """Convert secondary structure dictionary to a string without commas."""
        if not sec_struct_dict:
            return "N/A"
        
        # Sort by position
        positions = sorted(sec_struct_dict.keys())
        return ''.join(sec_struct_dict[pos] for pos in positions) 