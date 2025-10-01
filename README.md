# Proteome_Cleavage_Analyzer
A GUI tool for analyzing protein sequences in FASTA proteome files with integrated structure analysis and GO annotation capabilities.

## Features

- **Pattern Matching**: Search protein sequences using regular expressions
- **P1' Residue Filtering**: Filter hits based on specific P1' residues
- **C-Terminus Analysis**: Find patterns at or near the C-terminus
- **GO Annotation**: Integrate Gene Ontology annotations
- **Essential Gene Marking**: Identify essential genes
- **Structure Analysis**: Analyze protein structures including:
  - Secondary structure information
  - Solvent Accessible Surface Area (SASA)
  - Residue proximity analysis
- **Advanced Filtering**: Filter results using complex GO term expressions
- **Data Export**: Export results to Excel for further analysis

## Installation

### Prerequisites

- Python 3.7+
- Required Python packages (install via pip):
pip install tkinter openpyxl pandas biopython


### Setup

1. Clone this repository:
git clone https://github.com/yourusername/fasta-sequence-analyzer.git cd fasta-sequence-analyzer

2. Install dependencies:
pip install -r requirements.txt

3. Get proteome and GO annotation files:
- go to https://www.ncbi.nlm.nih.gov/datasets/genome/ 
- select the organism of interest
- select the proteome (preferable RefSeq)
- go to FTP
- download the proteome as ORGAINSM_translated_cds.faa
- download the GO annotation file as ORGAINSM_gene_ontology.gaf

4. Prepare your data directories:
- Place proteome FASTA files in the `genomes` folder
- Place GO annotation files in the `go_files` folder

## Usage

1. Run the application:
python gui_with_structure.py
2. Select a genome file or choose from the dropdown of available genomes
3. Enter a regex pattern to search for (e.g., `K[^P][^P]`)
4. Optionally specify P1' residues and other search parameters
5. Click "Search" to find matching sequences
6. Use the GO filter to refine results
7. Select entries and click "Structure Analysis" for detailed structural information
8. Export results to Excel for further analysis

## Example Searches

- **Caspase cleavage sites**: `D[EDST][^P][^P]`
- **Trypsin cleavage sites**: `[KR][^P]`
- **C-terminal motifs**: Check "At C-Terminus" and enter pattern like `[DE]EL$`

## GO Filter Examples

- `membrane` - Find all with 'membrane' in description
- `GO:0005886` - Find exact GO ID
- `membrane AND cytoplasm` - Require both terms
- `membrane OR transport` - Either term
- `membrane AND NOT nucleus` - Exclude term

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Gene Ontology Consortium for GO annotations
- UniProt for protein structure data
