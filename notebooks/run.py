
# Import required library
import os
import sys  
sys.path.insert(0, '../src/')
from utils import *

# Import required classes
from IPython.display import display, Markdown

# Create an instance of PDBJsonParser and parse the JSON file
parser = PDBJsonParser('../rcsb_pdb_custom_report_20230624064914.json')

# Get unique database accessions from the parsed JSON
unique_accessions = parser.get_unique_database_accessions()

unique_sequences = {}
for accession, sequence in unique_accessions.items():
    # Save a dictionary with the accession as key and the sequence as value
    unique_sequences[accession] = sequence

import concurrent.futures

# Set the directory where the protein data will be stored
data_dir = "proteins_data"

def process_protein(uniprot_id, data):
    # Create an instance of the Protein class for each unique sequence
    protein = Protein(uniprot_id=uniprot_id, 
                      pdb_id=data['pdb_id'],
                      data_dir=data_dir, 
                      sequence_rcsb=data['sequence'])
    
    # Print the protein sequence
    sequence = protein.sequence_rcsb
    print(f"Protein {uniprot_id} Sequence: {sequence}")

    # Download the protein structure from AF2
    pdb_structure_af2 = protein.download_pdb_structures("AF2")

    # Download the protein structure from RCSB
    pdb_structure_rcsb = protein.download_pdb_structures("RCSB")

    # Download the af2 fasta file
    test = protein.sequence_uniprot()
    
    # Predict the structure of the protein using ESMFold
    # Create FASTA file path
    fasta_file = os.path.join(protein.protein_dir, f"{protein.uniprot_id}.fasta")

    # Create output PDB file path for ESM predicted structure
    esm_pdb_file = os.path.join(protein.protein_dir, f"{protein.uniprot_id}_ESM")

    # Get predicted structure using ESM
    protein.predict_structure_with_esm(fasta_file, esm_pdb_file)

    # Now you should have the ESM predicted structure saved in the path defined by `esm_pdb_file`
    print(f"ESM predicted structure saved at: {esm_pdb_file}")

# Process Q96524 protein

process_protein("Q96524", unique_sequences["Q96524"])