from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.pairwise2 import align, format_alignment
from Bio.SubsMat.MatrixInfo import blosum62
import requests
import os
import json


class A_thaliana_Protein:
    def __init__(self, uniprot_id):
        self.uniprot_id = uniprot_id
        self.sequence = self.download_sequence()
  
    def download_sequence(self):
        """Download sequence from uniprot and return sequence as string"""
        url = f'https://www.uniprot.org/uniprot/{self.uniprot_id}.fasta'
        response = requests.get(url)
        sequence = "".join(response.text.split("\n")[1:])
        return sequence


class AF2_Structure:
    def __init__(self, pdb_id, download_dir):
        self.pdb_id = pdb_id
        self.download_dir = download_dir
  
    def download_structure(self):
        """Download structure and save to download_dir"""
        url = f'https://alphafold.ebi.ac.uk/entry/{self.pdb_id}'
        response = requests.get(url) 
        filename = os.path.join(self.download_dir, f'{self.pdb_id}.pdb')
        
        with open(filename, 'w') as file:
            file.write(response.text)

        return filename




class ProteinCompleteness:
    def __init__(self, pdb_protein_sequence, uniprot_sequence):
        self.pdb_protein_sequence = pdb_protein_sequence
        self.uniprot_sequence = uniprot_sequence

    def completeness_check(self, completeness_threshold):
        pdb_length = len(self.pdb_protein_sequence)
        uniprot_length = len(self.uniprot_sequence)

        completeness = (pdb_length / uniprot_length) * 100
        
        if completeness >= completeness_threshold:
            return True
        else:
            return False


class SequenceSimilarity:
    def __init__(self, protein_sequence1, protein_sequence2):
        self.protein_sequence1 = protein_sequence1
        self.protein_sequence2 = protein_sequence2

    def calculate_similarity(self, similarity_threshold):
        alignment = align.globalds(self.protein_sequence1, self.protein_sequence2, blosum62, -10, -0.5)

        aligned_sequence1 = alignment[0][0]
        aligned_sequence2 = alignment[0][1]

        matches = sum([1 for i, j in zip(aligned_sequence1, aligned_sequence2) if i == j])
        total_length = max(len(self.protein_sequence1), len(self.protein_sequence2))

        similarity = (matches / total_length) * 100
        
        if similarity >= similarity_threshold:
            return True
        else:
            return False

# For ProteinCompleteness
#pro_com = ProteinCompleteness(pdb_protein_sequence, uniprot_sequence)
#print(f'Protein Completeness: {pro_com.completeness_check(50)}')

# For SequenceSimilarity
#seq_sim = SequenceSimilarity(protein_sequence1, protein_sequence2)
#print(f'Sequence Similarity: {seq_sim.calculate_similarity(99)}')

class PDBJsonParser:
    def __init__(self, json_file):
        self.json_file = json_file
        self.data = self.load_json()
        
    def load_json(self):
        with open(self.json_file) as file:
            data = json.load(file)
        return data

    def extract_info(self):
        result = []
        for protein in self.data:
            protein_info = {}
            protein_info['identifier'] = protein.get('identifier')
            protein_info['deposit_date'] = protein.get('data', {}).get('rcsb_accession_info', {}).get('deposit_date')
            protein_info['initial_release_date'] = protein.get('data', {}).get('rcsb_accession_info', {}).get('initial_release_date')
            protein_info['entry_id'] = protein.get('data', {}).get('rcsb_entry_container_identifiers', {}).get('entry_id')
            protein_info['title'] = protein.get('data', {}).get('struct', {}).get('title')
            protein_info['polymer_entities'] = []

            for poly in protein.get('data', {}).get('polymer_entities', []):
                pol_info = {}
                pol_info['sequence'] = poly.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can')
                pol_info['organism'] = poly.get('rcsb_entity_source_organism', [{}])[0].get('ncbi_scientific_name')
                pol_info['entity_id'] = poly.get('rcsb_polymer_entity_container_identifiers', {}).get('entity_id')
                reference_sequence = poly.get('rcsb_polymer_entity_container_identifiers', {}).get('reference_sequence_identifiers', [{}])
                pol_info['database_accession'] = reference_sequence[0].get('database_accession')
                pol_info['database_name'] = reference_sequence[0].get('database_name')
                protein_info['polymer_entities'].append(pol_info)

            result.append(protein_info)

        return result

# Usage
#parser = PDBJsonParser('path_to_your_json_file.json')
#parsed_data = parser.extract_info()
