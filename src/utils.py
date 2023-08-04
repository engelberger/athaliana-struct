from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.pairwise2 import align, format_alignment
from Bio.SubsMat.MatrixInfo import blosum62
import requests
import os
import json
from collections import Counter


# test


class A_thaliana_Protein:
    def __init__(self, uniprot_id):
        self.uniprot_id = uniprot_id
        self.sequence = self.download_sequence()

    def download_sequence(self):
        """Download sequence from uniprot and return sequence as string"""
        url = f"https://www.uniprot.org/uniprot/{self.uniprot_id}.fasta"
        response = requests.get(url)
        sequence = "".join(response.text.split("\n")[1:])
        return sequence


class AF2_Structure:
    def __init__(self, pdb_id, download_dir):
        self.pdb_id = pdb_id
        self.download_dir = download_dir

    def download_structure(self):
        """Download structure and save to download_dir"""
        # Example https://alphafold.ebi.ac.uk/files/AF-Q8W3K0-F1-model_v4.pdb
        url = f"https://alphafold.ebi.ac.uk/files/AF-{self.pdb_id}-F1-model_v4.pdb"
        response = requests.get(url)
        filename = os.path.join(self.download_dir, f"{self.pdb_id}.pdb")

        with open(filename, "w") as file:
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
        alignment = align.globalds(
            self.protein_sequence1, self.protein_sequence2, blosum62, -10, -0.5
        )

        aligned_sequence1 = alignment[0][0]
        aligned_sequence2 = alignment[0][1]

        matches = sum(
            [1 for i, j in zip(aligned_sequence1, aligned_sequence2) if i == j]
        )
        total_length = max(len(self.protein_sequence1), len(self.protein_sequence2))

        similarity = (matches / total_length) * 100

        if similarity >= similarity_threshold:
            return True
        else:
            return False


class PDBJsonParser:
    def __init__(self, json_file):
        self.json_file = json_file

    def parse_json(self):
        with open(self.json_file) as file:
            data = json.load(file)
        result = []
        for entry in data:
            pdb_id = entry.get("identifier", "Not provided")
            rcsb_id = entry.get("data", {}).get("rcsb_id", "Not provided")
            rcsb_accession = (
                entry.get("data", {})
                .get("rcsb_accession_info", {})
                .get("deposit_date", "Not provided")
            )
            initial_release = (
                entry.get("data", {})
                .get("rcsb_accession_info", {})
                .get("initial_release_date", "Not provided")
            )
            title = entry.get("data", {}).get("struct", {}).get("title", "Not provided")

            polymer_entities = entry.get("data", {}).get("polymer_entities", [])
            for entity in polymer_entities:
                entity_id = entity.get(
                    "rcsb_polymer_entity_container_identifiers", {}
                ).get("entity_id", "Not provided")
                sequence = entity.get("entity_poly", {}).get(
                    "pdbx_seq_one_letter_code_can", "Not provided"
                )
                scientific_name = entity.get("rcsb_entity_source_organism", [{}])[
                    0
                ].get("ncbi_scientific_name", "Not provided")

                identifiers = entity.get(
                    "rcsb_polymer_entity_container_identifiers", {}
                ).get("reference_sequence_identifiers", None)
                if identifiers:
                    database_accession = identifiers[0].get(
                        "database_accession", "Not found"
                    )
                    database_name = identifiers[0].get("database_name", "Not found")
                else:
                    database_accession = "Not found"
                    database_name = "Not found"
                if (
                    scientific_name == "Arabidopsis thaliana"
                    and database_accession != "Not found"
                    and database_name != "Not found"
                ):
                    result.append(
                        {
                            "pdb_id": pdb_id,
                            "rcsb_id": rcsb_id,
                            "rcsb_accession": rcsb_accession,
                            "initial_release": initial_release,
                            "title": title,
                            "entity_id": entity_id,
                            "sequence": sequence,
                            "scientific_name": scientific_name,
                            "database_accession": database_accession,
                            "database_name": database_name,
                        }
                    )
        return result

    def save_to_markdown(self, output_file="summary.md"):
        parsed_data = self.parse_json()
        with open(output_file, "w") as md_file:
            for entry in parsed_data:
                md_file.write(f"## PDB ID: {entry['pdb_id']}\n")
                md_file.write(f"RCSB ID: {entry['rcsb_id']}\n")
                md_file.write(f"RCSB Accession: {entry['rcsb_accession']}\n")
                md_file.write(f"Initial Release: {entry['initial_release']}\n")
                md_file.write(f"Title: {entry['title']}\n")
                md_file.write(f"Entity ID: {entry['entity_id']}\n")
                md_file.write(f"Sequence: {entry['sequence']}\n")
                md_file.write(f"Scientific Name: {entry['scientific_name']}\n")
                md_file.write(f"Database Accession: {entry['database_accession']}\n")
                md_file.write(f"Database Name: {entry['database_name']}\n\n")
        print(f"Data has been saved to {output_file}")

    def get_repeated_database_accessions(self):
        data = self.parse_json()
        database_accessions = [entry["database_accession"] for entry in data]
        counter = Counter(database_accessions)
        repeated_accessions = {
            accession: count for accession, count in counter.items() if count > 1
        }
        return repeated_accessions

    def get_unique_database_accessions(self):
        """
        Return a dictionary of unique database accessions and their sequences

        """
        data = self.parse_json()
        unique_accessions = {}
        for entry in data:
            accession = entry["database_accession"]
            if accession not in unique_accessions:
                unique_accessions[accession] = entry["sequence"]
        return unique_accessions


class AF2ProteinDownloader:
    def __init__(self, data_dir="../outputs"):
        self.data_dir = data_dir
        self.pdb_dir = os.path.join(self.data_dir, "pdb")
        self.fasta_dir = os.path.join(self.data_dir, "fasta")

        # Create required directories if they don't exist
        if not os.path.exists(self.pdb_dir):
            os.makedirs(self.pdb_dir)
        if not os.path.exists(self.fasta_dir):
            os.makedirs(self.fasta_dir)

    def download_pdb(self, pdb_id):
        """Download PDB from AF2 and save to self.pdb_dir"""
        pdb = AF2_Structure(pdb_id, self.pdb_dir)
        pdb.download_structure()

    def download_fasta(self, uniprot_id):
        """Download sequence from uniprot and save to self.fasta_dir"""
        protein = A_thaliana_Protein(uniprot_id)
        sequence = protein.sequence
        filename = os.path.join(self.fasta_dir, f"{uniprot_id}.fasta")

        with open(filename, "w") as file:
            file.write(sequence)

        return filename
