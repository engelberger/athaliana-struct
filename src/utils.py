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


class RCSB_Structure:
    def __init__(self, pdb_id, download_dir):
        self.pdb_id = pdb_id
        self.download_dir = download_dir

    def download_structure(self):
        """Download structure from RCSB and save to download_dir"""
        pass
        #


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
        Return a dictionary of unique database accessions, their sequences and corresponding PDB ID

        """
        data = self.parse_json()
        unique_accessions = {}
        for entry in data:
            accession = entry["database_accession"]
            if accession not in unique_accessions:
                unique_accessions[accession] = {
                    "sequence": entry["sequence"],
                    "pdb_id": entry["pdb_id"],
                }
        return unique_accessions


class AF2ProteinDownloader:
    def __init__(self, data_dir="../outputs"):
        self.data_dir = data_dir
        self.pdb_dir_af2 = os.path.join(self.data_dir, "pdb_af2")
        self.pdb_dir_rcsb = os.path.join(self.data_dir, "pdb_rcsb")
        self.fasta_dir = os.path.join(self.data_dir, "fasta")

        # Create required directories if they don't exist
        if not os.path.exists(self.pdb_dir):
            os.makedirs(self.pdb_dir)
        if not os.path.exists(self.fasta_dir):
            os.makedirs(self.fasta_dir)
        if not os.path.exists(self.pdb_dir_rcsb):
            os.makedirs(self.pdb_dir_rcsb)

    def download_rcsb_pdb(self, pdb_id):
        """Download PDB from RCSB and save to self.pdb_dir_rcsb"""
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)
        filename = os.path.join(self.pdb_dir_rcsb, f"{pdb_id}.pdb")

        with open(filename, "w") as file:
            file.write(response.text)

        return filename

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


class Protein:
    """
    Class for a protein
        It
    """

    def __init__(self, uniprot_id: str, pdb_id: str, data_dir: str, sequence_rcsb: str):
        self.uniprot_id = uniprot_id
        self.data_dir = data_dir
        self.protein_dir = os.path.join(data_dir, self.uniprot_id)
        self.sequence_rcsb = sequence_rcsb
        self.pdb_id = pdb_id
        self.af_aligned_sequence = None
        self.rcsb_aligned_sequence = None
        self.date_uploaded_rcsb = None
        self.uniprot_first_deposited_pdb = None
        self.matchmaker_dict = None
        self.chimerax_path = "chimerax"
        os.makedirs(self.protein_dir, exist_ok=True)

    def sequence_uniprot(self):
        filename = os.path.join(self.protein_dir, f"{self.uniprot_id}.fasta")
        if os.path.exists(filename):
            with open(filename) as file:
                return file.read()
        else:
            return self.download_sequence("uniprot")

    def sequence_af2(self):
        filename = os.path.join(self.protein_dir, f"{self.uniprot_id}_AF2.fasta")
        if os.path.exists(filename):
            # If the file is not empty
            if os.path.getsize(filename) > 0:
                with open(filename) as file:
                    return file.read()
        else:
            return self.download_sequence("AF2")

    def download_sequence(self, source: str):
        if source == "uniprot":
            url = f"https://www.uniprot.org/uniprot/{self.uniprot_id}.fasta"
            response = requests.get(url)
            sequence = "".join(response.text.split("\n")[1:])
            filename = os.path.join(self.protein_dir, f"{self.uniprot_id}.fasta")
            with open(filename, "w") as file:
                file.write(f">{self.uniprot_id}\n")
                file.write(sequence)
            return sequence
        elif source == "RCSB":
            # https://www.rcsb.org/fasta/entry/
            url = f"https://www.rcsb.org/fasta/entry/{self.pdb_id}"
        elif source == "AF2":
            url = f"https://alphafold.ebi.ac.uk/files/AF-{self.uniprot_id}-F1-model_v4.fasta"

            response = requests.get(url)
            print(response.text)
            if response.status_code == 200:
                sequence = "".join(response.text.split("\n")[1:])
                filename = os.path.join(
                    self.protein_dir, f"{self.uniprot_id}_{source}.fasta"
                )
                with open(filename, "w") as file:
                    file.write(sequence)
                return sequence

    def predict_structure_with_esm(
        self,
        fasta_file,
        pdb_folder,
        model_dir=None,
        num_recycles=None,
        max_tokens_per_batch=None,
        chunk_size=None,
        cpu_only=False,
        cpu_offload=False,
    ):
        # esm-fold command path
        esm_fold_command = (
            "/home/sc.uni-leipzig.de/nq194gori/micromamba/envs/esmfold/bin/esm-fold"
        )

        # assemble the command
        cmd = [esm_fold_command, "-i", fasta_file, "-o", pdb_folder]

        # Add additional options if provided
        if model_dir is not None:
            cmd += ["-m", model_dir]
        if num_recycles is not None:
            cmd += ["--num-recycles", str(num_recycles)]
        if max_tokens_per_batch is not None:
            cmd += ["--max-tokens-per-batch", str(max_tokens_per_batch)]
        if chunk_size is not None:
            cmd += ["--chunk-size", str(chunk_size)]
        if cpu_only:
            cmd.append("--cpu-only")
        if cpu_offload:
            cmd.append("--cpu-offload")

        # Run the command and capture output
        import subprocess

        print(cmd)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        print(stdout.decode("utf-8"))
        print(stderr.decode("utf-8"))
        if process.returncode != 0:
            raise OSError(f"Command failed: {cmd}, error: {stderr.decode('utf-8')}")

        # create the pdb file path
        pdb_file = os.path.join(pdb_folder, f"{self.uniprot_id}.pdb")
        # Load the pdb file (if exists)
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        import shutil

        # Move and rename the pdb file to the protein directory
        shutil.move(
            pdb_file, os.path.join(self.protein_dir, f"{self.uniprot_id}_ESM.pdb")
        )
        # update the pdb file path
        pdb_file = os.path.join(self.protein_dir, f"{self.uniprot_id}_ESM.pdb")
        # Remove the pdb folder
        shutil.rmtree(pdb_folder)
        # Return the path of the pdb file
        return pdb_file

    @property
    def pdb_structures(self):
        filename = os.path.join(self.protein_dir, f"{self.uniprot_id}_PDB.json")
        if os.path.exists(filename):
            with open(filename) as file:
                return json.load(file)
        else:
            return self.download_pdb_structures()

    def download_pdb_structures(self, pdb_source):
        file_path = os.path.join(
            self.protein_dir, f"{self.uniprot_id}_{pdb_source}.cif"
        )
        if not os.path.exists(file_path):
            if pdb_source == "AF2":
                url = f"https://alphafold.ebi.ac.uk/files/AF-{self.uniprot_id}-F1-model_v4.cif"
            elif pdb_source == "RCSB":
                url = f"https://files.rcsb.org/download/{self.pdb_id}.cif"
            else:
                raise ValueError(f"Unrecognized pdb_source: {pdb_source}")
            response = requests.get(url)
            # Wait 1 second before retrying
            import time

            time.sleep(1)
            # Check if the response is valid
            if response.status_code != 200:
                error_log = "download_error.log"
                with open(error_log, "w") as file:
                    file.write(
                        f"Error downloading PDB structure {self.pdb_id} from {pdb_source}\n"
                    )
                    file.write(response.text)
                raise ValueError(
                    f"Error downloading PDB structure {self.pdb_id} from {pdb_source}"
                )

            with open(file_path, "w") as file:
                file.write(response.text)
            # Check if the file exists or is empty
            if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
                # Write error log
                error_log = "download_error.log"
                with open(error_log, "a") as file:
                    file.write(f"File {self.pdb_id} not created or is empty\n")
                raise ValueError("File not created or is empty")

        return file_path

    def calculate_rmsd_chimeraX(self, pdb_source: str):
        """
        Calculate RMSD between AF2 and RCSB structures using ChimeraX
        """
        import json

        # Download PDB structures
        if pdb_source == "ESM":
            # TODO update if ESM model does not exist it should be predicted
            second_pdb_path = os.path.join(
                self.protein_dir, f"{self.uniprot_id}_{pdb_source}.pdb"
            )
        if pdb_source == "AF2":
            second_pdb_path = self.download_pdb_structures("AF2")
        rcsb_pdb_path = self.download_pdb_structures("RCSB")

        # If the pdb files exist, calculate RMSD
        if os.path.exists(second_pdb_path) and os.path.exists(rcsb_pdb_path):
            # Calculate RMSD with matchmaker
            # chimerax --offscreen --script "matchmaker.py --dir /home/iwe30/Remote/nq194gori/github/athaliana-struct/notebooks/proteins_data/P93311 --file_1 P93311_AF2.cif --file_2 P93311_RCSB.cif --output matchmaker_full_P93311.html"
            cmd = f"{self.chimerax_path} --script '../src/matchmaker.py --dir {self.protein_dir} --file_1 {second_pdb_path} --file_2 {rcsb_pdb_path} --output matchmaker_full_{self.uniprot_id}.html'"
            # Save the command to a file
            with open(
                os.path.join(self.protein_dir, f"matchmaker_full_{self.uniprot_id}.sh"),
                "w",
            ) as file:
                file.write(cmd)

            # If there is already a matchmaker_full_Q95747.json file, don't calculate RMSD
            if not os.path.exists(
                os.path.join(
                    self.protein_dir, f"matchmaker_full_{self.uniprot_id}.json"
                )
            ):
                # Run command with subprocess
                import subprocess

                process = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                stdout, stderr = process.communicate()
                # If there is an error, print it
                if stderr:
                    print(stderr)
                html_file = os.path.join(
                    self.protein_dir, f"matchmaker_full_{self.uniprot_id}.html"
                )
                from bs4 import BeautifulSoup, NavigableString
                import json

                # Load html file
                with open(html_file) as f:
                    html_doc = f.read()

                # Parse html file
                soup = BeautifulSoup(html_doc, "html.parser")

                data = {}

                brs = soup.find_all("br")

                lines = [
                    br.next_sibling
                    for br in brs
                    if isinstance(br.next_sibling, NavigableString)
                ]

                lines = [str(line).strip() for line in lines]

                for i, line in enumerate(lines):
                    if line.startswith("Matchmaker"):
                        start = i
                        break

                end = lines.index("Residues:")

                relevant_lines = lines[start:end]

                for line in relevant_lines:
                    if "Matchmaker" in line:
                        chains = line.split(" with ")
                        data["Chain 1"] = (
                            chains[0].split(" ")[1].split(",")[0]
                        )  # keep only the relevant chain ID
                        data["Chain 2"] = (
                            chains[1].split(",")[0].strip()
                        )  # keep only the relevant chain ID
                        score = line.split("=")[-1].strip()
                        data["Sequence alignment score"] = float(score)

                    elif "Alignment identifier" in line:
                        data["Alignment identifier"] = int(line.split(" ")[-1].strip())

                    elif "Chains used in RMSD evaluation" in line:
                        # Chains used in RMSD evaluation for alignment 1: P93311_AF2.cif #1/A, P93311_RCSB.cif #2/Aa
                        # Save only the chain IDs ex. A and Aa
                        chains = line.split(":")[-1].strip().split(",")
                        data["Chain 1 for RMSD"] = chains[0].split("/")[-1]
                        data["Chain 2 for RMSD"] = chains[1].split("/")[-1]

                    elif "RMSD between" in line:
                        rmsd = line.split("is")[1].split(";")[0].strip()
                        # RMSD between 61 pruned atom pairs is 1.072 angstroms; (across all 96 pairs: 17.267)
                        data["RMSD pruned"] = float(rmsd.split(" ")[0].strip())
                        # Save the number of CA atoms used for RMSD pruned
                        data["Number of CA atoms used for RMSD pruned"] = str(
                            line.split(" ")[2].strip()
                        )
                        all_pairs = line.split(" ")[-1].replace(")", "").strip()
                        data["RMSD all pairs"] = float(all_pairs)
                        data["Number of CA atoms used for RMSD all"] = str(
                            line.split(" ")[-3].strip()
                        )

                    elif "Sequences:" in line:
                        seq_num = relevant_lines.index(line)
                        seq1 = relevant_lines[seq_num + 1].split(" ", 1)
                        seq2 = relevant_lines[seq_num + 2].split(" ", 1)
                        # Removing string "Chain A" by splitting and join remaining sequence back
                        # Same thing applied for sequence 2
                        seq1_clean = " ".join(seq1[1].split())
                        seq2_clean = " ".join(seq2[1].split())
                        data["Sequence 1"] = seq1_clean.split(" ")[-1]
                        data["Sequence 2"] = seq2_clean.split(" ")[-1]

                self.matchmaker_dict = data
                # Save it to json file
                with open(
                    os.path.join(
                        self.protein_dir, f"matchmaker_full_{self.uniprot_id}.json"
                    ),
                    "w",
                ) as fp:
                    json.dump(data, fp)
                return data
            else:
                # Load json file if it exists
                with open(
                    os.path.join(
                        self.protein_dir, f"matchmaker_full_{self.uniprot_id}.json"
                    ),
                    "r",
                ) as fp:
                    data = json.load(fp)
                    return data
