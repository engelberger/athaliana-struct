# A. thaliana Protein Structure Prediction Using AlphaFold

## Objective
The main aim is to use AlphaFold predictions to analyze a dataset of proteins from A. thaliana, aiming at a large-scale prediction and scoring process of non-synonymous protein variants.

## Description
This operation will involve the utilization of AlphaFold for the prediction of single-chain protein structures, rather than complex structures, for practical reasons. This approach simplifies calculations and eases the setup.

The process focuses on protein structures with a certain degree of completeness (e.g., 30% or 50% sequence consistency as per UniProt), avoiding the potential waste of resources on smaller fragments from larger proteins. 

Moreover, we propose sequence similarity cut-off (e.g., 99%) to reduce redundancies and elude multiple calculations for the same protein chain.

## Plan and Approach
1. Acquire the proper proteins: AlphaFold has a pre-computed structure database that should be used for the operation. Based on the information from PDB, we've found 830 proteins related to A. thaliana.

   1. Filter proteins that were made available on or after April 30, 2018, to ensure that they are not part of the AlphaFold's training set.

   2. Obtain single chain proteins over complexes to simplify the study.

   3. Choose proteins where the PDBs is a representation of a certain threshold of the entire protein, avoiding fragments.

   4. Eliminate any redundant protein chains with a sequence similarity cut-off.

2. Implement and run AlphaFold: Utilizing the acquired, filtered proteins as input, AlphaFold will then predict the protein structures.

3. Scoring and Analysis: Compare predicted structures with actual structures where they are available. Utilize the standard methods for comparing predicted and empirical structures.

The ultimate goal is to demonstrate the effectiveness of AlphaFold for protein-structure prediction specifically for A. thaliana proteins.  The output can be used for further research in the functional analysis for A. thaliana, or as a part of bioinformatic research on plant proteins. Overall, it will allow us to understand AlphaFold's performance when given a large-scale task on A. thaliana proteins, establishing a foundation for future research on similar topics.

## Milestones
1. Selection and download of relevant proteins from the pre computed AF2 database

2. Implementation of AlphaFold and generation of protein structure predictions

3. Comparison of predicted and actual structures to assess the success and reliability of AlphaFold predictions.

4. Comprehensive analysis and documentation of findings

With the plan and milestones clearly laid out, the prediction of the large-scale A. thaliana protein variants can be accomplished efficiently and in an organized manner.


# Protein Processing and Downloading Documentation

This Python script is designed for processing and downloading protein sequences and structures for the A. thaliana plant.

## Getting Started

### Prerequisites

To run this script, you'll need:

- Python3
- Internet connection
- Access to uniprot and AlphaFold databases.
- Python's `requests` and `os` libraries

# Implementation

Class Documentation

## 1. A_thaliana_Protein

### Usage:

```python
protein = A_thaliana_Protein(uniprot_id)
```

The `A_thaliana_Protein` class represents a protein object from the specie *A. thaliana*.

- **Attributes**:

    - `uniprot_id (str)`: The ID of the protein sequence from UniProt Database.
    - `sequence (str)`: The protein sequence downloaded from the UniProt API.

- **Methods**:

    - `_download_sequence() -> str`:
        Privately called to download sequence from UniProt API and return as a string.

- **Example**:

```python
protein = A_thaliana_Protein('Q8WZ42')
print(protein.sequence)
```

---

## 2. AF2_Structure

### Usage:

```python
structure = AF2_Structure(pdb_id, download_dir)
```

The `AF2_Structure` class represents a structure object in AlphaFold 2 (AF2) Database.

- **Attributes**:

    - `pdb_id (str)`: The ID of the protein structure from Protein Data Bank (PDB).
    - `download_dir (str)`: The directory where to save the downloaded structure.

- **Methods**:

    - `_download_structure() -> str`:
        Privately called to download protein structure from AF2 API and save to `download_dir`. Returns the file path of the downloaded structure.

- **Example**:

```python
structure = AF2_Structure('6LU7', '/path/to/your/directory')
file_path = structure.download_structure()
print(file_path)
```

---

## 3. ProteinCompleteness

### Usage:

```python
content = Protein_Completeness(pdb_protein_sequence, uniprot_sequence)
```

The `ProteinCompleteness` class represents an operational object for determining protein completeness.

- **Attributes**:

    - `pdb_protein_sequence (str)`: The sequence string of the protein from PDB.
    - `uniprot_sequence (str)`: The sequence string of the protein from UniProt.

- **Methods**:

    - `completeness_check(completeness_threshold: float) -> bool`:
        Calculates the completeness percentage of `pdb_protein_sequence` with respect to `uniprot_sequence` and returns whether it satisfies the `completeness_threshold`.

- **Example**:

```python
content = ProteinCompleteness('ACTG...', 'ACGTG...')
is_complete = content.completeness_check(50)
print(is_complete)
```

---

## 4. SequenceSimilarity

### Usage:

```python
similarity = SequenceSimilarity(protein_sequence1, protein_sequence2)
```

The `SequenceSimilarity` class is an operational object to determine the similarity between two protein sequences.

- **Attributes**:

    - `protein_sequence1 (str)`: The first protein sequence string for comparison.
    - `protein_sequence2 (str)`: The second protein sequence string for comparison.

- **Methods**:

    - `calculate_similarity(similarity_threshold: float) -> bool`:
        Calculates the sequence similarity between `protein_sequence1` and `protein_sequence2` and returns whether it satisfies the `similarity_threshold`.

- **Example**:

```python
similarity = SequenceSimilarity('ACTG...', 'ACGTG...')
is_similar = similarity.calculate_similarity(99)
print(is_similar)
```

## 5. PDBJsonParser

### Usage:

```python
parser = PDBJsonParser(json_file)
```

The `PDBJsonParser` class used to parse the JSON data from a local JSON file which contains PDB data.

- **Attributes**:

    - `json_file (str)`: The path to the JSON file.
    - `data (dict)`: Parsed JSON data. Initialized when the object is created.

- **Methods**:

    - `load_json() -> dict`: A private method called at initialization to parse and load JSON data from the provided `json_file`.
    
    - `extract_info() -> list`: Extracts key information from the JSON data and returns a list of dictionary objects. Each dictionary contains information of a protein.

- **Example**:

```python
parser = PDBJsonParser('path_to_your_json_file.json')
parsed_data = parser.extract_info()
print(parsed_data)
```

---

These classes are designed to facilitate the process of downloading and analyzing the protein data. Their functionality simplifies the Protein Data Bank and UniProt data extraction, protein completeness checking, and protein sequence similarity calculation.

The classes can be used in harmonious combination to download protein sequences, analyze completeness, and compare sequence similarities, saving considerable time and increasing code readability and accessibility.




