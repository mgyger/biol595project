"""---------------------------------------------------------------------------------------------------------------------
PDB database pull code

BIOL 595 Final Project
Morgan Gyger 15 04 2024
---------------------------------------------------------------------------------------------------------------------"""

import requests
import json
import py3Dmol


def predict_protein_structure(sequence):
    api_url = "https://www.predictprotein.org/api/v2/"
    endpoint = "predict/proteinstructure"

    # Request payload
    data = {
        "input_data": {
            "sequence": sequence,
            "format": "single_sequence",
            "submit": "Predict"
        }
    }

    # Send POST request
    response = requests.post(f"{api_url}{endpoint}", json=data)

    if response.status_code == 200:
        result = response.json()
        if result.get("result"):
            return result["result"]
        else:
            print("Error in prediction:", result.get("message", "Unknown error"))
            return None
    else:
        print("Error:", response.status_code)
        return None


# Read FASTA file and extract the first sequence
with open("blast_result.txt", "r") as file:
    fasta_data = file.read().strip().split("\n")

# Extract the first sequence
first_sequence = fasta_data[1]

# Predict protein structure
predicted_structure = predict_protein_structure(first_sequence)

if predicted_structure:
    print("Predicted protein structure:")
    print(json.dumps(predicted_structure, indent=2))


def fetch_pdb_ids(sequence):
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 0.001,
                "identity_cutoff": 95,
                "target": "pdb_protein_sequence",
                "value": sequence
            }
        },
        "request_options": {
            "pager": {
                "start": 0,
                "rows": 10
            },
            "scoring_strategy": "combined"
        },
        "return_type": "entry"
    }

    response = requests.post(search_url, json=query)
    if response.status_code == 200:
        data = response.json()
        pdb_ids = [entry["identifier"] for entry in data["result_set"]]
        return pdb_ids
    else:
        print("Error fetching PDB IDs:", response.status_code)
        return []


# Read FASTA file
with open("blast_result.txt", "r") as file:
    fasta_data = file.read().strip().split("\n")

# Extract sequences from FASTA file including the first one
sequences = ["".join(fasta_data[i]) for i in range(1, len(fasta_data), 2)]

# Fetch PDB IDs for each sequence and print
for sequence in sequences:
    pdb_ids = fetch_pdb_ids(sequence)
    print("PDB IDs for sequence:")
    print(pdb_ids)

    # Visualize PDB structures using py3Dmol
    if pdb_ids:
        viewer = py3Dmol.view(width=800, height=600)
        for pdb_id in pdb_ids:
            viewer.addModel(f"https://files.rcsb.org/download/{pdb_id}.pdb", "pdb")
        viewer.setStyle({"cartoon": {"color": "spectrum"}})
        viewer.zoomTo()
        viewer.show()
