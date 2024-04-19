import requests
import py3Dmol

# Fetch PDB file content by PDB ID
def fetch_pdb(pdb_id):
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(pdb_url)
    if response.status_code == 200:
        return response.text
    else:
        print("Error fetching PDB file:", response.status_code)
        return None

# Create model of PDB using py3Dmol
def create_pdb_model(pdb_id):
    pdb_content = fetch_pdb(pdb_id)
    if pdb_content:
        view = py3Dmol.view(width=800, height=400)
        view.addModel(pdb_content, "pdb")
        view.setStyle({"cartoon": {'color': 'spectrum'}})
        view.zoomTo()
        view.show()
    else:
        print("Failed to fetch PDB file for", pdb_id)

# Read blast results from text file and extract PDB IDs
def read_blast_results(file_path):
    with open(file_path, 'r') as file:
        blast_results = file.readlines()
    pdb_ids = [line.strip().split('\t')[1] for line in blast_results]  # Assuming the second column contains PDB IDs
    return pdb_ids

# Main function
def main(blast_results_file, fasta_sequence):
    # Display PDB models for BLAST hits
    blast_hits = read_blast_results(blast_results_file)
    print("PDB IDs from BLAST results:", blast_hits)
    print("Displaying PDB models:")
    for pdb_id in blast_hits:
        create_pdb_model(pdb_id)

    # Use AlphaFold to predict structure for FASTA sequence
    # Replace this with your actual AlphaFold prediction implementation
    alpha_fold_model = "https://alphafold.ebi.ac.uk/files/AF-Q8VCK6-F1-model_v4.cif"
    print("Displaying AlphaFold predicted structure:")
    view = py3Dmol.view(width=800, height=400)
    view.addModel(alpha_fold_model, "cif")
    view.setStyle({"cartoon": {'color': 'spectrum'}})
    view.zoomTo()
    view.show()

#inputs
blast_results_file = "blast_results.txt"

# Run the main function
main(blast_results_file, fasta_sequence)
