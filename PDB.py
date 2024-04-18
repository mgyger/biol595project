"""---------------------------------------------------------------------------------------------------------------------
PDB database pull code

BIOL 595 Final Project
Morgan Gyger 15 04 2024
---------------------------------------------------------------------------------------------------------------------"""

import requests
import py3Dmol

def send_query(query_structure):
    api_endpoint = "https://search.rcsb.org/rcsbsearch/v2/query"
    response = requests.post(api_endpoint, json=query_structure)
    if response.status_code == 200:
        return response.json()
    else:
        print("Error:", response.text)
        return None

def create_pdb_model(pdb_id):
    view = py3Dmol.view(width=800, height=400)
    view.addModel(pdb_file, "pdb")
    view.setStyle({"cartoon": {'color':'spectrum'}})
    view.zoomTo()
    view.show()

def main(fasta_string):
    # Parse FASTA sequences
    sequences = fasta_string.strip().split('>')[1:]

    # Iterate over the first two FASTA sequences
    for i, sequence in enumerate(sequences[:2], start=1):
        # Construct query payload based on the provided query structure
        query_structure = {
  "query": {
    "type": "terminal",
    "service": "text",
    "parameters": {
      "attribute": "rcsb_uniprot_protein.name.value",
      "operator": "exact_match",
      "value": "Free fatty acid receptor 2"
    }
  },
  "return_type": "entry",
  "request_options": {
    "results_content_type": [
      "computational",
      "experimental"
    ]
  }
}

        # Send query to PDB API
        print(f"Sending query for FASTA sequence {i}:")
        response_data = send_query(query_structure)
        if response_data:
            search_data = response_data.get('data', {}).get('search', {})
            pdb_ids = [entry['identifier'] for entry in search_data.get('edges', [])]
            print(f"Matching PDB IDs for FASTA sequence {i}: {pdb_ids}")
            if pdb_ids:
                # Create PDB model for the first matching PDB ID
                print("Creating model for the first matching PDB ID:")
                create_pdb_model(pdb_ids[0])
            else:
                print("No matching PDB IDs found.")
        else:
            print("No response data received.")

# Example FASTA string
fasta_string = """MSFRFGQHLIKPSVVFLKTELSFALVNRKPVVPGHVLVCPLRPVERFHDLRPDEVADLFQTTQRVGTVVEKHFHGTSLTFSMQDGPEAGQTVKHVHVHVLPRKAGDFHRNDSIYEELQKHDKEDFPASWRSEEEMAAEAAALRVYFQ"""

main(fasta_string)

