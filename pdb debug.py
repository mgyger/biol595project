import time
def fetch_pdb_ids(sequence, max_retries=3, retry_delay=1):
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

    for attempt in range(max_retries):
        response = requests.post(search_url, json=query)
        if response.status_code == 200:
            data = response.json()
            pdb_ids = [entry["identifier"] for entry in data["result_set"]]
            return pdb_ids
        elif response.status_code == 503:
            print(f"Retrying ({attempt + 1}/{max_retries}) after {retry_delay} seconds...")
            time.sleep(retry_delay)
        else:
            print("Error fetching PDB IDs:", response.status_code)
            return []

    print("Max retries exceeded. Unable to fetch PDB IDs.")
    return []
