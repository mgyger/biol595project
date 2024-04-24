import requests

web_url = "http://www.ebi.ac.uk/pdbe/search/pdb/select?"

def make_request(search_dict, num_rows=5):
    if 'rows' not in search_dict:
        search_dict['rows'] = num_rows
    search_dict['wt'] = 'json'

    response = requests.post(web_url, data=search_dict)

    if response.status_code == 200:
        return response.json()
    else:
        print("[No Data Retrieved - %s] %s")

    return {}


def format_seq_searchterms(sequence, filter_terms):
    params = {
        'json.nl': 'map',
        'start': '0',
        'sort': 'fasta(e_value) asc',
        'sequence': sequence,
        'q': '*:*'
    }
    if filter_terms:
        params['fl'] = ','.join(filter_terms)

    return params


def run_search(sequence, filter_terms=None, num_rows=5):
    search_dict = format_seq_searchterms(sequence=sequence, filter_terms=filter_terms)
    response = make_request(search_dict=search_dict, num_rows=num_rows)
    results = response.get('response', {}).get('docs', [])
    print('Number of results for sequence {}: {}'.format(sequence, len(results)))

    if results:
        pdb_ids = [result.get('pdb_id', '').lower() for result in results]
        return pdb_ids
    else:
        print("No results found for sequence {}.".format(sequence))
        return None


if __name__ == "__main__":
    filter_list = ['pfam_accession', 'pdb_id', 'molecule_name', 'ec_number',
                   'uniprot_accession_best', 'tax_id']

    with open('blast_result.txt', 'r') as file:
        data = file.read()

    # Split the data into individual records
    records = data.split('\n\n')

    fasta_sequences = []
    for record in records:
        lines = record.split('\n')
        fasta_sequence = None
        for line in lines:
            if line.startswith('Fasta Sequence:'):
                fasta_sequence = line.split(': ')[1].strip()
                fasta_sequences.append(fasta_sequence)
                break

    for sequence in fasta_sequences:
        pdbs = run_search(sequence, filter_terms=filter_list)
        if pdbs:
            print("PDB IDs for sequence {}: {}".format(sequence, pdbs))
