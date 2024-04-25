import requests
import json
import dash_bio as dashbio
from dash import Dash, html, Input, Output, callback
from dash_bio.utils import create_mol3d_style
from io import StringIO


#### Get top PDB ID matches for FASTA outputs ####
# get fasta from blast results
fasta_sequences = []
with open('blast_result.txt', 'r') as file:
    data = file.read()

# split the data
records = data.split('\n\n')

for record in records:
    lines = record.split('\n')
    for line in lines:
        if line.startswith('Fasta Sequence:'):
            fasta_sequence = line.split(': ')[1].strip()
            fasta_sequences.append(fasta_sequence)
            break

# set up to only get the first, closest result from pdb
top_identifiers = []

# set up search query and define
for fasta_sequence in fasta_sequences:
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1,
                "identity_cutoff": 0.9,
                "sequence_type": "protein",
                "value": fasta_sequence
            }
        },
        "request_options": {
            "scoring_strategy": "sequence"
        },
        "return_type": "entry"
    }

    # encode the query as JSON
    encoded_query = json.dumps(query)

    # construct the api url for searches of interest
    api_url = "https://search.rcsb.org/rcsbsearch/v2/query?json=" + encoded_query

    # set GET request to the API
    response = requests.get(api_url)

    # parse JSON response
    json_response = response.json()

    # extract the top identifier if match exists
    top_identifier = None
    if "result_set" in json_response and len(json_response["result_set"]) > 0:
        top_identifier = json_response["result_set"][0]["identifier"]

    # append top identifier to the list of all matches for all sequences
    top_identifiers.append(top_identifier)

#### Set up Dash ####
app = Dash(__name__)


# fetch PDB content
def fetch_pdb_content(pdb_id):
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for 4xx or 5xx status codes
        return response.text
    except requests.RequestException as e:
        print(f"Failed to fetch PDB content for {pdb_id}. Error: {e}")
        return None


# parse PDB content
def parse_pdb_content(pdb_content):
    # create file-like object from the PDB content string
    pdb_file = StringIO(pdb_content)

    # parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure("pdb_data", pdb_file)

    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # convert float32 coordinates to regular floats
                    coord = [float(coord) for coord in atom.coord]
                    atoms.append({
                        "elem": atom.element,
                        "chain": chain.id,
                        "residue_name": residue.resname,
                        "residue_number": residue.id[1],
                        "x": coord[0],
                        "y": coord[1],
                        "z": coord[2]
                    })

    return atoms


# fetch PDB content for the first result (initial query)
first_pdb_id = top_identifiers[0] if top_identifiers else None
pdb_content = fetch_pdb_content(first_pdb_id)

# parse PDB content
if pdb_content:
    atoms = parse_pdb_content(pdb_content)

else:
    atoms = []


#### timer set up ######
# callback to update the timer every second
@app.callback(
    Output('timer-output', 'children'),
    [Input('interval-component', 'n_intervals')]
)
def update_timer(n):
    return f"Page loaded {n} seconds ago."


# define the interval component to update the timer every 5 seconds
app.layout = html.Div([
    dcc.Interval(id='interval-component', interval=1000, n_intervals=0),
    html.Button('Update', id='update-button'),  # Button to trigger update
    html.Div(id='timer-output')  # container for displaying the timer
])

#### Molecule3DViewer ####
styles = create_mol3d_style()

app.layout = html.Div([
    dashbio.Molecule3dViewer(
        id='dashbio-default-molecule3d',
        modelData={},
        styles=styles
    ),
    "Selection data",
    html.Hr(),
    html.Div(id='default-molecule3d-output')
])

@app.callback(
    Output('default-molecule3d-output', 'children'),
    Input('dashbio-default-molecule3d', 'selectedAtomIds')
)

def show_selected_atoms(atom_ids):
    if atom_ids is None or len(atom_ids) == 0:
        return 'No atom has been selected. Click somewhere on the molecular \
        structure to select an atom.'
    return [html.Div([
        html.Div('Element: {}'.format(atoms['atoms'][atm]['elem'])),
        html.Div('Chain: {}'.format(atoms['atoms'][atm]['chain'])),
        html.Div('Residue name: {}'.format(atoms['atoms'][atm]['residue_name'])),
        html.Br()
    ]) for atm in atom_ids]


#### Main Code ####
if __name__ == '__main__':
    # print top identifiers
    for i, top_identifier in enumerate(top_identifiers):
        print(f"Top Identifier for sequence {i + 1}: {top_identifier if top_identifier else 'None'}")
        app.run_server(debug=True)
