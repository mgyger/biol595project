import requests
import re
from bs4 import BeautifulSoup
from Bio.Blast import NCBIWWW, NCBIXML
import sqlite3
import json
from urllib.error import HTTPError
from urllib.request import urlopen
from Bio.PDB import PDBParser
from io import StringIO

class SOPMA: # done by ANDY KECK
    """Class for running SOPMA analysis on protein sequences and storing results."""
    def write_to_database(self, accession_number, secondary_structure):
        """Write SOPMA analysis results to a SQLite database.

        Args:
            accession_number (str): The accession number of the protein sequence.
            secondary_structure (list): List of secondary structure information.

        Returns:
            None
        """
        conn = sqlite3.connect('outputs_data.db')
        cursor = conn.cursor()

        secondary_structure_text = '\n'.join(secondary_structure)

        cursor.execute("INSERT INTO sopma_results VALUES (?, ?)", (accession_number, secondary_structure_text))
        conn.commit()
        conn.close()

    def __init__(self):
        """Initialize SOPMA object with base URL for SOPMA analysis."""
        self.base_url = "https://npsa.lyon.inserm.fr/cgi-bin/npsa_automat.pl?page=/NPSA/npsa_sopma.html"

    def run_sopma(self, sequence, accession_number, callback):
        """Run SOPMA analysis on a protein sequence.

        Args:
            sequence (str): Protein sequence in FASTA format.
            accession_number (str): Accession number of the protein sequence.
            callback (function): Callback function to handle SOPMA results.

        Returns:
            None
        """
        # Make get request
        response = requests.get(self.base_url)

        # Check if request was successful
        if response.status_code == 200:
            print("SOPMA page retrieved successfully.")
            soup = BeautifulSoup(response.text, 'html.parser')

            form = soup.find('form')

            form_action = form.get('action')

            # Prepare data
            data = {
                'notice': sequence,
                'states': '3',
                'width': '17',
                'threshold': '8',
                'title': accession_number,
                'ali_width': '70'
            }

            # Make a POST request to the form action URL
            response = requests.post("https://npsa.lyon.inserm.fr" + form_action, data=data)

            # Check if request was successful
            if response.status_code == 200:
                print("SOPMA analysis successful.")
                # Extract the secondary structure information
                secondary_structure_text = self.extract_secondary_structure(response.text)
                callback(accession_number, secondary_structure_text)  # Pass the correct arguments here
            else:
                print(f"Error: Failed to retrieve data from SOPMA for {accession_number}")
        else:
            print("Error: Failed to retrieve SOPMA page")

    def extract_secondary_structure(self, html_content):
        """Extract secondary structure information from HTML content.

        Args:
            html_content (str): HTML content of the SOPMA analysis result page.

        Returns:
            list: List of secondary structure information.
        """
        structure_names = {
            'h': 'Alpha helix',
            'g': '310 helix',
            'i': 'Pi helix',
            'b': 'Beta bridge',
            'e': 'Extended strand',
            't': 'Beta turn',
            's': 'Bend region',
            'c': 'Random coil'
        }

        # HTML parsing
        soup = BeautifulSoup(html_content, 'html.parser')

        # Find pre tag
        pre_tag = soup.find('pre')

        secondary_structure_text = pre_tag.get_text()

        # Count lowercase letters
        counts = {}
        total_chars = 0
        for char in secondary_structure_text:
            if char.islower():
                counts[char] = counts.get(char, 0) + 1
                total_chars += 1

        # Calculate and return ratio
        results = []
        for char, count in counts.items():
            ratio = (count / total_chars) * 100
            structure_name = structure_names[char]
            results.append(f"{structure_name}: {count} is {ratio:.2f}%\n")

        return results

    def write_to_file(self, data, accession_number):
        """Write SOPMA analysis results to a text file.

        Args:
            data (list): List of secondary structure information.
            accession_number (str): Accession number of the protein sequence.

        Returns:
            None
        """
        with open('sopma_result.txt', 'a') as file:
            file.write(f"Accession Number: {accession_number}\n")
            for result in data:
                file.write(result)
            file.write('\n')

    def fasta_to_amino_acid(self, fasta_sequence):
        """Translate a protein sequence from FASTA format to amino acid sequence.

        Args:
            fasta_sequence (str): Protein sequence in FASTA format.

        Returns:
            str: Amino acid sequence.
        """

        codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }

        # Translate into amino acids
        amino_acids = ''
        for i in range(0, len(fasta_sequence), 3):
            codon = fasta_sequence[i:i+3]  # Get codon
            if codon in codon_table:
                amino_acid = codon_table[codon]
                amino_acids += amino_acid

        return amino_acids

    def read_sequence_file(self, filename):
        """Read protein sequences from a file.

        Args:
            filename (str): Name of the file containing protein sequences.

        Returns:
            list: List of tuples containing accession numbers and protein sequences.
        """
        sequences = []
        with open(filename, 'r') as file:
            accession_number = ''
            fasta_sequence = ''
            for line in file:
                line = line.strip()  # Remove whitespace
                if line.startswith("Accession Number:"):
                    accession_number = line.split(": ")[1]
                elif line.startswith("Fasta Sequence:"):
                    fasta_sequence = line.split(": ")[1]
                    fasta_sequence = fasta_sequence.replace("\n", "")  # Remove newline characters
                    sequences.append((accession_number, fasta_sequence))
        return sequences

class NCBI_BLAST: # done by both ANDY KECK and DAYE KWON
    """Class for running NCBI BLAST search and storing results."""
    def write_to_database(self, top_results):
        """Write NCBI BLAST search results to a SQLite database.

        Args:
            top_results (list): List of top BLAST search results.

        Returns:
            None
        """
        conn = sqlite3.connect('outputs_data.db')
        cursor = conn.cursor()
        for result in top_results:
            gene_name, accession, scientific_name, e_score, alignment_score, identity, fasta_sequence = result
            cursor.execute("INSERT INTO blast_results VALUES (?, ?, ?, ?, ?, ?, ?)",
                           (gene_name, accession, scientific_name, e_score, alignment_score, identity, fasta_sequence))
        conn.commit()
        conn.close()
    def run_ncbi_blast(self, sequence):
        """Run NCBI BLAST search.

        Args:
            sequence (str): Genomic sequence to search.

        Returns:
            list: List of top BLAST search results.
        """
        # Perform BLAST search on NCBI
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        blast_records = NCBIXML.parse(result_handle)

        top_matches = []
        for record in blast_records:
            # Iterate through each blast record
            for alignment in record.alignments[:5]:  # Limit to 5 results
                # Extract information from alignment
                gene_name = alignment.title.split("|")[4].split()[0]  # Extract gene name
                accession = alignment.accession  # Extract accession number
                scientific_name = alignment.hit_def.split("[")[0].strip()  # Extract scientific name
                e_score = alignment.hsps[0].expect  # Extract E score
                alignment_score = alignment.hsps[0].score  # Extract alignment score
                identity = alignment.hsps[0].identities / alignment.hsps[0].align_length * 100  # Calculate identity (%)
                fasta_sequence = alignment.hsps[0].sbjct  # Extract fasta sequence

                top_matches.append(
                    (gene_name, accession, scientific_name, e_score, alignment_score, identity, fasta_sequence))

        return top_matches

    def get_genomic_sequence_from_user(self):
        """Get genomic sequence input from the user.

        Returns:
            str: Genomic sequence input by the user.
        """
        sequence = input("Enter the genomic sequence: ")
        if all(nucleotide in ['A', 'G', 'C', 'T'] for nucleotide in sequence.upper()):
            return sequence
        else:
            print("Please input a sequence containing only A, G, C, or T nucleotides.")
        return sequence

    def main(self, sequence):
        """Main function to run NCBI BLAST search and store results.

        Args:
            sequence (str): Genomic sequence to search.

        Returns:
            None
        """
        # Run NCBI BLAST search and retrieve top 5 results
        top_results = self.run_ncbi_blast(sequence)

        # Write data to file
        with open("blast_result.txt", "w") as file:
            # Write user input sequence at the top
            file.write("User Input Sequence:\n")
            file.write(sequence + "\n\n")

            # Write blast results
            for i, top_result in enumerate(top_results):
                file.write(f"Result {i + 1}:\n")
                gene_name, accession, scientific_name, e_score, alignment_score, identity, fasta_sequence = top_result
                file.write("Gene Name: " + gene_name + "\n")
                file.write("Accession Number: " + accession + "\n")
                file.write("Scientific Name: " + scientific_name + "\n")
                file.write("E Score: " + str(e_score) + "\n")
                file.write("Alignment Score: " + str(alignment_score) + "\n")
                file.write("Identity (%): " + str(identity) + "\n")
                file.write("Fasta Sequence: " + fasta_sequence + "\n")
                file.write("\n")

        self.write_to_database(top_results)

        print("Data written to blast_result.txt.")
class UniprotFetcher: # done by both ANDY KECK and DAYE KWON
    """Class for fetching UniProt entries and storing results."""
    def __init__(self, blast_file_path):
        """Initialize UniprotFetcher object with blast file path."""
        self.blast_file_path = blast_file_path
        self.accession_numbers = []

    def _extract_accession_numbers(self):
        """Extract UniProt accession numbers from BLAST result file."""
        with open(self.blast_file_path, "r") as blast_file:
            for line in blast_file:
                accession_number_match = re.match(r'Accession Number: (.+)', line)
                if accession_number_match:
                    accession_number = accession_number_match.group(1)
                    if accession_number != "Entry":
                        self.accession_numbers.append(accession_number)

    def fetch_uniprot_entries(self, nucleotide_accession):
        """Fetch UniProt entries using nucleotide accession number.

        Args:
            nucleotide_accession (str): Nucleotide accession number.

        Returns:
            list: List of UniProt entries.
        """
        url = f"https://rest.uniprot.org/uniprotkb/search?query={nucleotide_accession}&format=tsv"
        response = requests.get(url)

        if response.status_code == 200:
            lines = response.text.split('\n')
            entries = [line.split('\t') for line in lines]
            return entries
        else:
            print(f"Error: {response.status_code} - {response.reason}")
            return None

    def fetch_and_write_uniprot_entries(self, output_file_path="uniprot_entry_codes.txt"):
        """Fetch and write UniProt entries to a file.

        Args:
            output_file_path (str, optional): Path to write UniProt entries file. Defaults to "uniprot_entry_codes.txt".

        Returns:
            None
        """
        self._extract_accession_numbers()

        all_uniprot_entries = set()  # Use a set to store unique entry codes
        for accession_number in self.accession_numbers:
            uniprot_entries = self.fetch_uniprot_entries(accession_number)
            if uniprot_entries:
                for entry in uniprot_entries:
                    if entry and entry[0].strip() and entry[0].strip() != "Entry":  # Exclude 'Entry'
                        all_uniprot_entries.add(entry[0])

        with open(output_file_path, "w") as f:
            for entry_code in all_uniprot_entries:
                f.write("UniProt Entry Code: " + entry_code + "\n")
        print("UniProt entry codes written to", output_file_path)

class DeepGOPredictor: # done by DAYE KWON
    """Class for predicting gene ontology (GO) functions using DeepGO."""
    def __init__(self):
        """Initialize DeepGOPredictor object with API URL and headers."""
        self.url = 'https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create'
        self.headers = {'Content-Type': 'application/json'}

    def predict_go_functions(self, sequence, accession_number, threshold=0.4):
        """Predict GO functions using DeepGO.

        Args:
            sequence (str): Protein sequence in FASTA format.
            accession_number (str): Accession number of the protein sequence.
            threshold (float, optional): Threshold for filtering predictions. Defaults to 0.4.

        Returns:
            list: List of filtered GO function predictions.
        """
        data = {
            "version": "1.0.18",
            "data_format": "fasta",  # must be small letters
            "data": sequence,
            "threshold": threshold
        }

        # Send post request
        response = requests.post(self.url, json=data, headers=self.headers)

        # Check response
        if response.status_code == 200 or response.status_code == 201:
            result = response.json()
            predictions = result['predictions'][0]['functions']

            # Filter predictions with score >= 0.5
            filtered_predictions = [(category['name'], [(go_id, description, round(score, 3)) for go_id, description, score in category['functions'] if score >= 0.5]) for category in predictions]

            # Write to file if there are filtered predictions
            if any(filtered_predictions):
                with open('deepgo_result.txt', 'a') as file:
                    file.write(f"Accession Number: {accession_number}\n")
                    for category_name, functions in filtered_predictions:
                        if functions:
                            file.write(category_name + '\n')
                            for func in functions:
                                go_id, description, score = func
                                file.write(f"{go_id} - {description} - {score}\n")
                            file.write('\n')
            return filtered_predictions
        else:
            print(f'Failed to retrieve data for {accession_number}:', response.status_code)
            print(response.text)
            return None

    def main(self, sequences):
        """Main function to predict GO functions for protein sequences and store results.

        Args:
            sequences (list): List of tuples containing accession numbers and protein sequences.

        Returns:
            None
        """
        for accession_number, fasta_sequence in sequences:
            amino_acid_sequence = SOPMA().fasta_to_amino_acid(fasta_sequence)
            print("Accession Number:", accession_number)
            print("Amino Acid Sequence:", amino_acid_sequence)
            predictions = self.predict_go_functions(amino_acid_sequence, accession_number)
            if predictions:
                self.write_to_database(accession_number, predictions)
            else:
                print("No predictions found for accession number:", accession_number)

        print("Data written to deepgo_results table in database.")

    def write_to_database(self, accession_number, predictions):
        """Write DeepGO predictions to a SQLite database.

        Args:
            accession_number (str): Accession number of the protein sequence.
            predictions (list): List of GO function predictions.

        Returns:
            None
        """
        if predictions is not None:  # Add a check for None
            conn = sqlite3.connect('outputs_data.db')
            cursor = conn.cursor()
            for category_name, functions in predictions:
                for go_id, description, score in functions:
                    cursor.execute("INSERT INTO deepgo_results VALUES (?, ?, ?, ?, ?)",
                                   (accession_number, category_name, go_id, description, score))
            conn.commit()
            conn.close()
class InterProDataLoader: # done by DAYE KWON
    """Class for loading InterPro data and storing results."""
    def __init__(self, uniprot_file_path, database_file_path):
        """Initialize InterProDataLoader object with UniProt file path and database file path."""
        self.uniprot_file_path = uniprot_file_path
        self.database_file_path = database_file_path
        self.api_url = "https://www.ebi.ac.uk/interpro/api"

    def get_accession_number_from_file(self):
        """Get UniProt accession number from file.

        Returns:
            str: UniProt accession number.
        """
        with open(self.uniprot_file_path, "r") as f:
            entry_code_line = f.readline().strip()

        entry_code = entry_code_line.split("UniProt Entry Code: ")[-1].strip()
        return entry_code

    def query_interpro_api(self, accession_number):
        """Query InterPro API for protein data.

        Args:
            accession_number (str): UniProt accession number.

        Returns:
            list: List of InterPro data.
        """
        url = f"{self.api_url}/entry/all/protein/UniProt/{accession_number}/?page_size=200&extra_fields=hierarchy,short_name"

        try:
            with urlopen(url) as res:
                data = json.loads(res.read().decode("utf-8"))
            return data["results"]
        except HTTPError as e:
            print(f"Error accessing the InterPro API: {e}")
            return []

    def insert_into_database(self, table_name, data):
        """Insert InterPro data into SQLite database.

        Args:
            table_name (str): Name of the table to insert data into.
            data (list): List of InterPro data.

        Returns:
            None
        """
        conn = sqlite3.connect(self.database_file_path)
        cursor = conn.cursor()

        insert_query = f"INSERT INTO {table_name} VALUES ({', '.join(['?' for _ in data[0]])})"
        cursor.executemany(insert_query, data)

        conn.commit()
        conn.close()

    def fetch_and_store_interpro_data(self, table_name):
        """Fetch and store InterPro data in the database.

        Args:
            table_name (str): Name of the table to store data in.

        Returns:
            None
        """
        accession_number = self.get_accession_number_from_file()
        interpro_data = self.query_interpro_api(accession_number)

        if interpro_data:
            table_data = []

            for m in interpro_data:
                meta = m["metadata"]
                if "proteins" in m:
                    protein = m["proteins"][0]
                else:
                    protein = {"accession": "-", "protein_length": 0, "entry_protein_locations": []}

                signatures = ",".join(
                    [sig for db in meta.get("member_databases", {}).values() for sig in db.keys()]) if meta.get(
                    "member_databases") else "-"
                go_terms = ",".join([t["identifier"] for t in meta["go_terms"]]) if meta.get("go_terms") else "-"

                locations = [f"{f['start']}..{f['end']}" for l in protein["entry_protein_locations"] for f in
                             l["fragments"]]

                table_data.append([
                    meta["accession"],
                    meta.get("name", "-"),
                    meta["source_database"],
                    meta["type"],
                    meta.get("integrated", "-"),
                    signatures,
                    go_terms,
                    protein["accession"].upper(),
                    str(protein["protein_length"]),
                    ",".join(locations)
                ])

            self.database_file_path = 'outputs_data.db'

            self.insert_into_database(table_name, table_data)
            print("Data inserted into the output_data database successfully.")
class PDB: # done by MORGAN GYGER
    """Class for fetching and parsing PDB data."""
    def __init__(self):
        """Initialize PDB object."""
        pass

    def fetch_top_identifiers(self):
        """Fetch top PDB identifiers from BLAST result file.

        Returns:
            list: List of top PDB identifiers.
        """
        # Get fasta from blast results
        fasta_sequences = []
        try:
            with open('flask_file/blast_result.txt', 'r') as file:
                data = file.read()
                records = data.split('\n\n')

                for record in records:
                    lines = record.split('\n')
                    for line in lines:
                        if line.startswith('Fasta Sequence:'):
                            fasta_sequence = line.split(': ')[1].strip()
                            fasta_sequences.append(fasta_sequence)
                            break
        except FileNotFoundError:
            print("Blast result file not found.")
            return []

        # Set up to only get the first, closest result from PDB
        top_identifiers = []

        # Set up search query and define
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

            # Encode the query as JSON
            encoded_query = json.dumps(query)

            # Construct the API URL for searches of interest
            api_url = "https://search.rcsb.org/rcsbsearch/v2/query?json=" + encoded_query

            # Set GET request to the API
            response = requests.get(api_url)

            # Parse JSON response
            json_response = response.json()

            # Extract the top identifier if match exists
            top_identifier = None
            if "result_set" in json_response and len(json_response["result_set"]) > 0:
                top_identifier = json_response["result_set"][0]["identifier"]

            # Append top identifier to the list of all matches for all sequences
            top_identifiers.append(top_identifier)

        for i, top_identifier in enumerate(top_identifiers):
            print(f"Top Identifier for sequence {i + 1}: {top_identifier if top_identifier else 'None'}")

        return top_identifiers

    def fetch_pdb_content(self, pdb_id):
        """Fetch PDB content for a given PDB ID.

        Args:
            pdb_id (str): PDB ID.

        Returns:
            str: PDB content.
        """
        url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an exception for 4xx or 5xx status codes
            return response.text
        except requests.RequestException as e:
            print(f"Failed to fetch PDB content for {pdb_id}. Error: {e}")
            return None

    def parse_pdb_content(self, pdb_content):
        """Parse PDB content and extract atom information.

        Args:
            pdb_content (str): PDB content.

        Returns:
            list: List of atom information.
        """
        # Create file-like object from the PDB content string
        pdb_file = StringIO(pdb_content)

        # Parse the PDB file
        parser = PDBParser()
        structure = parser.get_structure("pdb_data", pdb_file)

        atoms = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        # Convert float32 coordinates to regular floats
                        coord = [float(coord) for coord in atom.coord]
                        atoms.append({
                            "serial": atom.serial_number,
                            "name": atom.name,
                            "residue": residue.resname,
                            "chain": chain.id,
                            "resSeq": residue.id[1],
                            "x": coord[0],
                            "y": coord[1],
                            "z": coord[2],
                            "occupancy": atom.occupancy,
                            "tempFactor": atom.bfactor,
                            "element": atom.element
                        })

        return atoms

    def write_pdb(self, atoms, filename):
        """Write PDB data to a file.

        Args:
            atoms (list): List of atom information.
            filename (str): Name of the output file.

        Returns:
            None
        """
        with open(filename, 'w') as f:
            for atom in atoms:
                pdb_line = f"ATOM  {atom['serial']:>5} {atom['name']:<4} {atom['residue']:>3} {atom['chain']:1} \
{atom['resSeq']:>4}    {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}{atom['occupancy']:>6.2f}{atom['tempFactor']:>6.2f}{atom['element']:>2}\n"
                f.write(pdb_line)



def main(): # done by MORGAN GYGER, DAYE KWON, and ANDY KECK
    """Main function to execute the entire workflow."""
    conn = sqlite3.connect('outputs_data.db')
    cursor = conn.cursor()

    # Create tables for each class's output
    cursor.execute('''CREATE TABLE IF NOT EXISTS blast_results (
                        gene_name TEXT,
                        accession TEXT,
                        scientific_name TEXT,
                        e_score REAL,
                        alignment_score REAL,
                        identity REAL,
                        fasta_sequence TEXT
                    )''')

    cursor.execute('''CREATE TABLE IF NOT EXISTS sopma_results (
                        accession_number TEXT,
                        secondary_structure TEXT
                    )''')

    cursor.execute('''CREATE TABLE IF NOT EXISTS deepgo_results (
                        accession_number TEXT,
                        category TEXT,
                        go_id TEXT,
                        description TEXT,
                        score REAL
                    )''')
    cursor.execute('''CREATE TABLE IF NOT EXISTS InterPro_Data (
                         accession TEXT,
                         name TEXT,
                         source_database TEXT,
                         type TEXT,
                         integrated TEXT,
                         signatures TEXT,
                         go_terms TEXT,
                         accession_upper TEXT,
                         protein_length TEXT,
                         locations TEXT
                     )''')
    # Commit changes
    conn.commit()
    ncbi_blast = NCBI_BLAST()
    sequence = ncbi_blast.get_genomic_sequence_from_user()  # Get the genomic sequence from the user

    # Run NCBI BLAST search and retrieve top 5 results
    ncbi_blast.main(sequence)
    pdb = PDB()
    pdb.fetch_top_identifiers()
    sopma = SOPMA()
    sopma.conn = conn
    fasta_sequences = sopma.read_sequence_file('blast_result.txt')
    for accession_number, fasta_sequence in fasta_sequences:
        amino_acid_sequence = sopma.fasta_to_amino_acid(fasta_sequence)
        print("Accession Number:", accession_number)
        print("Amino Acid Sequence:", amino_acid_sequence)
        sopma.run_sopma(amino_acid_sequence, accession_number, sopma.write_to_database)

    predictor = DeepGOPredictor()
    predictor.conn = conn
    predictor.main(fasta_sequences)

    uniprot_fetcher = UniprotFetcher("blast_result.txt")
    uniprot_fetcher.fetch_and_write_uniprot_entries()

    loader = InterProDataLoader("uniprot_entry_codes.txt", "outputs_data.db")
    loader.fetch_and_store_interpro_data("InterPro_Data")

    conn.close()



if __name__ == "__main__":
    main()
