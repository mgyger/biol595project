import requests
from bs4 import BeautifulSoup
from Bio.Blast import NCBIWWW, NCBIXML
class SOPMA:
    def __init__(self):
        self.base_url = "https://npsa.lyon.inserm.fr/cgi-bin/npsa_automat.pl?page=/NPSA/npsa_sopma.html"

    def run_sopma(self, sequence, accession_number, callback):
        # Make a GET request to fetch the SOPMA page
        response = requests.get(self.base_url)

        # Check if request was successful
        if response.status_code == 200:
            print("SOPMA page retrieved successfully.")
            # Parse the HTML content using BeautifulSoup
            soup = BeautifulSoup(response.text, 'html.parser')

            # Find the form on the SOPMA page
            form = soup.find('form')

            # Extract the action attribute of the form
            form_action = form.get('action')

            # Prepare the data to be sent
            data = {
                'notice': sequence,
                'states': '3',  # assuming default parameters for simplicity
                'width': '17',  # assuming default parameters for simplicity
                'threshold': '8',  # assuming default parameters for simplicity
                'title': accession_number,  # using accession number as the title
                'ali_width': '70'  # output width
            }

            # Make a POST request to the form action URL
            response = requests.post("https://npsa.lyon.inserm.fr" + form_action, data=data)

            # Check if request was successful
            if response.status_code == 200:
                print("SOPMA analysis successful.")
                # Extract the secondary structure information
                secondary_structure_text = self.extract_secondary_structure(response.text)
                # Call the callback function to handle file writing
                callback(secondary_structure_text, accession_number)
            else:
                print(f"Error: Failed to retrieve data from SOPMA for {accession_number}")
        else:
            print("Error: Failed to retrieve SOPMA page")

    def extract_secondary_structure(self, html_content):
        # Mapping lowercase letters to their corresponding names in the secondary structure table
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

        # Parse HTML content using BeautifulSoup
        soup = BeautifulSoup(html_content, 'html.parser')

        # Find the <PRE> tag containing secondary structure information
        pre_tag = soup.find('pre')

        # Extract text from <PRE> tag
        secondary_structure_text = pre_tag.get_text()

        # Count occurrences of lowercase letters
        counts = {}
        total_chars = 0
        for char in secondary_structure_text:
            if char.islower():
                counts[char] = counts.get(char, 0) + 1
                total_chars += 1

        # Calculate and return the ratio of each letter's proportion
        results = []
        for char, count in counts.items():
            ratio = (count / total_chars) * 100
            structure_name = structure_names[char]
            results.append(f"{structure_name}: {count} is {ratio:.2f}%\n")

        return results

    def write_to_file(self, data, accession_number):
        with open('sopma_result.txt', 'a') as file:
            file.write(f"Accession Number: {accession_number}\n")
            for result in data:
                file.write(result)
            file.write('\n')

    def fasta_to_amino_acid(self, fasta_sequence):
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

        # Translate DNA or RNA sequence into amino acids
        amino_acids = ''
        for i in range(0, len(fasta_sequence), 3):  # Iterate over the sequence with step 3
            codon = fasta_sequence[i:i+3]  # Get the codon (3 characters)
            if codon in codon_table:
                amino_acid = codon_table[codon]
                amino_acids += amino_acid

        return amino_acids

    def read_sequence_file(self, filename):
        sequences = []
        with open(filename, 'r') as file:
            accession_number = ''
            fasta_sequence = ''
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace
                if line.startswith("Accession Number:"):
                    accession_number = line.split(": ")[1]
                elif line.startswith("Fasta Sequence:"):
                    fasta_sequence = line.split(": ")[1]
                    fasta_sequence = fasta_sequence.replace("\n", "")  # Remove newline characters
                    sequences.append((accession_number, fasta_sequence))
        return sequences

class NCBI_BLAST:
    def run_ncbi_blast(self, sequence):
        # Perform BLAST search on NCBI
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        blast_records = NCBIXML.parse(result_handle)

        top_matches = []
        for record in blast_records:
            # Iterate through each BLAST record
            for alignment in record.alignments[:5]:  # Limit to top 5 results
                # Extract information from alignment
                gene_name = alignment.title.split("|")[4].split()[0]  # Extract gene name
                accession = alignment.accession  # Extract accession number
                scientific_name = alignment.hit_def.split("[")[0].strip()  # Extract scientific name
                e_score = alignment.hsps[0].expect  # Extract E score
                alignment_score = alignment.hsps[0].score  # Extract alignment score
                identity = alignment.hsps[0].identities / alignment.hsps[0].align_length * 100  # Calculate identity (%)
                fasta_sequence = alignment.hsps[0].sbjct  # Extract fasta sequence

                # Append data to top_matches list
                top_matches.append(
                    (gene_name, accession, scientific_name, e_score, alignment_score, identity, fasta_sequence))

        return top_matches

    def get_genomic_sequence_from_user(self):
        sequence = input("Enter the genomic sequence: ")
        return sequence

    def main(self, sequence):  # Adjusted to accept sequence as argument
        # Run NCBI BLAST search and retrieve top 5 results
        top_results = self.run_ncbi_blast(sequence)

        # Write data to a text file
        with open("blast_result.txt", "w") as file:
            # Write user input sequence at the top
            file.write("User Input Sequence:\n")
            file.write(sequence + "\n\n")

            # Write BLAST results
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

        print("Data written to blast_result.txt.")



class DeepGOPredictor:
    def __init__(self):
        self.url = 'https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create'
        self.headers = {'Content-Type': 'application/json'}

    def predict_go_functions(self, sequence, threshold=0.4):
        data = {
            "version": "1.0.18",
            "data_format": "fasta",  # must be small letters
            "data": sequence,
            "threshold": threshold
        }

        # Send the POST request
        response = requests.post(self.url, json=data, headers=self.headers)

        # Check the response and parse it if successful
        if response.status_code == 200 or response.status_code == 201:
            result = response.json()
            predictions = result['predictions'][0]['functions']

            # Print formatted output
            for category in predictions:
                print(category['name'])
                for func in category['functions']:
                    go_id, description, score = func
                    print(f"{go_id} - {description} - {round(score, 3)}")
                print()  # Blank line for visual separation
        else:
            print('Failed to retrieve data:', response.status_code)
            print(response.text)

    def main(self, sequence):
        # Predict GO functions for the protein sequence
        self.predict_go_functions(sequence)


def main():
    ncbi_blast = NCBI_BLAST()
    sequence = ncbi_blast.get_genomic_sequence_from_user()  # Get the genomic sequence from the user

    # Run NCBI BLAST search and retrieve top 5 results
    ncbi_blast.main(sequence)

    sopma = SOPMA()
    fasta_sequences = sopma.read_sequence_file('blast_result.txt')
    for accession_number, fasta_sequence in fasta_sequences:
        amino_acid_sequence = sopma.fasta_to_amino_acid(fasta_sequence)
        print("Accession Number:", accession_number)
        print("Amino Acid Sequence:", amino_acid_sequence)
        sopma.run_sopma(amino_acid_sequence, accession_number, sopma.write_to_file)

    predictor = DeepGOPredictor()
    predictor.main(sequence)  # Use the same sequence for DeepGO prediction


if __name__ == "__main__":
    main()