import requests
from bs4 import BeautifulSoup

def run_sopma(sequence, accession_number, callback):
    # Define the SOPMA URL
    base_url = "https://npsa.lyon.inserm.fr/cgi-bin/npsa_automat.pl?page=/NPSA/npsa_sopma.html"

    # Make a GET request to fetch the SOPMA page
    response = requests.get(base_url)

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
            secondary_structure_text = extract_secondary_structure(response.text)
            # Call the callback function to handle file writing
            callback(secondary_structure_text, accession_number)
        else:
            print(f"Error: Failed to retrieve data from SOPMA for {accession_number}")
    else:
        print("Error: Failed to retrieve SOPMA page")


def extract_secondary_structure(html_content):
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





def write_to_file(data, accession_number):
    with open('sopma_result.txt', 'a') as file:
        file.write(f"Accession Number: {accession_number}\n")
        for result in data:
            file.write(result)
        file.write('\n')

def fasta_to_amino_acid(fasta_sequence):
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


# Read FASTA sequences from file
def read_sequence_file(filename):
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
                print("Accession Number:", accession_number)  # Debugging
                print("Fasta Sequence:", fasta_sequence)  # Debugging
                fasta_sequence = fasta_sequence.replace("\n", "")  # Remove newline characters
                sequences.append((accession_number, fasta_sequence))
    return sequences




# Example usage:
fasta_sequences = read_sequence_file('CLICKHERE/flask_file/blast_result.txt')
for accession_number, fasta_sequence in fasta_sequences:
    amino_acid_sequence = fasta_to_amino_acid(fasta_sequence)
    print("Accession Number:", accession_number)
    print("Amino Acid Sequence:", amino_acid_sequence)
    run_sopma(amino_acid_sequence, accession_number, write_to_file)
