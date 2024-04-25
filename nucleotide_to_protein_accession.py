import requests


def fetch_uniprot_entries(nucleotide_accession):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={nucleotide_accession}&format=tsv"

    # Send HTTP GET request to UniProt REST API
    response = requests.get(url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Parse the TSV response
        lines = response.text.split('\n')
        entries = [line.split('\t') for line in lines]
        return entries
    else:
        # Print error message if request was not successful
        print(f"Error: {response.status_code} - {response.reason}")
        return None


if __name__ == "__main__":
    nucleotide_accession = input("Enter the nucleotide accession number: ")

    uniprot_entries = fetch_uniprot_entries(nucleotide_accession)
    if uniprot_entries:
        # Write the UniProt entry codes to a text file
        with open("uniprot_entry_codes.txt", "w") as f:
            for i, entry in enumerate(uniprot_entries):
                if i == 0 or i == 2:  # Skip the first and third lines
                    continue
                entry_code = entry[0]  # UniProt ID is typically the first column
                f.write("UniProt Entry Code: " + entry_code + "\n")
        print("UniProt entry codes written to uniprot_entry_codes.txt")