import requests

class UniprotFetcher:
    def __init__(self, blast_result_file):
        self.blast_result_file = blast_result_file

    def fetch_uniprot_entries(self, accession_numbers):
        uniprot_entries = []

        for accession_number in accession_numbers:
            url = f"https://rest.uniprot.org/uniprotkb/search?query={accession_number}&format=tsv"
            response = requests.get(url)

            if response.status_code == 200:
                lines = response.text.split('\n')
                entries = [line.split('\t') for line in lines]
                uniprot_entries.extend(entries[1:])  # Skip the header line
            else:
                print(f"Error: {response.status_code} - {response.reason}")

        return uniprot_entries

    def extract_accession_numbers(self):
        accession_numbers = []

        with open(self.blast_result_file, "r") as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip()
                if "Accession Number:" in line:
                    accession_number = line.split(":")[1].strip()
                    accession_numbers.append(accession_number)

        return accession_numbers

    def write_uniprot_entries_to_file(self, uniprot_entries):
        with open("uniprot_entry_codes.txt", "w") as f:
            entry_codes_set = set()  # To ensure unique entry codes
            for entry in uniprot_entries:
                entry_code = entry[0]
                if entry_code not in entry_codes_set:
                    f.write("UniProt Entry Code: " + entry_code + "\n")
                    entry_codes_set.add(entry_code)

        print("UniProt entry codes written to uniprot_entry_codes.txt")

    def fetch_and_write_entries(self):
        accession_numbers = self.extract_accession_numbers()
        uniprot_entries = self.fetch_uniprot_entries(accession_numbers)
        self.write_uniprot_entries_to_file(uniprot_entries)


if __name__ == "__main__":
    uniprot_fetcher = UniprotFetcher("CLICKHERE/flask_file/blast_result.txt")
    uniprot_fetcher.fetch_and_write_entries()