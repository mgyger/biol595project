import json
from urllib.error import HTTPError
from urllib.request import urlopen
import sqlite3

def get_accession_number_from_file(file_path):
    # Open the file and read the accession number
    with open(file_path, "r") as f:
        entry_code_line = f.readline().strip()  # Read the first line

    # Extract the entry code from the line
    entry_code = entry_code_line.split("UniProt Entry Code: ")[-1].strip()

    print("Extracted UniProt Entry Code:", entry_code)
    return entry_code

def main():
    # Read the accession number from the file
    accession_number = get_accession_number_from_file("uniprot_entry_codes.txt")

    # Query the InterPro API with the accession number
    api_url = "https://www.ebi.ac.uk/interpro/api"
    url = f"{api_url}/entry/all/protein/UniProt/{accession_number}/?page_size=200&extra_fields=hierarchy,short_name"

    try:
        with urlopen(url) as res:
            data = json.loads(res.read().decode("utf-8"))

        # Initialize a list to store the data to be inserted into the SQL table
        table_data = []

        for i, m in enumerate(data["results"]):
            meta = m["metadata"]
            protein = m["proteins"][0]

            if meta["member_databases"]:
                dbs = meta["member_databases"].values()
                signatures = ",".join([sig for db in dbs for sig in db.keys()])
            else:
                signatures = "-"

            if meta["go_terms"]:
                go_terms = ",".join([t["identifier"] for t in meta["go_terms"]])
            else:
                go_terms = "-"

            locations = []
            for l in protein["entry_protein_locations"]:
                for f in l["fragments"]:
                    locations.append(f"{f['start']}..{f['end']}")

            table_data.append([
                meta["accession"],
                meta["name"] or "-",
                meta["source_database"],
                meta["type"],
                meta["integrated"] or "-",
                signatures,
                go_terms,
                protein["accession"].upper(),
                str(protein["protein_length"]),
                ",".join(locations)
            ])

        # Create SQL table schema
        table_name = "InterPro_Data"
        columns = ["accession", "name", "source_database", "type", "integrated", "signatures", "go_terms", "protein_accession", "protein_length", "locations"]
        create_table_query = f"CREATE TABLE IF NOT EXISTS {table_name} ({', '.join([f'{col} TEXT' for col in columns])});"

        # Connect to SQLite3 database
        conn = sqlite3.connect('interpro_data.db')
        cursor = conn.cursor()

        # Execute the SQL command to create the table
        cursor.execute(create_table_query)

        # Insert data into the created table
        insert_query = f"INSERT INTO {table_name} VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
        cursor.executemany(insert_query, table_data)

        # Commit the changes and close the connection
        conn.commit()
        conn.close()

        print("SQL table created successfully.")

    except HTTPError as e:
        print(f"Error accessing the InterPro API: {e}")

if __name__ == "__main__":
    main()
