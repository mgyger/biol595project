import json
from urllib.error import HTTPError
from urllib.request import urlopen
import sqlite3
# Enter the UniProt accession: P49789
def main():
    # Ask the user for the UniProt accession
    query = input("Enter the UniProt accession: ")

    api_url = "https://www.ebi.ac.uk/interpro/api"
    url = f"{api_url}/entry/all/protein/UniProt/{query}/"
    url += "?page_size=200&extra_fields=hierarchy,short_name"

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
    create_table_query = f"CREATE TABLE {table_name} ({', '.join([f'{col} TEXT' for col in columns])});"

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

if __name__ == "__main__":
    main()
