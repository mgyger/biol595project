from flask import Flask, render_template
import sqlite3

app = Flask(__name__)

def create_sql_table():
    try:
        # Connect to SQLite3 database
        conn = sqlite3.connect('interpro_data.db')
        cursor = conn.cursor()

        # Create SQL table schema
        table_name = "InterPro_Data"

        columns = ["Domain_accession", "name", "source_database", "type", "integrated", "signatures", "go_terms",
                   "protein_accession", "protein_length", "locations"]
        create_table_query = f"CREATE TABLE IF NOT EXISTS {table_name} ({', '.join([f'{col} TEXT' for col in columns])});"

        # Execute the SQL command to create the table
        cursor.execute(create_table_query)

        # Commit the changes and close the connection
        conn.commit()
        conn.close()

        print("SQL table created successfully.")

    except Exception as e:
        print(f"Error creating SQL table: {e}")

def insert_data_into_sql_table(table_data):
    try:
        # Connect to SQLite3 database
        conn = sqlite3.connect('interpro_data.db')
        cursor = conn.cursor()

        # Insert data into the created table
        insert_query = f"INSERT INTO InterPro_Data VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
        cursor.executemany(insert_query, table_data)

        # Commit the changes and close the connection
        conn.commit()
        conn.close()

        print("Data inserted into SQL table successfully.")

    except Exception as e:
        print(f"Error inserting data into SQL table: {e}")

def generate_html_table_from_sqlite(database_file, output_file):
    try:
        # Connect to SQLite database
        conn = sqlite3.connect(database_file)
        cursor = conn.cursor()

        # Fetch data from the SQL table
        cursor.execute("SELECT * FROM InterPro_Data")
        table_data = cursor.fetchall()

        # Close the connection
        conn.close()

        # Generate HTML table
        html_table = "<!DOCTYPE html>\n<html lang='en'>\n<head>\n<meta charset='UTF-8'>\n<title>SQL Table</title>\n<style>\ntable {\nborder-collapse: collapse;\nwidth: 100%;\n}\nth, td {\nborder: 1px solid #dddddd;\ntext-align: left;\npadding: 8px;\n}\nth {\nbackground-color: #f2f2f2;\n}\n</style>\n</head>\n<body>\n<h1>SQL Table</h1>\n<table>\n<thead>\n<tr>\n<th>Protein Accession</th>\n<th>Domain Accession</th>\n<th>Name</th>\n<th>Source Database</th>\n<th>Type</th>\n<th>Integrated</th>\n<th>Signatures</th>\n<th>GO Terms</th>\n<th>Protein Length</th>\n<th>Locations</th>\n</tr>\n</thead>\n<tbody>\n"

        for row in table_data:
            html_table += "<tr>"
            for col in row:
                html_table += f"<td>{col}</td>"
            html_table += "</tr>\n"

        html_table += "</tbody>\n</table>\n</body>\n</html>"

        # Write HTML table to output file
        with open(output_file, "w") as f:
            f.write(html_table)

        print(f"HTML table successfully created and saved to {output_file}")

    except Exception as e:
        print(f"Error generating HTML table: {e}")

@app.route('/')
def display_table():
    try:
        # Connect to SQLite3 database
        conn = sqlite3.connect('interpro_data.db')
        cursor = conn.cursor()

        # Fetch data from the SQL table
        cursor.execute("SELECT * FROM InterPro_Data")
        table_data = cursor.fetchall()

        # Close the connection
        conn.close()

        # Render HTML template with table data
        return render_template('table.html', table_data=table_data)

    except Exception as e:
        return f"Error displaying table: {e}"

if __name__ == "__main__":
    create_sql_table()

    # Call the function to generate the HTML table from the SQLite database
    generate_html_table_from_sqlite("interpro_data.db", "templates/table.html")

    # Run the Flask application
    app.run(debug=True)
