import requests
import sqlite3

# API endpoint
url = 'https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create'

# Prompt the user for the protein sequence
user_sequence = input("Please enter the protein sequence in fasta format: ")

# Data payload with user's input
data = {
    "version": "1.0.18",
    "data_format": "fasta",
    "data": user_sequence,
    "threshold": 0.4
}

def create_sql_table():
    try:
        # Connect to SQLite3 database
        conn = sqlite3.connect('DeepGoWeb.db')
        cursor = conn.cursor()

        # Create SQL table schema
        table_name = "ProteinPredictions"

        columns = ["cellular_component TEXT", "molecular_function TEXT", "biological_function TEXT"]
        create_table_query = f"CREATE TABLE IF NOT EXISTS {table_name} ({', '.join(columns)});"

        # Execute the SQL command to create the table
        cursor.execute(create_table_query)

        # Commit the changes and close the connection
        conn.commit()
        conn.close()

        print("SQL table created successfully.")

    except Exception as e:
        print(f"Error creating SQL table: {e}")

def insert_predictions_into_sql(predictions):
    try:
        # Connect to SQLite3 database
        conn = sqlite3.connect('DeepGoWeb.db')
        cursor = conn.cursor()

        # Insert predictions into the created table
        insert_query = "INSERT INTO ProteinPredictions (cellular_component, molecular_function, biological_function) VALUES (?, ?, ?)"
        cursor.executemany(insert_query, predictions)

        # Commit the changes and close the connection
        conn.commit()
        conn.close()

        print("Predictions inserted into SQL table successfully.")

    except Exception as e:
        print(f"Error inserting predictions into SQL table: {e}")

# Headers specify that the payload is JSON
headers = {
    'Content-Type': 'application/json'
}

# Send the POST request
response = requests.post(url, json=data, headers=headers)

# Check the response and parse it if successful
if response.status_code == 200 or response.status_code == 201:
    result = response.json()
    predictions = result['predictions'][0]['functions']

    # Prepare predictions for insertion into SQL table
    prediction_data = [(func[1], func[0], func[2]) for category in predictions for func in category['functions']]

    # Print formatted output
    for category in predictions:
        print(category['name'])
        for func in category['functions']:
            go_id, description, score = func
            print(f"{go_id} - {description} - {round(score, 3)}")
        print()  # Blank line for visual separation

    # Call the function to insert predictions into the SQL table
    insert_predictions_into_sql(prediction_data)
else:
    print('Failed to retrieve data:', response.status_code)
    print(response.text)

if __name__ == "__main__":
    # Call the function to create the SQL table
    create_sql_table()
