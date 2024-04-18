import requests

# API endpoint
url = 'https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create'

# Prompt the user for the protein sequence
user_sequence = input("Please enter the protein sequence in fasta format: ")

# Data payload with user's input
data = {
    "version": "1.0.18",
    "data_format": "fasta", # must be small letters
    "data": user_sequence,
    "threshold": 0.4
}

# Headers specify that the payload is JSON
# This is required for the server to understand the data
headers = {
    'Content-Type': 'application/json'
}

# Send the POST request
# json=data will automatically convert the data dictionary to JSON
response = requests.post(url, json=data, headers=headers)

# Check the response and parse it if successful
# 200 and 201 are the HTTP status codes for success
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
