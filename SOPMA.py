import requests
from bs4 import BeautifulSoup
import re

def run_sopma(sequence):
    # define the SOPMA URL
    base_url = "https://npsa.lyon.inserm.fr/cgi-bin/npsa_automat.pl?page=/NPSA/npsa_sopma.html"

    # Make a POST request to fetch the SOPMA page
    response = requests.get(base_url)

    # Check if request was successful
    if response.status_code == 200:
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
            'title': '',  # optional sequence title
            'ali_width': '70'  # output width
        }

        # Make a POST request to the form action URL
        response = requests.post("https://npsa.lyon.inserm.fr" + form_action, data=data)

        # Check if request was successful
        if response.status_code == 200:
            # Extract the secondary structure information
            extract_secondary_structure(response.text)
        else:
            print("Error: Failed to retrieve data from SOPMA")
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

    # Calculate and print the ratio of each letter's proportion
    for char, count in counts.items():
        ratio = (count / total_chars) * 100
        structure_name = structure_names[char]
        print(f"{structure_name}: {count} is {ratio:.2f}%")

# Example usage
sequence = "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG"
run_sopma(sequence)


