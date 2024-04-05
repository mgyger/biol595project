import requests
import re
from bs4 import BeautifulSoup

def run_sopma(sequence):
    # Define the SOPMA URL
    base_url = "https://npsa.lyon.inserm.fr"
    url = base_url + "/cgi-bin/npsa_automat.pl?page=/NPSA/npsa_sopma.html"

    # Make a GET request to fetch the SOPMA page
    response = requests.get(url)

    # Check if request was successful
    if response.status_code == 200:
        # Parse the HTML content using BeautifulSoup
        soup = BeautifulSoup(response.text, 'html.parser')

        # Find all forms on the page
        forms = soup.find_all('form')

        # Search for the form containing the paste_seq input
        for form in forms:
            if form.find('textarea', {'name': 'notice'}):
                form_action = form.get('action')
                break
        else:
            print("Error: Form for sequence input not found on SOPMA page")
            return

        # Prepare the data to be sent
        data = {
            'notice': sequence,
            'result_submit': 'Send'
        }

        # Make a POST request to the form action URL
        response = requests.post(base_url + form_action, data=data)

        # Check if request was successful
        if response.status_code == 200:
            # Extract the URL of the secondary structure table from the response URL
            secondary_structure_url = 'https://npsa.lyon.inserm.fr/cgi-bin/secpred_sopma.pl'

            # Make a GET request to fetch the secondary structure table
            response = requests.get(secondary_structure_url)

            # Check if request was successful
            if response.status_code == 200:
                # Use regular expressions to extract secondary structure composition
                match = re.search(r'Alpha helix.*?(\d+).*?(\d+\.\d+)%.*?3_10 helix.*?(\d+).*?(\d+\.\d+)%.*?Pi helix.*?(\d+).*?(\d+\.\d+)%.*?Beta bridge.*?(\d+).*?(\d+\.\d+)%.*?Extended strand.*?(\d+).*?(\d+\.\d+)%.*?Beta turn.*?(\d+).*?(\d+\.\d+)%.*?Bend region.*?(\d+).*?(\d+\.\d+)%.*?Random coil.*?(\d+).*?(\d+\.\d+)%', response.text)

                # Check if match was found
                if match:
                    # Print secondary structure composition
                    print("Alpha helix (Hh):", match.group(1), "is", match.group(2), "%")
                    print("3_10 helix (Gg):", match.group(3), "is", match.group(4), "%")
                    print("Pi helix (Ii):", match.group(5), "is", match.group(6), "%")
                    print("Beta bridge (Bb):", match.group(7), "is", match.group(8), "%")
                    print("Extended strand (Ee):", match.group(9), "is", match.group(10), "%")
                    print("Beta turn (Tt):", match.group(11), "is", match.group(12), "%")
                    print("Bend region (Ss):", match.group(13), "is", match.group(14), "%")
                    print("Random coil (Cc):", match.group(15), "is", match.group(16), "%")
                else:
                    print("Error: Failed to find secondary structure composition table on SOPMA page")
            else:
                print("Error: Failed to retrieve data from SOPMA")
        else:
            print("Error: Failed to retrieve SOPMA page")
    else:
        print("Error: Failed to retrieve SOPMA page")

# Example usage
sequence = "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG"

run_sopma(sequence)
