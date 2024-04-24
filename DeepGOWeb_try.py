# Requires Python 3.6 to install DeepGOWeb

import subprocess
import json

def predict_protein_functions(sequence):
    # Execute DeepGOPlus command with the input sequence
    command = f"echo '{sequence}' | deepgoplus"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Decode the byte string to text string
    stdout = result.stdout.decode('utf-8')

    # Parse DeepGOPlus output to extract predictions
    predictions = parse_predictions(stdout)

    # Return predictions
    return predictions

def parse_predictions(output):
    try:
        # Parse the JSON output
        data = json.loads(output)

        # Extract relevant information
        predictions = {
            "uuid": data["uuid"],
            "version": data["version"],
            "predictions": []
        }

        # Iterate over predictions
        for prediction in data["predictions"]:
            sequence = prediction["sequence"]
            functions = []

            # Iterate over function categories
            for category in prediction["functions"]:
                category_name = category["name"]
                category_functions = []

                # Iterate over functions within the category
                for func in category["functions"]:
                    go_id, go_name, confidence = func
                    category_functions.append({"id": go_id, "name": go_name, "confidence": confidence})

                # Add category if there are functions
                if category_functions:
                    functions.append({"name": category_name, "functions": category_functions})

            # Add prediction if there are functions
            if functions:
                predictions["predictions"].append({"sequence": sequence, "functions": functions})

        return predictions

    except json.JSONDecodeError as e:
        print("Error:", e)
        return {}

def save_to_json(data, filename):
    # Save data to a JSON file
    with open(filename, 'w') as file:
        json.dump(data, file, indent=4)

    print(f"DeepGOPlus predictions saved to '{filename}' in JSON format")

if __name__ == "__main__":
    # Get protein sequence input from user
    sequence = input("Enter the protein sequence: ")

    # Predict protein functions using DeepGOPlus
    predictions = predict_protein_functions(sequence)

    # Save predictions to a JSON file
    save_to_json(predictions, "deepgoplus_predictions.json")
