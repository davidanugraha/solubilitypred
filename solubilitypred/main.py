import argparse
import csv
import multiprocessing
import os

import pandas as pd

from .utils import *
from .predict import predict_validated

def _validate_and_prepare_input_file(input_file):
    data = []
    with open(input_file, 'r') as file:
        if input_file.endswith('.csv'):
            # Assume SMILES is in the first column
            csv_reader = csv.reader(file)

            # Read the first row (header) or return None if empty file
            header_row = next(csv_reader, None)
            if header_row is not None and is_valid_smiles(header_row[0]):
                data.append((smiles, ))

            # Start from line 2
            for line_num, row in enumerate(csv_reader, start=2):  
                smiles = row[0].strip()
                if is_valid_smiles(smiles):
                    data.append((smiles, ))
                else:
                    raise ValueError(f"Invalid SMILES at line {line_num}: {smiles}")
        elif input_file.endswith('.txt'):  # Assume TXT file with 1 SMILES per line
            for line_num, line in enumerate(file, start=1):
                smiles = line.strip()
                if is_valid_smiles(smiles):
                    data.append((smiles, ))
                else:
                    raise ValueError(f"Invalid SMILES at line {line_num}: {smiles}")
        else:
            raise ValueError(f"File type is not supported! Current input file: {input_file}")

    df = pd.DataFrame(data, columns=[SMILES_SOLUTE])
    return df

def _validate_output_file(output_file):
    # File needs to end with CSV or TXT format
    if output_file is not None and not output_file.endswith('.csv') and not output_file.endswith('.txt'):
        raise ValueError(f"File type is not supported! Current output file: {output_file}")

def _parse_args():
    parser = argparse.ArgumentParser(description="Run aqueous solubility predictor")

    parser.add_argument("--input",
                        type=str,
                        required=True,
                        help="Path to the file containing the SMILES. Assumes the content is 1 SMILE per line. Accepts CSV or TXT format.")
    parser.add_argument("--output",
                        type=str,
                        default="output.csv",
                        help="Name of the output file. Defaults to output.csv in current directory. Accepts CSV or TXT format. Note this will overwrite the content of the output file.")
    parser.add_argument('--cpus',
                        default=multiprocessing.cpu_count(),
                        type=int,
                        help="Number of CPU cores to be used. Defaults to use all available cores.",
                        choices=range(1, multiprocessing.cpu_count() + 1))

    args = parser.parse_args()

    return args

# Main function, the entry point
def main():
    args = _parse_args()
    
    # Validate the input file
    try:
        df = _validate_and_prepare_input_file(args.input)
    except ValueError as e:
        print(str(e))
        return None
    
    # Validate the output (if there is any)
    try:
        _validate_output_file(args.output)
    except ValueError as e:
        print(str(e))
        return None
    
    # Perform prediction
    try:
        predictions = predict_validated(df, num_workers=args.cpus)
    except Exception as e:
        print(str(e))
        return None
    
    # Write predictions to output file
    if args.output.endswith(".csv"):
        print(f"Writing prediction to CSV file: {args.output}.")
        predictions.to_csv(args.output, index=False)
    elif args.output.endswith(".txt"):
        print(f"Writing prediction to TXT file: {args.output}.")
        predictions[SMILES_SOLUTE].to_csv(args.output, index=False, header=False)
    else:
        # Shouldn't go here, we have validated the output earlier
        print(f"Invalid output file extension. Please use a CSV or TXT extension for the output file. Output file was: {args.output}.")
        
    print("Solubility prediction has been completed successfully.")
        
if __name__ == '__main__':
    main()