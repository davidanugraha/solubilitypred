# For systems
import os
from pathlib import Path
import joblib
import json
import inspect
from concurrent.futures import ThreadPoolExecutor, as_completed

# For data processing
import pandas as pd

# For cheminformatics library
from rdkit import Chem
from rdkit.Chem import Descriptors
from mordred import Calculator, descriptors

# For files
from .utils import *

def _load_features():
    features_json_path = os.path.join(DIR_PATH, "features.json")
    try:
        with open(features_json_path, "r") as f:
            features = json.load(f)
        return features
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: '{features_json_path}' file not found. Please ensure the file is present before running the program.")
    except json.JSONDecodeError:
        raise ValueError(f"Error: Unable to load '{features_json_path}'. Please make sure the file contains valid JSON data.")
    
def _calculate_property_rdkit(smiles, property_name):
    molecule = Chem.MolFromSmiles(smiles)
    property_value = getattr(Descriptors, property_name)(molecule)
    return property_value

def _calculate_all_properties_rdkit(smiles, properties):
    return {prop: _calculate_property_rdkit(smiles, prop) for prop in properties}

def _remove_duplicate_features(df):
    # Remove duplicated columns (potentially from RDKit and Mordred overlap)
    duplicate_columns = df.columns[df.columns.duplicated()]
    if len(duplicate_columns) > 0:
        # Create a new DataFrame with only the duplicated columns
        duplicated_df = df[duplicate_columns]

        # Check if the values in the duplicated columns are the same for each row
        if not duplicated_df.apply(lambda col: col.nunique(), axis=1).eq(1).all():
            # Resolve by averaging for each set of duplicated columns
            averaged_columns = {col_name: duplicated_df[col_name].mean(axis=1) for col_name in duplicated_df}
            df = df.drop(columns=duplicate_columns)
            for col_name, averaged_values in averaged_columns.items():
                df[col_name] = averaged_values
    
    return df

def _add_selected_features(df, features, num_workers):
    rdkit_features = []
    
    # Find RDKit features, if can't find the attr_name, then it must be a mordred feature
    for attr_name in features:
        for _, module in inspect.getmembers(Chem, inspect.ismodule):
            if hasattr(module, attr_name):
                rdkit_features.append(attr_name)
                break
    
    # Calculate RDKit properties
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        rdkit_result = list(executor.map(_calculate_all_properties_rdkit, df[SMILES_SOLUTE], [rdkit_features] * len(df)))
    rdkit_df = pd.DataFrame(rdkit_result, columns=rdkit_features)
    df = pd.concat([df, rdkit_df], axis=1, join='outer')
    
    # Calculate mordred properties and filter them
    calc = Calculator(descriptors)
    smiles_data = [Chem.MolFromSmiles(smi) for smi in df[SMILES_SOLUTE].to_list()[:]]
    mordred_df = calc.pandas(smiles_data, nproc=num_workers, quiet=True)
    mordred_features = [feat for feat in features if feat in mordred_df.columns]
    mordred_df = mordred_df[mordred_features]
    df = pd.concat([df, mordred_df], axis=1, join='outer')
    
    # Remove duplicated features from RDKit and mordred (potentially)
    df = _remove_duplicate_features(df)
    
    # Reorder based on features
    df = df.reindex(columns=([SMILES_SOLUTE] + features))
    
    return df

def _scale_the_features(df):
    scaler_path = os.path.join(DIR_PATH, "data_scaler.pkl")
    if not os.path.exists(scaler_path):
        raise FileNotFoundError(f"Scaler file not found: {scaler_path}")
    try:
        data_scaler = joblib.load(scaler_path)
        final_input = data_scaler.transform(df)
        return final_input
    except Exception as e:
        raise ValueError(f"Failed to load and use the scaler: {str(e)}")

def preprocess_input(input_df, num_workers):
    try:
        # Get features to be used
        features = _load_features()
        
        # Add selected features through RDKit and mordred
        modified_df = _add_selected_features(input_df, features, num_workers)
        
        # Drop SMILES
        modified_df = modified_df.drop(columns=[SMILES_SOLUTE])
        
        # Scale the features
        final_input = _scale_the_features(modified_df)

        return final_input
    except Exception as e:
        raise RuntimeError(f"Failed to preprocess the input: {str(e)}")
