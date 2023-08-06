import os
import joblib
import pandas as pd
import multiprocessing

from .utils import *
from .data_preprocess import preprocess_input

# Predict based on input_df (note input_df has already been validated)
def predict_validated(input_df, num_workers=multiprocessing.cpu_count()):
    # Preprocess the input first
    try:
        modified_input = preprocess_input(input_df, num_workers)
    except Exception as e:
        raise RuntimeError(f"Failed to preprocess the input SMILES: {str(e)}")

    # Load the model and do predictions
    model_path = os.path.join(DIR_PATH, "solubilitypred_model.pkl")
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")
    try:
        model = joblib.load(model_path)
    except Exception as e:
        raise ValueError(f"Failed to load the model: {str(e)}")
    
    # Perform predictions
    try:
        predictions = model.predict(modified_input)
    except Exception as e:
        raise RuntimeError(f"Failed to make predictions: {str(e)}")

    # Concat input columns with solubility output
    predictions_df = pd.concat([input_df, pd.Series(predictions, name=SOLUBILITY)], axis=1)
    return predictions_df

# Predict based on input (can only accept List, Series, or DataFrame)
def predict(input, num_workers=multiprocessing.cpu_count()):
    # Convert input to DataFrame
    if isinstance(input, list):
        input_df = pd.DataFrame({SMILES_SOLUTE: input})
    elif isinstance(input, pd.DataFrame):
        if len(input.columns) != 1:
            raise ValueError("DataFrame should have only one column for SMILES.")
        input_df = input.rename(columns={input.columns[0]: SMILES_SOLUTE})
    elif isinstance(input, pd.Series):
        input_df = input.to_frame()
        input_df = input_df.rename(columns={input_df.columns[0]: SMILES_SOLUTE})
    else:
        raise ValueError("Input should be either a list of SMILES, a Series of SMILES, or a DataFrame with a single column of SMILES.")
    
    # Do Validation to ensure that SMILES strings are valid
    for smiles in input_df[SMILES_SOLUTE]:
        if not is_valid_smiles(smiles):
            raise ValueError(f"Found invalid SMILES string `{smiles}`.")
    
    return predict_validated(input_df, num_workers)