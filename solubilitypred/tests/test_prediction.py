import numpy as np
import pytest

from .test_config import *
from ..predict import *

def test_valid_predict_dataframe():
    expected_path = os.path.join(TEST_DATA_DIR_PATH, "test_expected_predict.npy")
    expected_output = np.load(expected_path)

    input_df = pd.read_csv(TEST_INPUT_PATH)
    result = predict(input_df, num_workers=5)[SOLUBILITY].values

    # Perform the assertion to check if the results match
    np.testing.assert_almost_equal(result, expected_output, decimal=TOLERANCE)
    
def test_valid_predict_series():
    expected_path = os.path.join(TEST_DATA_DIR_PATH, "test_expected_predict.npy")
    expected_output = np.load(expected_path)

    input_df = pd.read_csv(TEST_INPUT_PATH)
    input_df = input_df.rename(columns={SMILES_SOLUTE: "RandomName"})
    input_series = input_df["RandomName"]
    result = predict(input_series, num_workers=5)[SOLUBILITY].values

    # Perform the assertion to check if the results match
    np.testing.assert_almost_equal(result, expected_output, decimal=TOLERANCE)
    
def test_valid_predict_list():
    expected_path = os.path.join(TEST_DATA_DIR_PATH, "test_expected_predict.npy")
    expected_output = np.load(expected_path)

    input_df = pd.read_csv(TEST_INPUT_PATH)
    input_list = input_df[SMILES_SOLUTE].to_list()
    result = predict(input_list, num_workers=5)[SOLUBILITY].values

    # Perform the assertion to check if the results match
    np.testing.assert_almost_equal(result, expected_output, decimal=TOLERANCE)

def test_invalid_predict():
    # Test with invalid input, should raise ValueError
    with pytest.raises(ValueError, match="Input should be either a list of SMILES, a Series of SMILES, or a DataFrame with a single column of SMILES."):
        predict(None, num_workers=5)
        
    # Test with invalid DataFrame, should raise ValueError
    invalid_df = pd.DataFrame({SMILES_SOLUTE: ["CO", "O=O"], "Info": [1, 2]})
    with pytest.raises(ValueError, match="DataFrame should have only one column for SMILES."):
        predict(invalid_df, num_workers=5)

    # Test with invalid SMILES, should raise ValueError
    invalid_smiles_list = ["invalid1", "invalid2", "invalid3"]
    with pytest.raises(ValueError, match=r"Found invalid SMILES string `invalid\d`."):
        predict(invalid_smiles_list, num_workers=5)
    
if __name__ == '__main__':
    test_valid_predict_dataframe()
    test_valid_predict_series()
    test_valid_predict_list()
    test_invalid_predict()