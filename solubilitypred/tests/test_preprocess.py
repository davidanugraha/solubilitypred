import numpy as np
import pandas as pd
import pytest

from .test_config import *
from ..data_preprocess import preprocess_input

def test_preprocess_input():
    expected_path = os.path.join(TEST_DATA_DIR_PATH, "test_expected_preprocess.npy")
    expected_output = np.load(expected_path)

    input_df = pd.read_csv(TEST_INPUT_PATH)
    result = preprocess_input(input_df, num_workers=5)

    # Perform the assertion to check if the results match
    np.testing.assert_almost_equal(result, expected_output, decimal=TOLERANCE)

if __name__ == '__main__':
    test_preprocess_input()