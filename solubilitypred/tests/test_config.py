import os

TEST_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
TEST_DATA_DIR_PATH = os.path.join(TEST_DIR_PATH, "test_data")
TEST_INPUT_PATH = os.path.join(TEST_DATA_DIR_PATH, "test_input.csv")
TOLERANCE = 10