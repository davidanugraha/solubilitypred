# For systems
from rdkit import Chem
import os

CUR_FILE_PATH = os.path.abspath(__file__)
DIR_PATH = os.path.dirname(CUR_FILE_PATH)

RANDOM_SEED = 42

SMILES_SOLUTE = "SMILES"
SMILES_SOLVENT = "SMILES_Solvent"
TEMPERATURE = "Temperature"
SOLUBILITY = "Solubility"
OCCURRENCES = "Occurrences" # Temporary column to average them during preprocessing

DEFAULT_TEMPERATURE = 298.15
DEFAULT_SOLVENT = "O" # Water SMILES String

def is_valid_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return True
    except:
        return False
