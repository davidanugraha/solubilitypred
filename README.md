# SolubilityPred

This is a solubility prediction project where, based on the SMILES input, predictions of solubility in logS will be generated as output.
Multiple models were experimented with, and `XGBoost` was identified as the best performer for this dataset.
For more detailed information on how the dataset and each model were explored, please refer to `explore_data.ipynb`.

## Requirements

#### Main Requirements
- Python 3.6+
- RDKit 2023.3.1+
- Mordred 1.2.0+
- Numpy 1.21.6+
- Pandas 1.3.5+
- Joblib 1.2.0+

#### Optional: for Jupyter notebook
- Jupyter_client 7.4.9+
- Jupyter_core 4.12.0+
- Kaggle 1.5.16+
- Matplotlib 3.5.3+
- Seaborn 0.12.2+
- Opendatasets 0.1.22+
- TensorFlow 2.11.0+
- Keras 2.11.0+
- Scikit-Learn 1.0.2+
- Scipy 1.4.1+
- XGBoost 1.6.2+

## Installation

Simply do `pip install .` in the same directory as `setup.py`

## Usage

### Command line tool

After successful pip installation, a command line tool will be installed.
To generate the predictions for SMILES provided in `input_smiles.csv` and store them into `my_output.csv`, you can do:

```bash
solubilitypred input_smiles.csv my_output.csv
```

You can see all of the options available for the command line tool:

```bash
usage: solubilitypred [-h] --input INPUT [--output OUTPUT] [--cpus CPUS]

Run aqueous solubility predictor

required arguments:
  --input INPUT        Path to the file containing the SMILES. Assumes the content is 1 SMILE per line. Accepts CSV or TXT format.

optional arguments:
  -h, --help           show this help message and exit
  --output OUTPUT      Name of the output file. Defaults to output.csv. Accepts CSV or TXT format. Note this will overwrite the content of the output file.
  --cpus CPUS          Number of CPU cores to be used. Defaults to use all available cores. Must be between 1 and the number of available CPU cores.
```

### In a Python environment

`solubilitypred` also supports integration in a Python 3 environment. Example usage can be seen as followed:

```python
import pandas as pd
import solubilitypred as solpred

# Using List as an input
my_smiles = ["[H+].[H+].[Cl-].[Cl-].NNCc1ccccc1", "CC(=O)OC(C1=CC=CC=C1)C(Cl)(Cl)Cl", "CO[P](=O)(OC)OC=C(Cl)Cl", "[Na+].[Cl-]"]
predictions = solpred.predict(my_smiles)

# Using DataFrame as an input
my_smiles_df = pd.DataFrame({"SMILES": my_smiles})
predictions = solpred.predict(my_smiles)
```