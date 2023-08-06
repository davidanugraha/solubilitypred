from setuptools import setup, find_packages
from pathlib import Path

dir_path = Path('__file__').parent.absolute()
md_path = os.path.join(dir_path, "README.md")
long_description = Path(md_path).read_text(encoding='utf-8')

setup(
    name="solubilitypred",
    version="1.0.0",
    description="A solubility prediction model in aqueous solution based on SMILES",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/davidanugraha/solubilitypred",
    author="David Anugraha",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering :: Chemistry :: Chemoinformatics",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    keywords="ML, solubility",
    packages=find_packages(),
    python_requires=">=3.6, <4",
    install_requires=[
        "joblib>=1.2.0",
        "jupyter_client>=7.4.9",
        "jupyter_core>=4.12.0",
        "kaggle>=1.5.16",
        "keras==2.11.0",
        "matplotlib>=3.5.3",
        "mordred>=1.2.0",
        "numpy>=1.21.6",
        "opendatasets>=0.1.22",
        "pandas>=1.3.5",
        "rdkit>=2023.3.1",
        "scikit-learn==1.0.2",
        "scipy>=1.4.1",
        "seaborn>=0.12.2",
        "tensorflow==2.11.0",
        "xgboost==1.6.2",
    ],
    package_data={
        "solubilitypred":["solubilitypred_model.pkl", "data_scaler.pkl", "features.json"],
    },
    entry_points={
        "console_scripts":["solubilitypred=solubilitypred.main:main"],
    },
    project_urls={
        "Bug Reports": "https://github.com/davidanugraha/solubilitypred/issues",
        "Source": "https://github.com/davidanugraha/solubilitypred",
    },
)