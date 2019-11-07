# SYBA
SYnthetic BAyesian classifier is Python package for binary classification to easy- and hard-to-synthesize structures.
SYBA determines class by the difference in the occurences of ECFP-like fragments in traning compounds from both classes.
Precomputed count file is provided within this package however can be recalculated. Workflow is provided in the form of
Jupyter notebooks and can be used even for other classification problems with different training sets.

## Instalation
### Prerequisities
#### Supported platforms:
* All platforms

#### Dependencies
* RDKit

### Installation with Anaconda
SYBA is distributed as a conda package. At the moment, this is the preferred way to install and use the library.
All you need to do is get the full Anaconda distribution or its lightweight variant, Miniconda. It is essentially a
Python distribution, package manager and virtual environment in one and makes setting up a development environment
for any project very easy. After installing Anaconda/Miniconda (and environment preparing) you can run the following in
the Linux terminal (probably in your environment):
```bash
conda install -c rdkit -c lich syba
```

### Installation with setup.py
Anyway you have installed RDKit, you can download/clone Syba and install it 
from its directory with:
```bash
python setup.py install
```

## Quick start
SYBA can be applied
Input could be a file or STDIN and should be in a csv-like format ID,SMILES,OTHER_COLUMNS. Other columns are skipped.
Output will be generated in the format ID,SMILES,SYBA_SCORE. As is defined, SYBA score shows its confidancy of prediction
to one or the other class. In default behaivior, negative values mean hard-to-synthesize and positive easy-to-synthesize.

To run the script just write: 

```bash
python -m syba.syba [INPUT_FILE [OUTPUT_FILE]]
```
## Use in Python script
### Basic usage
```python
from rdkit import Chem
from syba import SybaClassifier

syba = SybaClassifier()
syba.fitDefaultScore()
smi = "O=C(C)Oc1ccccc1C(=O)O"
syba.predict(smi)
# syba works also with RDKit RDMol objects
mol = Chem.MolFromSmiles(smi)
syba.predict(mol=mol)
# syba.predict is actually method with two keyword parameters "smi" and "mol", if both provided score is calculated for compound defined in "smi" parameter has the priority
syba.predict(smi=smi, mol=mol)
```

## SYBA workflow
How were fragment scores prepared for SYBA is accessible in `docs\notebooks\prepare_fragment_counts.ipynb` and calculation of SYBA score in Python script altogether with SAScore and SCScore is in `docs\notebooks\prepare_results.ipynb`.
To run this notebooks, you will need to install Jupyter, preferably with conda `conda install jupyter`.

