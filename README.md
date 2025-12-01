# AABY

Automated AMBER System Builder

## Installation

AABY can now be installed via pip, so we recommend that you use mamba/conda to create a virtual environment to install the necessary packages in which can be done in a single step using the command below.

```
micromamba env create -n AABY -f environment.yml
```



## Usage
Example:
```
python -m aaby.scripts.AABY -f INPUT.pdb --chains A --ph 7.4 --coby-yaml coby_input.yaml --antechamber True --mol2 UNK --nc 0
```
For help use
```
python -m aaby.scripts.AABY -h
```

Example of input coby file (coby_input.yaml)
```
- -box
- 15
- 15
- 15
- -membrane
- lipid:POPC:6
- lipid:CHL:4
- -protein
- file:{protein_pdb}
- charge:0
```
