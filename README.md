# AABY

Automated AMBER System Builder Workflow

## Installation

   git clone https://github.com/AmandaStange/AABY.git
   cd AABY
   micromamba env create -n AABY -f environment.yml

## Usage
    Example:

    python Scripts/AABY.py -f INPUT.pdb --chains A --ph 7.4 --coby-yaml coby_input.yaml --antechamber True --mol2 UNK --nc 0

