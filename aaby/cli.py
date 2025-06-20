import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="AABY: End-to-end AMBER system builder from PDB")
    parser.add_argument('-f', '--input', required=True, help='Input PDB file')
    parser.add_argument('--chains', default="A,B,C,D", help='Comma-separated list of protein chains for ACE/NME')
    parser.add_argument('--ssbond', action='store_true', help='Renumber SSBOND.txt, convert CYS to CYX, insert ssbonds into tleap.in')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--ph', type=float, help='pH for propka (protonation by predicted pKa)')
    group.add_argument('--protlist', type=str, help='Residue list file for direct protonation changes')
    parser.add_argument('--coby-yaml', type=str, help='YAML file specifying COBY arguments (use {protein_pdb} as placeholder)')
    parser.add_argument('--coby-args', nargs=argparse.REMAINDER, help='Arguments passed directly to COBY after this flag (use {protein_pdb} as placeholder)')
    return parser.parse_args()
