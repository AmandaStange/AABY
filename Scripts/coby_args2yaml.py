#!/usr/bin/env python3

import argparse
import shlex
import yaml
from pathlib import Path

#python coby_args2yaml.py "-box 16 16 25 -membrane lipid:POPC:6:params:Amber center:0:0:-7 -protein file:{protein_pdb} charge:-31" -o my_coby.yaml
#python coby_args2yaml.py -box 16 16 25 -membrane lipid:POPC:6:params:Amber center:0:0:-7 -protein file:{protein_pdb} charge:-31


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a COBY command-line to a YAML file for AABY"
    )
    parser.add_argument(
        'command',
        nargs=argparse.REMAINDER,
        help="COBY command-line (everything after python3 -m COBY ...)"
    )
    parser.add_argument(
        '-o', '--output',
        default='coby_input.yaml',
        help="Output YAML filename (default: coby_input.yaml)"
    )
    return parser.parse_args()

def main():
    args = parse_args()
    # If user provided a quoted string, split with shlex
    if len(args.command) == 1:
        coby_args = shlex.split(args.command[0])
    else:
        coby_args = args.command

    # Write as YAML list
    with open(args.output, 'w') as f:
        yaml.dump(coby_args, f, default_flow_style=False)
    print(f"COBY YAML written to {args.output} ({len(coby_args)} arguments)")

if __name__ == '__main__':
    main()

