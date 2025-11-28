#!/usr/bin/env python3

import subprocess
import sys
import re
import os
import shutil
import argparse
from pathlib import Path


# Package root: .../site-packages/aaby
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
RESOURCES_DIR = PACKAGE_ROOT / "resources"
MDPS_DIR = RESOURCES_DIR / "mdps"



def gen_rest(path, itp):
    with open(path + itp) as f:
        lines = f.readlines()

    if 'protein' in itp:
        atoms = False
        res_BB = []
        res_SC = []
        for line in lines:
            if atoms and line[0] != ';':
                if line == '\n':
                    break
                l = line.split()
                if l[4] in ['N','CA', 'C', 'O']:
                    # print(line, end='')
                    res_BB.append([l[0], l[4]])
                elif l[4][0] != 'H':
                    # print(line, end='')
                    res_SC.append([l[0], l[4]])

            if '[ atoms ]' in line:
                atoms = True
        restraint = ''
        restraint += '#ifdef POSRES\n'
        restraint +='[ position_restraints ]\n'
        for i in res_BB:
            restraint += f'{i[0]:>5}     1     POSRES_FC_BB    POSRES_FC_BB    POSRES_FC_BB\n'
        for i in res_SC:
            restraint += f'{i[0]:>5}     1     POSRES_FC_SC    POSRES_FC_SC    POSRES_FC_SC\n'
        restraint += '#endif\n'
        with open(path + itp, 'a') as f:
            f.write(restraint)

    
    elif itp.split('.')[0] in ['CHL', 'DLPC', 'DMPC', 'DPPC', 'DSPC', 'DOPC', 'POPC', 'POPE', 'DLPG', 'DMPG', 'DPPG', 'DSPG', 'DOPG', 'POPG', 'DOPS', 'POPS', 'POPA', 'DAPC', 'SDPC', 'PSM', 'SSM']:
        atoms = False
        res = []
        for line in lines:
            if atoms and line[0] != ';':
                if line == '\n':
                    break
                l = line.split()
                if l[4] in ['P31', 'O1']:
                    # print(l)
                    res.append([l[0], l[4]])
            if '[ atoms ]' in line:
                atoms = True
                
        restraint = ''
        restraint += '#ifdef POSRES\n'
        restraint += '[ position_restraints ]\n'
        for i in res:
            restraint += f'{i[0]}     1     0.0             0.0            POSRES_FC_LIPID\n'
        restraint += '#endif\n'
        
        
        with open(path + itp, 'a') as f:
            f.write(restraint)

def main():
    path = 'toppar/'
    itps = os.listdir(path)
    for itp in itps:
        gen_rest(path, itp)
    print(f"[DONE] Restraints added to itp files")

if __name__ == "__main__":
    main()

    