import sys
import os
import re
from parmed import load_file
import shutil

def parmed_convert(prmtop_file, inpcrd_file, top_file, gro_file):
    amber = load_file(prmtop_file, inpcrd_file)
    amber.save(gro_file, format="gro", overwrite=True)
    amber.save(top_file, format="gromacs", overwrite=True)
    print(f"Wrote {gro_file} and {top_file}")

def split_top_to_itps(top_file, toppar_dir="toppar", output_top='topol.top'):
    with open(top_file, 'r') as f:
        lines = f.readlines()

    # Gather includes and preamble up to first [ moleculetype ]
    preamble = []
    i = 0
    while i < len(lines):
        if re.match(r'^\s*\[\s*moleculetype\s*\]', lines[i], re.IGNORECASE):
            break
        preamble.append(lines[i])
        i += 1

    # Extract each moleculetype block into its own .itp
    itp_blocks = []
    itp_names = []
    while i < len(lines):
        if re.match(r'^\s*\[\s*moleculetype\s*\]', lines[i], re.IGNORECASE):
            block = []
            while i < len(lines):
                if (re.match(r'^\s*\[\s*moleculetype\s*\]', lines[i], re.IGNORECASE) and block) or \
                   (re.match(r'^\s*\[\s*system\s*\]', lines[i], re.IGNORECASE)):
                    break
                block.append(lines[i])
                i += 1
            # Get moleculetype name for filename and [ molecules ] block
            for j, line in enumerate(block):
                if j == 0:
                    continue
                if not line.strip().startswith(';') and line.strip():
                    mtname = line.split()[0]
                    itp_names.append(mtname)
                    itp_blocks.append(block)
                    break
        else:
            i += 1

    # Find [ molecules ] section to extract molecule counts
    mol_section = []
    found_mol = False
    for idx, line in enumerate(lines):
        if re.match(r'^\s*\[\s*molecules\s*\]', line, re.IGNORECASE):
            found_mol = True
            continue
        if found_mol:
            if line.strip().startswith('['):
                break
            if line.strip() and not line.strip().startswith(';'):
                mol_section.append(line)

    # Write ITP files to toppar directory
    os.makedirs(toppar_dir, exist_ok=True)
    for name, block in zip(itp_names, itp_blocks):
        fname = os.path.join(toppar_dir, f"{name}.itp")
        with open(fname, "w") as fout:
            fout.writelines(block)
        print(f"Wrote {fname}")

    # Write new topol.top with #include "toppar/NAME.itp"
    with open(output_top, "w") as f:
        f.writelines(preamble)
        for name in itp_names:
            f.write(f'#include "{toppar_dir}/{name}.itp"\n')
        f.write("\n[ system ]\n; Name\nGeneric title\n\n")
        f.write("[ molecules ]\n; Compound       #mols\n")
        for line in mol_section:
            f.write(line)
    print(f"Wrote {output_top}")

def split_ff_from_topol(topol_file, ff_file="ff.itp", new_topol="topol.top"):
    with open(topol_file, "r") as f:
        lines = f.readlines()

    ff_block = []
    rest_block = []
    found_first_include = False

    for idx, line in enumerate(lines):
        if not found_first_include and line.strip().startswith("#include"):
            found_first_include = True
            rest_block = lines[idx:]  # include all from this include onward
            break
        if not found_first_include:
            ff_block.append(line)

    # Write ff.itp
    with open(ff_file, "w") as fff:
        fff.writelines(ff_block)
    print(f"Wrote {ff_file}")

    # Write new topol.top (include ff.itp at top, then the rest)
    with open(new_topol, "w") as ftop:
        ftop.write(f'#include "{ff_file}"\n')
        for line in rest_block:
            ftop.write(line)
    print(f"Wrote {new_topol}")

def main(prmtop_file, inpcrd_file):
    top_file = "system.top"
    gro_file = "system.gro"
    parmed_convert(prmtop_file, inpcrd_file, top_file, gro_file)
    split_top_to_itps(top_file, toppar_dir="toppar", output_top='topol.top.temp')
    split_ff_from_topol('topol.top.temp', ff_file="ff.itp", new_topol="topol.top")
    os.remove("topol.top.temp")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_and_split_full.py file.prmtop file.inpcrd")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
