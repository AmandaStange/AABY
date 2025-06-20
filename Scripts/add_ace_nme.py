import sys
from pathlib import Path
import random

def make_kicked_coords(x, y, z, direction=1):
    """Add a random small displacement in each direction, direction = +1 or -1 for NME/ACE."""
    # direction: +1 for NME, -1 for ACE
    dx = direction * (1.2 + random.uniform(0.1, 0.5))
    dy = random.uniform(-0.4, 0.4)
    dz = random.uniform(-0.4, 0.4)
    return x + dx, y + dy, z + dz

def make_ace_atom(atom_num, chain, first_atom_line):
    # Extract coordinates from the first atom
    x = float(first_atom_line[30:38])
    y = float(first_atom_line[38:46])
    z = float(first_atom_line[46:54])
    xn, yn, zn = make_kicked_coords(x, y, z, direction=-1)
    return (
        f"ATOM  {atom_num:5d}  C   ACE {chain}{first_atom_line[22:26]}{first_atom_line[26:30]}"
        f"{xn:8.3f}{yn:8.3f}{zn:8.3f}  1.00  0.00\n"
    )

def make_nme_atom(atom_num, chain, last_atom_line):
    # Extract coordinates from the last atom
    x = float(last_atom_line[30:38])
    y = float(last_atom_line[38:46])
    z = float(last_atom_line[46:54])
    xn, yn, zn = make_kicked_coords(x, y, z, direction=1)
    return (
        f"ATOM  {atom_num:5d}  C   NME {chain}{last_atom_line[22:26]}{last_atom_line[26:30]}"
        f"{xn:8.3f}{yn:8.3f}{zn:8.3f}  1.00  0.00\n"
    )

def main():
    if len(sys.argv) != 3:
        print("Usage: python add_ace_nme.py <pdb_file> <chainids_comma_separated>")
        sys.exit(1)

    pdb_file = Path(sys.argv[1])
    chains = set(sys.argv[2].split(","))

    with open(pdb_file) as f:
        lines = f.readlines()

    # Find highest atom number in the PDB to continue numbering
    max_atom_num = 0
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            try:
                max_atom_num = max(max_atom_num, int(l[6:11]))
            except Exception:
                pass

    # Map: chain -> [atom_lines_indices]
    chain_indices = {c: [] for c in chains}
    for i, l in enumerate(lines):
        if l.startswith(("ATOM", "HETATM")):
            c = l[21]
            if c in chains:
                chain_indices[c].append(i)

    # Prepare insertions
    insertions = []
    for chain, indices in chain_indices.items():
        if not indices:
            print(f"Warning: chain {chain} not found.")
            continue
        first_idx = indices[0]
        last_idx = indices[-1]
        max_atom_num += 1
        ace_atom = make_ace_atom(max_atom_num, chain, lines[first_idx])
        max_atom_num += 1
        nme_atom = make_nme_atom(max_atom_num, chain, lines[last_idx])

        insertions.append((first_idx, ace_atom))       # ACE goes *before* first atom
        insertions.append((last_idx + 1, nme_atom))    # NME goes *after* last atom

    # Insert from back so indices don't shift
    for idx, new_line in sorted(insertions, reverse=True):
        lines.insert(idx, new_line)

    # Write output
    out_file = pdb_file.with_name(pdb_file.stem + "_ACE_NME.pdb")
    with open(out_file, "w") as f:
        f.writelines(lines)
    print(f"Wrote {out_file}")

if __name__ == "__main__":
    main()
