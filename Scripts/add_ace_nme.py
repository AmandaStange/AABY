import sys
from pathlib import Path
import random


def jitter(val, scale=0.20):
    return val + random.uniform(-scale, scale)


def parse_coords(atom_line):
    return (
        float(atom_line[30:38]),
        float(atom_line[38:46]),
        float(atom_line[46:54]),
    )


def parse_resseq(atom_line):
    return atom_line[22:26]


def parse_icode(atom_line):
    return atom_line[26] if atom_line[26].strip() else " "


def format_pdb_atom(atom_num, atom_name, res_name, chain, resseq, icode, x, y, z, element):
    """
    Format a PDB ATOM line.
    Keeps residue number the same as the terminal residue.
    """
    atom_field = f"{atom_name:>4s}"
    return (
        f"ATOM  {atom_num:5d} {atom_field}{res_name:>4s} {chain}"
        f"{resseq}{icode}"
        f"   {x:8.3f}{y:8.3f}{z:8.3f}"
        f"  1.00  0.00          {element:>2s}\n"
    )


def make_ace_atoms(start_atom_num, chain, first_atom_line):
    """
    Create ACE heavy atoms only:
      CH3, C, O

    Same residue number as first residue in chain.
    Coordinates are approximate only.
    """
    x, y, z = parse_coords(first_atom_line)
    resseq = parse_resseq(first_atom_line)
    icode = parse_icode(first_atom_line)

    # Rough placement "before" the N-terminus
    ch3_x, ch3_y, ch3_z = x - 2.4, y, z
    c_x,   c_y,   c_z   = x - 1.2, y + 0.1, z
    o_x,   o_y,   o_z   = x - 0.2, y + 0.3, z

    atoms = [
        format_pdb_atom(start_atom_num + 0, "CH3", "ACE", chain, resseq, icode,
                        jitter(ch3_x), jitter(ch3_y), jitter(ch3_z), "C"),
        format_pdb_atom(start_atom_num + 1, "C",   "ACE", chain, resseq, icode,
                        jitter(c_x),   jitter(c_y),   jitter(c_z),   "C"),
        format_pdb_atom(start_atom_num + 2, "O",   "ACE", chain, resseq, icode,
                        jitter(o_x),   jitter(o_y),   jitter(o_z),   "O"),
    ]
    return atoms


def make_nme_atoms(start_atom_num, chain, last_atom_line):
    """
    Create NME heavy atoms only:
      N, C

    Same residue number as last residue in chain.
    Coordinates are approximate only.
    """
    x, y, z = parse_coords(last_atom_line)
    resseq = parse_resseq(last_atom_line)
    icode = parse_icode(last_atom_line)

    # Rough placement "after" the C-terminus
    n_x, n_y, n_z = x + 1.3, y, z
    c_x, c_y, c_z = x + 2.6, y + 0.2, z

    atoms = [
        format_pdb_atom(start_atom_num + 0, "N", "NME", chain, resseq, icode,
                        jitter(n_x), jitter(n_y), jitter(n_z), "N"),
        format_pdb_atom(start_atom_num + 1, "C", "NME", chain, resseq, icode,
                        jitter(c_x), jitter(c_y), jitter(c_z), "C"),
    ]
    return atoms


def main():
    if len(sys.argv) != 3:
        print("Usage: python add_ace_nme.py <pdb_file> <chainids_comma_separated>")
        sys.exit(1)

    pdb_file = Path(sys.argv[1])
    chains = set(c.strip() for c in sys.argv[2].split(",") if c.strip())

    with open(pdb_file) as f:
        lines = f.readlines()

    # Find highest atom number in the PDB
    max_atom_num = 0
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            try:
                max_atom_num = max(max_atom_num, int(l[6:11]))
            except ValueError:
                pass

    # Collect atom line indices for requested chains
    chain_indices = {c: [] for c in chains}
    for i, l in enumerate(lines):
        if l.startswith(("ATOM", "HETATM")):
            chain_id = l[21]
            if chain_id in chains:
                chain_indices[chain_id].append(i)

    insertions = []

    for chain, indices in chain_indices.items():
        if not indices:
            print(f"Warning: chain {chain} not found.")
            continue

        first_idx = indices[0]
        last_idx = indices[-1]

        ace_atoms = make_ace_atoms(max_atom_num + 1, chain, lines[first_idx])
        max_atom_num += len(ace_atoms)

        nme_atoms = make_nme_atoms(max_atom_num + 1, chain, lines[last_idx])
        max_atom_num += len(nme_atoms)

        # ACE before first atom of chain
        for offset, atom_line in enumerate(ace_atoms):
            insertions.append((first_idx + offset, atom_line))

        # NME after last atom of chain
        for offset, atom_line in enumerate(nme_atoms):
            insertions.append((last_idx + 1 + offset, atom_line))

    # Insert from back to front so indices remain valid
    for idx, new_line in sorted(insertions, key=lambda x: x[0], reverse=True):
        lines.insert(idx, new_line)

    out_file = pdb_file.with_name(pdb_file.stem + "_ACE_NME.pdb")
    with open(out_file, "w") as f:
        f.writelines(lines)

    print(f"Wrote {out_file}")


if __name__ == "__main__":
    main()