from pathlib import Path

pdb_file = Path("COBY_TER.pdb")
output_file = Path("COBY_TER_lipid21.pdb")

# Define POPC block structure
ranges = [(1, 46, "PA"), (47, 84, "PC"), (85, 134, "OL")]
atoms_per_popc = 134

def get_resname(index_in_block):
    for start, end, resname in ranges:
        if start <= index_in_block <= end:
            return resname
    raise ValueError(f"Atom index {index_in_block} not in any range")

def split_popc_and_add_ter(pdb_path, output_path):
    lines_out = []
    popc_atom_counter = 0
    current_resid = None

    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")) and line[17:21].strip() == "POPC":
            popc_atom_counter += 1
            index_in_block = (popc_atom_counter - 1) % atoms_per_popc + 1

            # Save original residue number
            resid_str = line[22:26]
            if index_in_block == 1:
                current_resid = resid_str

            new_resname = get_resname(index_in_block)

            # Apply new resname, keep old resid
            line = line[:17] + new_resname.ljust(4) + line[21:22] + current_resid + line[26:]

            # Insert TER after every 134 atoms (i.e. at end of each POPC block)
            if index_in_block == atoms_per_popc:
                lines_out.append(line)
                lines_out.append(f"TER\n")
                continue

        lines_out.append(line)

    with open(output_path, "w") as f:
        f.writelines(lines_out)

split_popc_and_add_ter(pdb_file, output_file)
print(f"✅ POPC split into PA/PC/OL (residues preserved, TERs inserted) → {output_file}")
