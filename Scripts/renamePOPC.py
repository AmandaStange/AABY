from pathlib import Path

pdb_file = Path("COBY_TER.pdb")
output_file = Path("COBY_TER_lipid21.pdb")

# Define POPC block structure

#ranges_per_lipid = {'POPC': (1, 46, "PA"), (47, 84, "PC"), (85, 134, "OL")], 
#'POPE': [(1, 46, "PA"), (47, 75, "PC"), (76, 125, "OL")]}

ranges_per_lipid = {
"DAPC": [(1, 50, "AR"), (51, 88, "PC"), (89, 138, "AR")],
"DLPC": [(1, 34, "LAL"), (35, 72, "PC"), (73, 106, "LAL")],
"DLPG": [(1, 34, "LAL"), (35, 65, "PGR"), (66, 99, "LAL")],
"DMPC": [(1, 40, "MY"), (41, 78, "PC"), (79, 118, "MY")],
"DMPG": [(1, 40, "MY"), (41, 71, "PGR"), (72, 111, "MY")],
"DOPC": [(1, 50, "OL"), (51, 88, "PC"), (89, 138, "OL")],
"DOPG": [(1, 50, "OL"), (51, 81, "PGR"), (82, 131, "OL")],
"DOPS": [(1, 50, "OL"), (51, 81, "PS"), (82, 131, "OL")],
"DPPC": [(1, 46, "PA"), (47, 84, "PC"), (85, 130, "PA")],
"DPPG": [(1, 46, "PA"), (47, 77, "PGR"), (78, 123, "PA")],
"DSPC": [(1, 52, "ST"), (53, 90, "PC"), (91, 142, "ST")],
"DSPG": [(1, 52, "ST"), (53, 83, "PGR"), (84, 135, "ST")],
"POPA": [(1, 46, "PA"), (47, 66, "PH-"), (67, 116, "OL")],
"POPC": [(1, 46, "PA"), (47, 84, "PC"), (85, 134, "OL")],
"POPE": [(1, 46, "PA"), (47, 75, "PE"), (76, 125, "OL")],
"POPG": [(1, 46, "PA"), (47, 77, "PGR"), (78, 127, "OL")],
"POPS": [(1, 46, "PA"), (47, 77, "PS"), (78, 127, "OL")],
"PSM": [(1, 46, "PA"), (47, 83, "SPM"), (84, 127, "SA")],
"SDPC": [(1, 52, "ST"), (53, 90, "PC"), (91, 142, "DHA")],
"SSM": [(1, 52, "ST"), (53, 89, "SPM"), (90, 133, "SA")]}

#atoms_per_lipid = {'POPC': 134, 'POPE': 125}

atoms_per_lipid = {'DAPC': 138,
'DLPC': 106,
'DLPG': 99,
'DMPC': 118,
'DMPG': 111,
'DOPC': 138,
'DOPG': 131,
'DOPS': 131,
'DPPC': 130,
'DPPG': 123,
'DSPC': 142,
'DSPG': 135,
'POPA': 116,
'POPC': 134,
'POPE': 125,
'POPG': 127,
'POPS': 127,
'PSM': 127,
'SDPC': 142,
'SSM': 133}

lipids = ['DAPC', 'DLPC', 'DLPG', 'DMPC', 'DMPG', 'DOPC', 'DOPG', 'DOPS', 'DPPC', 'DPPG', 'DSPC', 'DSPG', 'POPA', 'POPC', 'POPE', 'POPG', 'POPS', 'PSM', 'SDPC', 'SSM']

def get_resname(index_in_block, lipid):
    for start, end, resname in ranges_per_lipid[lipid]:       
        if start <= index_in_block <= end:
            return resname
    raise ValueError(f"Atom index {index_in_block} not in any range")

def split_popc_and_add_ter(pdb_path, output_path):
    lines_out = []
    popc_atom_counter = 0
    current_resid = None

    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")) and line[17:21].strip() in lipids:
            lipid = line[17:21].strip()
            popc_atom_counter += 1
            index_in_block = (popc_atom_counter - 1) % atoms_per_lipid[lipid] + 1

            # Save original residue number
            resid_str = line[22:26]
            if index_in_block == 1:
                current_resid = resid_str

            new_resname = get_resname(index_in_block, lipid)

            # Apply new resname, keep old resid
            line = line[:17] + new_resname.ljust(4) + line[21:22] + current_resid + line[26:]

            # Insert TER after every 134 atoms (i.e. at end of each POPC block)
            if index_in_block ==  atoms_per_lipid[lipid]:
                lines_out.append(line)
                lines_out.append(f"TER\n")
                popc_atom_counter = 0
                continue

        lines_out.append(line)

    with open(output_path, "w") as f:
        f.writelines(lines_out)

split_popc_and_add_ter(pdb_file, output_file)
print(f"✅ Lipids split into lipid21 naming (residues preserved, TERs inserted) → {output_file}")
