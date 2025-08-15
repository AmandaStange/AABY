from pathlib import Path
import string

input_file = Path("COBY_TER_lipid21.pdb")
output_file = input_file.with_name(input_file.stem + "_chains.pdb")

protein_residues = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "HIE", "HID", "HIP", "CYX",
    "LYN", "GLH", "ASH", "ACE", "NME","Na", "Cl","K", "Ca", "Mg", "Zn","NA", "CL","K", "CA", "MG", "ZN"
}

lipid_residues = ['SA', 'PGR', 'AR', 'LAL', 'ST', 'DHA', 'PE', 'MY', 'OL', 'PH-', 'PS', 'PC', 'SPM', 'PA','CHL']


chain_labels = list(string.ascii_uppercase)
current_chain = 0
in_protein_block = False

with open(input_file) as inp, open(output_file, "w") as out:
    for line in inp:
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if resname not in lipid_residues:
                line = line[:21] + chain_labels[current_chain] + line[22:]
                in_protein_block = True
            else:
                line = line[:21] + "X" + line[22:]
                in_protein_block = False
            out.write(line)
        elif line.startswith("TER"):
            # Only modify if TER line is long enough (has chain info)
            if len(line.rstrip()) == 3:  # 'TER\n' or 'TER\r\n'
                out.write(line)
            elif len(line) >= 22:
                if in_protein_block:
                    line = line[:21] + chain_labels[current_chain] + line[22:]
                    current_chain += 1
                else:
                    line = line[:21] + "X" + line[22:]
                in_protein_block = False
                out.write(line)
            else:
                out.write(line)
        else:
            out.write(line)

print(f"Wrote: {output_file}")
