from parmed import load_file

# Load AMBER topology + coordinates
amber = load_file("opc_0.15mM.prmtop", "opc_0.15mM.inpcrd")

# Save structure
amber.save("opc_0.15mM.gro", format="gro", overwrite=True)

# Save full GROMACS topology for reference
amber.save("opc_0.15mM.top", format="gromacs", overwrite=True)

# Identify unique residue names
residue_names = sorted(set(res.name for res in amber.residues))
print(f"Molecules found: {residue_names}")

# Save each molecule type to separate .itp
for molname in residue_names:
    atom_indices = [atom.idx for atom in amber.atoms if atom.residue.name == molname]
    mol_sel = amber[atom_indices]  # atom slicing returns a Structure
    filename = f"{molname.lower()}.itp"
    mol_sel.save(filename, format="gromacs", combine=None, overwrite=True)
    print(f"Wrote {filename}")

# Count residue occurrences
from collections import Counter
counts = Counter(res.name for res in amber.residues)

# Write combined topol.top
with open("topol.top", "w") as f:
    f.write("; Combined GROMACS topology\n\n")
    for molname in residue_names:
        f.write(f'#include "{molname.lower()}.itp"\n')
    f.write('\n[ system ]\nOPC box with 0.15 mM NaCl\n\n')
    f.write('[ molecules ]\n')
    for molname in residue_names:
        f.write(f"{molname} {counts[molname]}\n")

