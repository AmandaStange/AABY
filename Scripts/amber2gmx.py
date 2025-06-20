import sys
from parmed import load_file
from collections import Counter
import os

# if len(sys.argv) != 2:
#     print("Usage: python amber2gmx_with_posres.py <prefix>")
#     sys.exit(1)

# top_name = sys.argv[1]
top_name = '7SL6'
prmtop = f"{top_name}.prmtop"
inpcrd = f"{top_name}.inpcrd"

# === LOAD STRUCTURE ===
amber = load_file(prmtop, inpcrd)

# === SAVE GRO + TOP ===
amber.save(f"{top_name}.gro", format="gro", overwrite=True)
amber.save(f"{top_name}.top", format="gromacs", overwrite=True)



# === GET RESIDUE TYPES ===
residue_names = sorted(set(res.name for res in amber.residues))
print(f"Molecules found: {residue_names}")

# === WRITE ITP FILES ===
for molname in residue_names:
    atom_indices = [atom.idx for atom in amber.atoms if atom.residue.name == molname]
    mol_sel = amber[atom_indices]
    filename = f"{molname.lower()}.itp"
    mol_sel.save(filename, format="gromacs", combine=None, overwrite=True)
    print(f"Wrote {filename}")

    # Write position restraints for heavy atoms (only first molecule)
    posre_lines = ["[ position_restraints ]\n; atom  type  fx  fy  fz\n"]
    for atom in mol_sel.atoms:
        if atom.atomic_number > 1:
            posre_lines.append(f"{atom.idx + 1:6d}    1   1000 1000 1000\n")
    with open("posre.itp", "w") as f:
        f.writelines(posre_lines)
    print(f"Wrote posre.itp for {molname}")
    break  # One posre is enough

# === COUNT RESIDUES ===
counts = Counter(res.name for res in amber.residues)

# === WRITE topol.top ===
with open("topol.top", "w") as f:
    f.write("; GROMACS topology\n\n")
    for molname in residue_names:
        f.write(f'#include "{molname.lower()}.itp"\n')
    f.write(f'#include "posre.itp"\n')
    f.write('\n[ system ]\n; Name\n' + top_name + ' system\n\n')
    f.write('[ molecules ]\n; Compound  #mols\n')
    for molname in residue_names:
        f.write(f"{molname} {counts[molname]}\n")

print("✅ Done. Files written:")
print(f"- {top_name}.gro")
print(f"- {top_name}.top")
print(f"- topol.top")
print(f"- posre.itp")
print(f"- {[f'{r.lower()}.itp' for r in residue_names]}")
