
import argparse
from pathlib import Path

LIPID_RESNAMES = {"CHL", "PA", "PC", "OL"}
KNOWN_BLOCKS = {
    "CHL": {
        "POSRES": [
            "#ifdef POSRES",
            "[ position_restraints ]",
            "   73     1     0.0             0.0            POSRES_FC_LIPID",
            "#endif"
        ],
        "DIHRES": [
            "#ifdef DIHRES",
            "[ dihedral_restraints ]",
            "#endif"
        ]
    },
    "POPC": {
        "POSRES": [
            "#ifdef POSRES",
            "[ position_restraints ]",
            "   59     1     0.0             0.0            POSRES_FC_LIPID",
            "#endif"
        ],
        "DIHRES": [
            "#ifdef DIHRES",
            "[ dihedral_restraints ]",
            "   50    55    53    82     1    120.0      2.5       DIHRES_FC",
            "  103   106   108   110     1      0.0      0.0       DIHRES_FC",
            "#endif"
        ]
    }
}

def extract_resname_and_atoms(itp_file):
    molname = None
    resnames = set()
    atoms = []
    inside_atoms = False

    with open(itp_file) as f:
        lines = iter(f)
        for line in lines:
            if line.strip().startswith("[ moleculetype ]"):
                while True:
                    mol_line = next(lines).strip()
                    if mol_line and not mol_line.startswith(";"):
                        molname = mol_line.split()[0]
                        break
            elif line.strip().startswith("[ atoms ]"):
                inside_atoms = True
                continue
            elif inside_atoms and line.strip().startswith("["):
                break
            elif inside_atoms and line.strip() and not line.strip().startswith(";"):
                parts = line.split()
                if len(parts) >= 5:
                    atom_idx = int(parts[0])
                    atom_name = parts[4]
                    resname = parts[3]
                    resnames.add(resname)
                    atoms.append((atom_idx, atom_name))

    is_lipid = any(r in LIPID_RESNAMES for r in resnames)
    return molname, is_lipid, atoms, resnames

def insert_inline_posres_block(itp_file, atom_lines, is_lipid=False):
    lines_out = ["#ifdef POSRES", "[ position_restraints ]", "; atom  type  fx                  fy                  fz"]
    for atom_idx, atom_name in atom_lines:
        if atom_name.startswith("H"):  # Skip hydrogens
            continue
        if is_lipid:
            fc = "POSRES_FC_LIPID"
        else:
            fc = "POSRES_FC_BB" if atom_name in ("N", "CA", "C") else "POSRES_FC_SC"
        lines_out.append(f"{atom_idx:6d}     1   {fc:<20} {fc:<20} {fc:<20}")
    lines_out.append("#endif")

    with open(itp_file, "r") as f:
        lines = f.readlines()

    insert_at = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("[ moleculetype ]"):
            insert_at = i
            while insert_at < len(lines) and not lines[insert_at].strip().startswith("[", 1):
                insert_at += 1
            break

    lines[insert_at:insert_at] = [l + "\n" for l in lines_out] + ["\n"]

    with open(itp_file, "w") as f:
        f.writelines(lines)

def insert_inline_dihres_block(itp_file, dih_entries=None):
    lines_out = ["#ifdef DIHRES", "[ dihedral_restraints ]", "; ai  aj  ak  al  type  phi  dphi  kfac"]
    if dih_entries is None:
        dih_entries = [(1, 2, 3, 4, 180.0, 30.0)]
    for a1, a2, a3, a4, phi, dphi in dih_entries:
        lines_out.append(f"{a1:5d} {a2:5d} {a3:5d} {a4:5d}     1  {phi:6.1f}  {dphi:5.1f}  DIHRES_FC")
    lines_out.append("#endif")

    with open(itp_file, "r") as f:
        lines = f.readlines()

    insert_at = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("[ moleculetype ]"):
            insert_at = i
            while insert_at < len(lines) and not lines[insert_at].strip().startswith("[", 1):
                insert_at += 1
            break

    lines[insert_at:insert_at] = [l + "\n" for l in lines_out] + ["\n"]

    with open(itp_file, "w") as f:
        f.writelines(lines)

def insert_known_blocks(itp_file, posres_block, dihres_block):
    with open(itp_file, "r") as f:
        lines = f.readlines()

    insert_at = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("[ moleculetype ]"):
            insert_at = i
            while insert_at < len(lines) and not lines[insert_at].strip().startswith("[", 1):
                insert_at += 1
            break

    block_lines = [l + "\n" for l in posres_block + [""] + dihres_block + [""]]
    lines[insert_at:insert_at] = block_lines

    with open(itp_file, "w") as f:
        f.writelines(lines)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--itp", required=True, help="Path to molecule .itp file")
    args = parser.parse_args()

    molname, is_lipid, atoms, resnames = extract_resname_and_atoms(args.itp)

    if "CHL" in resnames:
        insert_known_blocks(args.itp, KNOWN_BLOCKS["CHL"]["POSRES"], KNOWN_BLOCKS["CHL"]["DIHRES"])
        print(f"✅ Inserted CHL POSRES/DIHRES blocks into {args.itp}")
    elif any(r in {"PA", "PC", "OL"} for r in resnames):
        insert_known_blocks(args.itp, KNOWN_BLOCKS["POPC"]["POSRES"], KNOWN_BLOCKS["POPC"]["DIHRES"])
        print(f"✅ Inserted POPC POSRES/DIHRES blocks into {args.itp}")
    else:
        insert_inline_posres_block(args.itp, atoms, is_lipid=is_lipid)
        print(f"✅ Inserted inline POSRES into {args.itp}")

if __name__ == "__main__":
    main()
