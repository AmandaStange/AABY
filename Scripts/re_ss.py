import re
import sys
from pathlib import Path

def extract_residue_mapping(old_pdb, new_pdb):
    """Map (resname, chain, resid) in old to new resid."""
    def get_residues(pdb_path):
        residues = []
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    resname = line[17:20].strip()
                    chain = line[21]
                    resid = int(line[22:26])
                    atom = line[12:16].strip()
                    key = (resname, chain, atom)
                    residues.append((key, resid))
        return residues

    old = get_residues(old_pdb)
    new = get_residues(new_pdb)

    mapping = {}
    for (okey, old_resid), (_, new_resid) in zip(old, new):
        resname, chain, atom = okey
        mapping[(resname, chain, old_resid)] = new_resid
    return mapping

def renumber_ssbond(ssbond_txt, mapping, output_txt):
    with open(ssbond_txt) as f:
        lines = f.readlines()

    updated_lines = []
    for line in lines:
        match = re.match(r"(SSBOND\s+\d+\s+CYS\s+)(\w)\s+(\d+)(\s+CYS\s+)(\w)\s+(\d+)", line)
        if match:
            chain1, res1 = match.group(2), int(match.group(3))
            chain2, res2 = match.group(5), int(match.group(6))
            new_res1 = mapping.get(("CYS", chain1, res1), res1)
            new_res2 = mapping.get(("CYS", chain2, res2), res2)
            new_line = (f"{match.group(1)}{chain1:>2}{new_res1:>4}{match.group(4)}"
                        f"{chain2:>2}{new_res2:>4}{line[match.end():]}")
            updated_lines.append(new_line)
        else:
            updated_lines.append(line)

    with open(output_txt, "w") as f:
        f.writelines(updated_lines)

def parse_renumbered_ssbond(ssbond_renumbered_txt):
    # returns set of (chain, resid) for all CYX (formerly CYS) in ssbonds
    ss_cyx = set()
    with open(ssbond_renumbered_txt) as f:
        for line in f:
            if line.startswith('SSBOND'):
                chain1 = line[16].strip()
                res1 = int(line[17:21].strip())
                chain2 = line[30].strip()
                res2 = int(line[31:35].strip())
                ss_cyx.add((chain1, res1))
                ss_cyx.add((chain2, res2))
    return ss_cyx

def write_cyx_pdb(pdbfile, cyx_set, outname):
    with open(pdbfile) as f:
        lines = f.readlines()
    with open(outname, 'w') as f:
        for line in lines:
            if (line.startswith("ATOM") or line.startswith("HETATM")):
                resname = line[17:20]
                chain = line[21]
                resid = int(line[22:26])
                # print(cyx_set)
                if resname.strip() == "CYS" and (chain, resid) in cyx_set:
                    line = line[:17] + "CYX" + line[20:]
            f.write(line)

# Usage
if len(sys.argv) < 2:
    print("Usage: python re_ss.py <pdb_file>")
    sys.exit(1)

pdb_file = sys.argv[1]
pdb_file_re = pdb_file.replace(".pdb", "_re.pdb")

mapping = extract_residue_mapping(
    pdb_file,
    pdb_file_re
)

renumber_ssbond("SSBOND.txt", mapping, "SSBOND_renumbered.txt")

with open('SSBOND_renumbered.txt') as f:
    lines = f.readlines()

with open('tleap_SSBONDs.txt', 'w') as f:
    for line in lines:
        if not line.startswith("SSBOND"):
            continue
        res1 = line[17:21].strip()
        res2 = line[31:35].strip()
        f.write(f'bond ramp.{res1}.SG ramp.{res2}.SG\n')

def insert_ssbonds_into_tleap(tleap_in_path='Scripts/tleap.in', ssbond_path='tleap_SSBONDs.txt'):
    tleap_in = Path(tleap_in_path)
    ssbonds = Path(ssbond_path)
    if not tleap_in.exists() or not ssbonds.exists():
        print(f"Warning: {tleap_in} or {ssbonds} does not exist. Skipping SSBOND insertion.")
        return

    with tleap_in.open() as f:
        lines = f.readlines()

    # Find the line with loadpdb (should be only one)
    for idx, line in enumerate(lines):
        if 'loadpdb' in line:
            loadpdb_idx = idx
            break
    else:
        print("Could not find 'loadpdb' in tleap.in!")
        return

    with ssbonds.open() as f:
        ssbond_lines = f.readlines()

    # Insert after loadpdb line
    new_lines = lines[:loadpdb_idx + 1] + ssbond_lines + lines[loadpdb_idx + 1:]

    # Overwrite tleap.in
    with tleap_in.open('w') as f:
        f.writelines(new_lines)

    print(f"Inserted {ssbonds} after loadpdb in {tleap_in}")

# --- New Part: Write *_re_ssbond.pdb with CYX ---

cyx_set = parse_renumbered_ssbond("SSBOND_renumbered.txt")
pdb_file_re_ssbond = Path(pdb_file).with_name(Path(pdb_file_re).stem + "_ssbond.pdb")
write_cyx_pdb(pdb_file_re, cyx_set, pdb_file_re_ssbond)
print(f"Wrote {pdb_file_re_ssbond}")
