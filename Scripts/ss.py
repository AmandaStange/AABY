def extract_cys_resids(pdb_path):
    """Extract CYX residue (chain, resid) in order from a PDB file."""
    resids = []
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[17:20] == "CYX":
                resid = int(line[22:26])
                chain = line[21]
                key = (chain, resid)
                if key not in seen:
                    seen.add(key)
                    resids.append(key)
    return resids

def make_resid_mapping(old_pdb, new_pdb):
    """Map (chain, old_resid) -> new_resid based on CYX order."""
    old = extract_cys_resids(old_pdb)
    new = extract_cys_resids(new_pdb)
    if len(old) != len(new):
        raise ValueError(f"CYX count mismatch: {len(old)} vs {len(new)}")
    return {old[i]: new[i][1] for i in range(len(old))}

def renumber_ssbond_file(ssbond_txt, mapping, output_txt):
    """Update SSBOND.txt with new residue numbers."""
    import re
    with open(ssbond_txt) as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        match = re.match(r"(SSBOND\s+\d+\s+CYS\s+)(\w)\s+(\d+)(\s+CYS\s+)(\w)\s+(\d+)", line)
        if match:
            chain1, res1 = match.group(2), int(match.group(3))
            chain2, res2 = match.group(5), int(match.group(6))
            new_res1 = mapping.get((chain1, res1), res1)
            new_res2 = mapping.get((chain2, res2), res2)
            new_line = (
                f"{match.group(1)}{chain1:>2}{new_res1:>4}"
                f"{match.group(4)}{chain2:>2}{new_res2:>4}{line[match.end():]}"
            )
            new_lines.append(new_line)
        else:
            new_lines.append(line)

    with open(output_txt, "w") as f:
        f.writelines(new_lines)

# === Usage ===
old_pdb = "7SL6_medoid_cluster_1_AA_renumbered_SSBOND_fix_heavy.pdb"
new_pdb = "7SL6_medoid_cluster_1_AA_renumbered_SSBOND_fix_heavy_re.pdb"
ssbond_txt = "SSBOND.txt"
output_txt = "SSBOND_renumbered.txt"

mapping = make_resid_mapping(old_pdb, new_pdb)
renumber_ssbond_file(ssbond_txt, mapping, output_txt)

