import subprocess
from pathlib import Path

def run(cmd, shell=True, check=True):
    if isinstance(cmd, list):
        print(f"Running: {' '.join(str(x) for x in cmd)}")
    else:
        print(f"Running: {cmd}")
    subprocess.run(cmd, shell=shell, check=check)

def autodetect_amber_files():
    prmtop_candidates = sorted(Path('.').glob('*.prmtop'), key=lambda f: f.stat().st_mtime, reverse=True)
    inpcrd_candidates = sorted(Path('.').glob('*.inpcrd'), key=lambda f: f.stat().st_mtime, reverse=True)
    if not prmtop_candidates or not inpcrd_candidates:
        raise RuntimeError("tleap failed to produce .prmtop or .inpcrd")
    prmtop = prmtop_candidates[0]
    inpcrd = inpcrd_candidates[0]
    print(f"Using prmtop: {prmtop}")
    print(f"Using inpcrd: {inpcrd}")
    return prmtop, inpcrd

def substitute_protein_pdb(args_list, protein_pdb):
    return [str(x).replace("{protein_pdb}", str(protein_pdb)) for x in args_list]

def make_coby_ter_pdb(coby_pdb, ter_pdb, coby_ter_pdb):
    # (Insert logic here, as previously shown)
    pass

def insert_ssbonds_into_tleap(tleap_in_path, ssbond_path):
    # (Insert logic here, as previously shown)
    pass
