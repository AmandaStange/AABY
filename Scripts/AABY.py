#!/usr/bin/env python3
import time
start_time = time.time()

import argparse
import subprocess
from pathlib import Path
import sys
import os
import shutil
import yaml

def run(cmd, shell=True, check=True):
    if isinstance(cmd, list):
        print(f"Running: {' '.join(str(x) for x in cmd)}")
    else:
        print(f"Running: {cmd}")
    subprocess.run(cmd, shell=shell, check=check)

def insert_ssbonds_into_tleap(tleap_in_path='tleap.in', ssbond_path='tleap_SSBONDs.txt'):
    tleap_in = Path(tleap_in_path)
    ssbonds = Path(ssbond_path)
    if not tleap_in.exists() or not ssbonds.exists():
        print(f"Warning: {tleap_in} or {ssbonds} does not exist. Skipping SSBOND insertion.")
        return
    with tleap_in.open() as f:
        lines = f.readlines()
    for idx, line in enumerate(lines):
        if 'loadpdb' in line:
            loadpdb_idx = idx
            break
    else:
        print("Could not find 'loadpdb' in tleap.in!")
        return
    with ssbonds.open() as f:
        ssbond_lines = f.readlines()
    new_lines = lines[:loadpdb_idx + 1] + ssbond_lines + lines[loadpdb_idx + 1:]
    with tleap_in.open('w') as f:
        f.writelines(new_lines)
    print(f"Inserted {ssbonds} after loadpdb in {tleap_in}")

def make_coby_ter_pdb(coby_pdb, pre_coby_pdb, coby_ter_pdb):
    """
    Insert TER records into the COBY PDB wherever a TER occurred in the pre-COBY PDB,
    by matching residue numbers. Assumes all chains become A in COBY output.
    """
    # 1. Find residue numbers before each TER in the pre-COBY file
    ter_resids = []
    prev_resid = None
    with open(pre_coby_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                prev_resid = int(line[22:26])
            elif line.startswith("TER"):
                if prev_resid is not None:
                    ter_resids.append(prev_resid)
                prev_resid = None

    # 2. Parse COBY PDB, insert TER after each matching resid
    output_lines = []
    last_resid = None
    last_atomnum = None
    last_resname = None

    for i, line in enumerate(open(coby_pdb)):
        output_lines.append(line)
        if line.startswith(("ATOM", "HETATM")):
            resid = int(line[22:26])
            atomnum = int(line[6:11])
            resname = line[17:20]
            if last_resid is not None and last_resid in ter_resids and resid != last_resid:
                # Insert TER *after* last atom of the pre-TER residue
                ter_line = f"TER   {last_atomnum:5d}      {last_resname} A{last_resid:>4}\n"
                output_lines.insert(-1, ter_line)  # Insert before the current line
            last_resid = resid
            last_atomnum = atomnum
            last_resname = resname

    # Also handle case where the last residue should get a TER (optional)
    # (Not usually needed; leave as is)

    with open(coby_ter_pdb, 'w') as f:
        f.writelines(output_lines)
    print(f"Wrote merged COBY/protein file with TERs: {coby_ter_pdb}")




def substitute_protein_pdb(args_list, protein_pdb):
    return [str(x).replace("{protein_pdb}", str(protein_pdb)) for x in args_list]

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

def main():
    parser = argparse.ArgumentParser(description="AABY: End-to-end AMBER system builder from PDB")
    parser.add_argument('-f', '--input', required=True, help='Input PDB file')
    parser.add_argument('--chains', default="A,B,C,D", help='Comma-separated list of protein chains for ACE/NME')
    parser.add_argument('-r', '--replicas', type=int, default=1, help='Number of resolvated replicas')
    parser.add_argument('--ssbond', action='store_true', help='Renumber SSBOND.txt, convert CYS to CYX, insert ssbonds into tleap.in')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--ph', type=float, help='pH for propka (protonation by predicted pKa)')
    group.add_argument('--protlist', type=str, help='Residue list file for direct protonation changes')
    parser.add_argument('--coby-yaml', type=str, help='YAML file specifying COBY arguments (use {protein_pdb} as placeholder)')
    parser.add_argument('--coby-args', nargs=argparse.REMAINDER, help='Arguments passed directly to COBY after this flag (use {protein_pdb} as placeholder)')
    args = parser.parse_args()


    # Copy tleap.in to local working directory
    tleap_template = Path('Scripts/tleap.in')
    tleap_local = Path('tleap.in')
    shutil.copy(tleap_template, tleap_local)

    pdb = Path(args.input)
    base = pdb.stem

    # 1. Add caps
    out_ace_nme = pdb.with_name(base + '_ACE_NME.pdb')
    run(f'python Scripts/add_ace_nme.py {pdb} {args.chains}')
    assert out_ace_nme.exists(), "add_ace_nme failed"

    # 2. Renumber
    out_renum = pdb.with_name(base + '_ACE_NME_re.pdb')
    run(f'gmx editconf -f {out_ace_nme} -resnr 1 -o {out_renum}')
    assert out_renum.exists(), "editconf renumber failed"

    # 3. SSBOND handling if requested
    if args.ssbond:
        run(f'python Scripts/re_ss.py {out_ace_nme}')
        ssbonded_pdb = out_ace_nme.with_name(out_renum.stem + '_ssbond.pdb')
        assert ssbonded_pdb.exists(), "re_ss failed"
        working_pdb = ssbonded_pdb
        run(f'cp Scripts/tleap.in .')
        insert_ssbonds_into_tleap('tleap.in', 'tleap_SSBONDs.txt')
    else:
        working_pdb = out_renum

    # 4. Delete H (heavy atom only)
    out_heavy = working_pdb.with_name(working_pdb.stem + '_heavy.pdb')
    run(f'echo 2 | gmx trjconv -f {working_pdb} -s {working_pdb} -o {out_heavy}')
    assert out_heavy.exists(), "trjconv failed"

    # 5. Protonation (by pH/pKa or residue list)
    if args.ph is not None:
        pka_file = out_heavy.with_suffix('.pka')
        run(f'python Scripts/pka.py {out_heavy} {args.ph}')
        assert pka_file.exists(), "pka calculation failed"
        run(f'python Scripts/change_prot.py {pka_file} --ph {args.ph} -o mutate_residues.sh')
        shutil.copy(out_heavy, out_heavy.with_name(out_heavy.stem + '_renamed.pdb'))
        renamed_pdb = out_heavy.with_name(out_heavy.stem + '_renamed.pdb')
        run(f'bash mutate_residues.sh {renamed_pdb}')
    elif args.protlist:
        run(f'python Scripts/change_prot.py {args.protlist} --list -o mutate_residues.sh')
        shutil.copy(out_heavy, out_heavy.with_name(out_heavy.stem + '_renamed.pdb'))
        renamed_pdb = out_heavy.with_name(out_heavy.stem + '_renamed.pdb')
        run(f'bash mutate_residues.sh {renamed_pdb}')
    else:
        renamed_pdb = out_heavy.with_name(out_heavy.stem + '_renamed.pdb')
        shutil.copy(out_heavy, renamed_pdb)

    # 6. Rename CD ILE and HIS
    run(f'sed -i "s/CD  ILE/CD1 ILE/g" {renamed_pdb}')
    run(f'sed -i "s/HIS/HIE/g" {renamed_pdb}')

    # 7. Add TER
    ter_pdb = renamed_pdb.with_name(renamed_pdb.stem + '_TER.pdb')
    run(f'python Scripts/add_ter.py {renamed_pdb} {ter_pdb}')
    assert ter_pdb.exists(), "add_ter failed"

    # Always remove END lines
    with open(str(ter_pdb), 'r') as f_in:
        ter_cleaned = [line for line in f_in if not line.startswith('END')]
    with open(str(ter_pdb), 'w') as f_out:
        f_out.writelines(ter_cleaned)

    # Prepare for COBY or direct to TLEAP
    last_pdb_before_coby = ter_pdb
    coby_input_pdb = Path("COBY_input.pdb")
    shutil.copy(last_pdb_before_coby, coby_input_pdb)

    # 8. Run COBY if and only if coby-yaml or coby-args are given
    coby_should_run = bool(args.coby_yaml or args.coby_args)
    if coby_should_run:
        if args.coby_yaml:
            with open(args.coby_yaml) as f:
                coby_args = yaml.safe_load(f)
            if not isinstance(coby_args, list):
                print(f"YAML file must be a list of arguments, got {type(coby_args)}")
                sys.exit(1)
            coby_args = substitute_protein_pdb(coby_args, coby_input_pdb)
            coby_cmd = ["python3", "-m", "COBY"] + coby_args
        else:
            coby_args = substitute_protein_pdb(args.coby_args, coby_input_pdb)
            coby_cmd = ["python3", "-m", "COBY"] + coby_args
        run(coby_cmd, shell=False)
        coby_pdb = Path('COBY.pdb')
        assert coby_pdb.exists(), "COBY membrane build failed"
        # 9. Add TER to COBY, then rename POPC and chains
        make_coby_ter_pdb('COBY.pdb', str(coby_input_pdb), 'COBY_TER.pdb')
        run(f'python Scripts/renamePOPC.py')
        run(f'python Scripts/renameChain.py')
        # Always add TER after POPC/chain renaming if needed
        # You can re-use the existing add_ter.py for this purpose if needed
    else:
        # Not COBY: (If you want to split POPC for standalone systems, you could add it here)
        pass

    # 10. Run tleap
    run('tleap -f tleap.in')
    prmtop, inpcrd = autodetect_amber_files()

    # 11. Run convert_and_split
    run(f'python Scripts/convert_and_split.py {prmtop} {inpcrd}')
    run(f'python Scripts/resolvate_replicas.py --base {base}_AABY --replicas {args.replicas}')
    # if args.replicas == 1:
    #     run(f'python Scripts/prepare_for_simulations.py {base}_AABY')
    # else:
    #     for rep in range(1, args.replicas + 1):
    #         ros.chdir(f'r{rep}')
    #         run(f'python Scripts/prepare_for_simulations.py {base}_AABY')
    #         os.chdir("..")

    total_time = time.time() - start_time
    print(f"\n✅ All done. Your system is ready for GROMACS.")
    print(f"Total build time: {total_time:.1f} seconds ({total_time/60:.2f} minutes)")

if __name__ == "__main__":
    main()
