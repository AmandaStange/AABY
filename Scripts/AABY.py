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

def add_ter_rename_chains(pdb, out_ter_rename_chains):

    run(f"gmx editconf -f {pdb} -o {str(pdb).split('.')[0]}_gmx.pdb")
    with open(f'{str(pdb).split('.')[0]}_gmx.pdb') as f:
        lines = f.readlines()
    resids = []
    resids_unique = []
    chainids = []
    new_chain = []
    TERs = []
    idxs = []
    isATOM = False
    i = 0
    for idx, line in enumerate(lines):
        if line[:4] in ['ATOM', 'HETA']:
            resid = int(line[22:27])
            chainid = line[21]
            if i == 0:
                resids_unique.append(resid)
                new_chain.append([chainid, idx])
                i += 1
            if resids_unique[-1] != resid:
                resids_unique.append(resid)
            if new_chain[-1][0] != chainid:
                new_chain.append([chainid, idx])
            resids.append(resid)
            chainids.append(chainid)
            idxs.append(idx)

        if line[:3] in ['TER']:
            TERs.append(idx)
    new_TERs = []

    for i in range(1,len(resids)):
        if abs(resids[i] - resids[i-1]) > 1:
            if idxs[i] - idxs[i - 1] != 2:

                new_TERs.append(idxs[i] + len(new_TERs))
                print(resids[i-1], resids[i], chainids[i-1], chainids[i], idxs[i-1], idxs[i])

    for idx in new_TERs:
        lines.insert(idx, 'TER\n')



    TERs = []
    for idx, line in enumerate(lines):
        if line == 'TER\n':
            TERs.append(idx)



    chainids = {}
    chains = range(65, 65+28)
    TERs = [0] + TERs

    for idx in range(len(TERs)-1):
        chainids[chr(chains[idx])] = range(TERs[idx]+1, TERs[idx+1]+1)


    new_lines = ''
    for idx, line in enumerate(lines):
        if line[:4] in ['ATOM', 'HETA']:
            for k, v in chainids.items():
                if idx in v:
                    chain = k
                    break
            new_lines += f'{line[:21]}{chain}{line[22:]}'
        else:
            new_lines += line
    with open(out_ter_rename_chains, "w") as f:
        f.write("".join(new_lines))





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


def auto_detect_types(water_model='OPC'):
    #input4amber.pdb
    # source leaprc.protein.ff19SB
    # source leaprc.water.opc
    # source leaprc.lipid21
    #sed -i '1s/^/task goes here\n/' todo.txt

    forcefields = {'protein': 'protein.ff19SB', 'lipid': 'lipid21', 'dna': 'DNA.OL24', 'rna': 'RNA.OL3', 'water': {'OPC': 'water.opc', 'TIP3P': 'water.tip3p', 'TIP4PEW': 'water.tip4pew'}}

    protein = ['CYS','ASP','SER','GLN','LYS','ILE','PRO','THR','PHE','ASN','GLY','HIS','LEU','ARG','TRP','ALA','VAL','GLU','TYR','MET', 'GLH','NME','HIE','ACE', 'HID', 'CYX']
    lipid = ['PC', 'PA', 'OL', 'CHL']
    lipid = ['PA', 'PH-', 'AR', 'PS', 'ST', 'MY', 'PC', 'DHA', 'SA', 'SPM', 'LAL', 'PGR', 'OL', 'PE']
    nucleic = ['A','C','G']
    rna = ['U']
    dna = ['T']
    with open('input4amber.pdb', 'r') as f:
        lines = f.readlines()

    resnames = []

    ff_types = []

    for line in lines:

        l = line.split()

        if l[0] == 'END':
            break

        if l[0] in ['ATOM', 'HETATM']:
            resnames.append(l[3])

    resnames = list(set(resnames))



    for resname in resnames:
        if resname in protein:
            if 'protein' not in ff_types:
                ff_types.append('protein')
        elif resname in lipid:
            if 'lipid' not in ff_types:
                ff_types.append('lipid')
        # elif resname in nucleic:
        #     if 'nucleic' not in ff_types:
        #         ff_types.append('nucleic')
        elif resname in rna:
            if 'rna' not in ff_types:
                ff_types.append('rna')
        elif resname in dna:
            if 'dna' not in ff_types:
                ff_types.append('dna')
        else:
            print('type not found', resname, ff_types)

    leaprc = ''

    for ff_type in sorted(ff_types, key=lambda x: len(x), reverse=True):
        leaprc += f'source leaprc.{forcefields[ff_type]}\\n'

    leaprc += f'source leaprc.{forcefields['water'][water_model.upper()]}\\n'

    return leaprc

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

def infer_element(line):
    TWO_LETTER_ELEMS = {
    "Cl","Br","Na","Mg","Al","Si","Li","Be","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu",
    "Zn","Ga","Ge","As","Se","Kr","Rb","Sr","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "In","Sn","Sb","Te","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy",
    "Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi",
    "Po","At","Rn"
}
    """Prefer PDB element field (cols 77–78). Else infer from atom name."""
    elem = line[76:78].strip() if len(line) >= 78 else ""
    if elem:
        e = elem.capitalize()
        return e if e in TWO_LETTER_ELEMS else e[:1]
    # Fallback from atom name field (cols 13–16)
    an = (line[12:16] if len(line) >= 16 else "").strip()
    if not an:
        return "X"
    if an[:2].capitalize() in ("Cl","Br"):
        return an[:2].capitalize()
    if an[:2].capitalize() in TWO_LETTER_ELEMS:
        return an[:2].capitalize()
    return an[0].upper()

def format_pdb_atom_name(name, element):
    """Return a 4-char atom-name field with PDB alignment rules."""
    if len(element) == 2:
        return name.ljust(4)[:4]      # 2-letter element -> left-justified
    else:
        return name.rjust(4)[-4:]      # 1-letter element -> right-justified

def rename_pdb4antechamber(
    inp, outp, cap=4, include_resnames=None, case_insensitive=True,
    renumber=False, renumber_start=1, update_conect=False
):
    """
    Rename atom names uniquely (per-residue) only for residues in include_resnames,
    and optionally renumber all atom serials in the full PDB (updates CONECT).

    Args
    ----
    inp : str        input PDB
    outp: str        output PDB
    cap : int        max atom-name length (default 4)
    include_resnames: {str} or list/tuple/set of resnames to affect; None = all
    case_insensitive: bool  compare resnames case-insensitively
    renumber : bool  renumber serials across entire file
    renumber_start : int    starting serial
    update_conect : bool    update CONECT lines to new serials when renumbering
    """
    # normalize filter set
    if include_resnames is None:
        filter_set = None
    else:
        rs = include_resnames if isinstance(include_resnames, (list, tuple, set)) else [include_resnames]
        filter_set = set(r.upper() if case_insensitive else r for r in rs)

    # First pass: rename (only filtered residues), collect lines
    out_lines = []
    counts = {}       # per-residue & element counters
    model_id = 0

    def resname_pass(line):
        rn = line[17:20]
        return rn.strip().upper() if case_insensitive else rn.strip()

    for raw in open(inp, "r"):
        line = raw.rstrip("\n")
        rec = line[:6]

        if rec.startswith("MODEL"):
            # try parse model number, else increment
            try:
                model_id = int(line[10:14].strip())
            except Exception:
                model_id += 1
            out_lines.append(line)
            continue

        if rec.startswith("ATOM  ") or rec.startswith("HETATM"):
            # ensure padding for safe slicing
            if len(line) < 80:
                line = line + " " * (80 - len(line))

            resname = resname_pass(line)
            do_rename = (filter_set is None) or (resname in filter_set)

            if do_rename:
                chain   = line[21]
                resseq  = line[22:26]
                icode   = line[26]
                elem    = infer_element(line)

                # per-residue + element key
                key = (chain, resseq, icode, resname, model_id, elem)
                counts[key] = counts.get(key, 0) + 1

                new = f"{elem}{counts[key]}"
                if len(new) > cap:
                    new = elem + str(counts[key])[: max(1, cap - len(elem))]

                name_field = format_pdb_atom_name(new, elem)
                line = line[:12] + name_field + line[16:]

                # normalize element field
                elem_field = elem.rjust(2)[:2]
                line = line[:76] + elem_field + line[78:]

            out_lines.append(line)
        else:
            out_lines.append(line)

    # If no renumbering, write and return
    if not renumber:
        with open(outp, "w") as g:
            g.write("\n".join(out_lines) + ("\n" if out_lines and out_lines[-1] != "" else ""))
        return

    # Second pass: renumber serials (ATOM/HETATM) and update CONECT
    # Build mapping old_serial -> new_serial
    serial_map = {}
    next_serial = renumber_start
    renum_lines = []

    for line in out_lines:
        rec = line[:6]
        if rec.startswith("ATOM  ") or rec.startswith("HETATM"):
            if len(line) < 80:
                line = line + " " * (80 - len(line))
            old_serial = int(line[6:11])
            serial_map[old_serial] = next_serial
            # write new serial into cols 7-11
            new_serial_field = f"{next_serial:5d}"
            line = line[:6] + new_serial_field + line[11:]
            next_serial += 1
        renum_lines.append(line)

    # Update CONECT if requested
    if update_conect:
        final_lines = []
        for line in renum_lines:
            if line.startswith("CONECT"):
                # CONECT + up to five serials in 5-wide fields
                # columns: 7-11 (atom) then 12-16, 17-21, 22-26, 27-31
                # We'll parse ints and rewrite via mapping (skip if not present)
                fields = [line[:6]]  # "CONECT"
                nums = []
                # grab 5-wide fields after col 6
                for i in range(6, len(line), 5):
                    chunk = line[i:i+5]
                    if chunk.strip().isdigit():
                        nums.append(int(chunk))
                    else:
                        nums.append(None)
                if nums:
                    mapped = []
                    for n in nums:
                        if n is None:
                            mapped.append("     ")
                        else:
                            mapped.append(f"{serial_map.get(n, n):5d}")
                    line = fields[0] + "".join(mapped)
                    # ensure newline
                final_lines.append(line)
            else:
                final_lines.append(line)
    else:
        final_lines = renum_lines

    with open(outp, "w") as g:
        g.write("\n".join(final_lines) + ("\n" if final_lines and final_lines[-1] != "" else ""))



def antechamber(mol2=None, nc=None, input_pdb=None):
    if nc is None:
        sys.exit("--nc option is not set. This is required for using the --antechamber option!")
    ## Lines taken from amber_geostad gcif_to_mol2
    #res = mol2.split('_')[0]
    run(f'grep {mol2} {input_pdb} > {mol2}.pdb')
    # rename_pdb4antechamber(f'{mol2}.pdb', f'{mol2}_fix.pdb')
    #run(f'gmx editconf -f {mol2}.pdb -resnr 1 -o {mol2}.pdb')
    try:
        run(f'obabel {mol2}.pdb -O {mol2}_obabel.mol2')
    except:
        sys.exit("obabel not installed. This is required for using the --antechamber option! Run 'sudo apt install openbabel' to continue")
    # run(f'sed -i "s/{mol2}.pdb/{mol2}/" {mol2}_obabel.mol2')
    # run(f'sed -i "s/{mol2}1/{mol2} /" {mol2}_obabel.mol2')
    run(f"$AMBERHOME/bin/antechamber -i {mol2}_obabel.mol2 -fi mol2 -o {mol2}.mol2 -fo mol2 -bk comp_{mol2} -s 0 -dr no -nc {nc} -at gaff2 -c abcg2 -ek 'qm_theory=\"AM1\", maxcyc=1000, ndiis_attempts=700,'")
    run(f"$AMBERHOME/bin/parmchk2 -s 2 -i {mol2}.mol2 -f mol2 -o {mol2}.frcmod")
    run(f'/bin/rm -f ANTECH* ATOMTYP* sqm.*')
    tleap_ante = f"source leaprc.gaff2 \\n{mol2} = loadMol2 {mol2}.mol2 \\nloadAmberParams {mol2}.frcmod\\n"
    run(f'sed -i "1s/^/{tleap_ante}/" tleap.in')
    run(f'sed -i "1s/^/{tleap_ante}/" tleap_solv.in')



def main():
    parser = argparse.ArgumentParser(description="AABY: End-to-end AMBER system builder from PDB")
    parser.add_argument('-f', '--input', required=True, help='Input PDB file')
    parser.add_argument('--chains', default="A", help='Comma-separated list of protein chains for ACE/NME')
    parser.add_argument('-r', '--replicas', type=int, default=1, help='Number of resolvated replicas')
    parser.add_argument('--ssbond', action='store_true', help='Renumber SSBOND.txt, convert CYS to CYX, insert ssbonds into tleap.in')
    parser.add_argument('--water', default='OPC', help='Which water model to use (default: OPC) (Options: OPC, TIP3P, TIP4PEW)')
    parser.add_argument('--ions', default='Na+,Cl-', help='Which ions to use for solvation (first positive then negative)')
    parser.add_argument('--conc', default='0.15', help='Which ion concentration to build')
    parser.add_argument('--antechamber', default=False, help='Option to use antechamber to parameterise ligands')
    parser.add_argument('--mol2', default=None, help='Mol2 filename (without extension) to specify the ligand to be parameterised')
    parser.add_argument('--nc', default=None, help='Net charge of ligand')
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

    tleap_solv_template = Path('Scripts/tleap_solv.in')
    tleap_solv_local = Path('tleap_solv.in')
    shutil.copy(tleap_solv_template, tleap_solv_local)

    pdb = Path(args.input)
    base = pdb.stem

    if args.antechamber:
        out_renamedligandatoms = pdb.with_name(base + '_renameligandatoms.pdb')
        rename_pdb4antechamber(pdb, out_renamedligandatoms, include_resnames=args.mol2)
        # 0. Add TERs and rename chains
        out_ter_rename_chains = pdb.with_name(base + '_breaks.pdb')
        add_ter_rename_chains(out_renamedligandatoms, out_ter_rename_chains)
    else:
        # 0. Add TERs and rename chains
        out_ter_rename_chains = pdb.with_name(base + '_breaks.pdb')
        add_ter_rename_chains(pdb, out_ter_rename_chains)




    # 1. Add caps
    out_ace_nme = pdb.with_name(base + '_breaks_ACE_NME.pdb')
    run(f'python Scripts/add_ace_nme.py {out_ter_rename_chains} {args.chains}')
    assert out_ace_nme.exists(), "add_ace_nme failed"

    # 2. Renumber
    out_renum = pdb.with_name(base + '_breaks_ACE_NME_re.pdb')
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
        insert_ssbonds_into_tleap('tleap_solv.in', 'tleap_SSBONDs.txt')
    else:
        working_pdb = out_renum

    # 4. Delete H (heavy atom only)
    out_heavy = working_pdb.with_name(working_pdb.stem + '_heavy.pdb')
    run(f'printf "del 0-100\n!a H*\n q\n" | gmx make_ndx -f {working_pdb} -o heavy.ndx')
    file_size = os.path.getsize('heavy.ndx')
    if file_size == 0:
        run(f'cp {working_pdb} {out_heavy}')
    else:
        run(f'gmx trjconv -f {working_pdb} -s {working_pdb} -n heavy.ndx -o {out_heavy}')
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
    run(f'sed -i "s/HSD/HID/g" {renamed_pdb}')
    run(f'sed -i "s/HSE/HIE/g" {renamed_pdb}')

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
        apls = {'CHL': .41, 'DLPC': 0.61, 'DMPC': 0.60, 'DPPC': 0.62, 'DSPC': 0.60, 'DOPC': 0.67, 'POPC': 0.64, 'POPE': 0.56, 'DLPG': 0.66, 'DMPG': 0.65, 'DPPG': 0.68, 'DSPG': 0.67, 'DOPG': 0.71, 'POPG': 0.68, 'DOPS': 0.65, 'POPS': 0.62, 'POPA': 0.64, 'DAPC': 0.72, 'SDPC': 0.65, 'PSM': 0.58, 'SSM': 0.57}
        if args.coby_yaml:
            with open(args.coby_yaml) as f:
                coby_args = yaml.safe_load(f)
            mol_import = []
            ##NEW
            for idx, arg in enumerate(coby_args):
                if isinstance(arg, str):
                    if arg == "-membrane":
                        coby_args.insert(idx + 1, 'optimize_run:False')
                        coby_args.insert(idx + 1, 'grid_maker_algorithm:no_groups')
                    if arg.split(':')[0] == 'lipid':
                        mol_import.append(arg.split(':')[1])
                        coby_args[idx] = coby_args[idx] + f":apl:{apls[arg.split(':')[1]]}:params:Amber"
            for mol in mol_import:
                coby_args.append('-molecule_import')
                coby_args.append(f'file:models/{mol}.pdb')
                coby_args.append(f'moleculetype:{mol}')
                coby_args.append('params:Amber')

            coby_args.append('-out_sys')
            coby_args.append('COBY.pdb')
            coby_args.append('-out_top')
            coby_args.append('COBY.top')

            with open('coby_input_params.yaml', 'w') as f:
                f.write(yaml.dump(coby_args))

            ## newparams ended

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
        run(f'cp COBY_TER_lipid21_chains.pdb input4amber.pdb')
        # Always add TER after POPC/chain renaming if needed
        # You can re-use the existing add_ter.py for this purpose if needed
    else:
        # Not COBY: (If you want to split POPC for standalone systems, you could add it here)
        run(f'cp {coby_input_pdb} input4amber.pdb')
        pass

    # 10. Run tleap

    # 10.a Run antechamber
    if args.antechamber:
        antechamber(mol2=args.mol2, nc=args.nc, input_pdb=str(out_renamedligandatoms))

    leap = auto_detect_types(water_model=args.water)
    run(f"sed -i '1s/^/{leap}/' tleap.in")
    run(f"sed -i '1s/^/{leap}/' tleap_solv.in")

    run('tleap -f tleap.in')
    prmtop, inpcrd, top_file, gro_file, topol = "system.prmtop","system.inpcrd", "system.top", "system.gro", "topol"
    run(f'python Scripts/convert_and_split.py {prmtop} {inpcrd} {top_file} {gro_file} {topol}')
    run(f'gmx editconf -f input4amber.pdb -o input4amber.gro; am=$(tail -n 1 input4amber.gro); sm=$(tail -n 1 system.gro); sed -i "s/$sm/$am/" system.gro')


    # 11. Create topology that includes solvent
    ## instert water and ions insert-molecules
    run('gmx editconf -f input4amber.pdb -o input4amber_solvX.pdb')
    nr_molecules = 0
    nr_atoms = 0 #'OPC': 'water.opc', 'TIP3P': 'water.tip3p', 'TIP4PEW': 'water.tip4pew'
    nr_atoms_water = {'OPC': 4, 'TIP3P': 3, 'TIP4PEW': 4}
    run(f'gmx insert-molecules -f input4amber_solvX.pdb -ci models/{args.water.lower()}.gro -o input4amber_solv{nr_molecules}.pdb -nmol 1')
    nr_molecules += 1
    nr_atoms += nr_atoms_water[args.water]
    for ion in args.ions.split(','):
        run(f'sed "s/XX /{ion}/g" models/ion.pdb > models/tmp.pdb')
        run(f'gmx insert-molecules -f input4amber_solv{nr_molecules-1}.pdb -ci models/tmp.pdb -o input4amber_solv{nr_molecules}.pdb -nmol 1')
        nr_molecules += 1
        nr_atoms += 1

    # with open(f"input4amber_solv{nr_molecules-1}.pdb", "rb") as f:
    #     after_insert = sum(1 for _ in f)

    # with open(f"input4amber_solvX.pdb", "rb") as f:
    #     before_insert = sum(1 for _ in f)

    difference_lines = nr_atoms + 2
    inserted_lines = []
    with open(f"input4amber_solv{nr_molecules-1}.pdb", "r") as f:
        lines = f.readlines()
        for line in lines[-difference_lines:]:
            inserted_lines.append(line)



    with open(f"input4amber.pdb", "r") as f:
        og_lines = f.readlines()
        original_length = len(og_lines)

    # print("OG LIENS", og_lines)

    with open('input4amber_solv.pdb','w') as f:
        for line in og_lines[:-2]:

            f.write(line)
        for line in inserted_lines:

            f.write(line)



    run('tleap -f tleap_solv.in')
    prmtop, inpcrd, top_file, gro_file, topol = "system_solv.prmtop","system_solv.inpcrd", "system_solv.top", "system_solv.gro", "topol_solv"
    run(f'python Scripts/convert_and_split.py {prmtop} {inpcrd} {top_file} {gro_file} {topol}')

    run(f'grep include topol_solv.top > topol_nowater.top')
    run(f'grep -v include topol.top >> topol_nowater.top')

    # 12. Resolvate
    run(f'python Scripts/resolvate_replicas.py --base {base}_AABY --replicas {args.replicas} --water {args.water} --ions {args.ions} --conc {args.conc}')
    if args.replicas == 1:
        if args.mol2 is None:
            run(f'python Scripts/prepare_for_simulations.py {base}_AABY')
        else:
            run(f'python Scripts/prepare_for_simulations.py {base}_AABY {args.mol2}')

    else:
        for rep in range(1, args.replicas + 1):
            os.chdir(f'r{rep}')
            if args.mol2 is None:
                run(f'python Scripts/prepare_for_simulations.py {base}_AABY')
            else:
                run(f'python Scripts/prepare_for_simulations.py {base}_AABY {args.mol2}')

            os.chdir("..")

    total_time = time.time() - start_time
    print(f"\n✅ All done. Your system is ready for GROMACS.")
    print(f"Total build time: {total_time:.1f} seconds ({total_time/60:.2f} minutes)")

if __name__ == "__main__":
    main()
