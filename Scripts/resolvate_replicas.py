#!/usr/bin/env python3

import subprocess
import sys
import os
import shutil
import argparse

# ---- USER SETTINGS ----
input_gro = "system.gro"
input_top = "topol.top"
input_top = "topol_nowater.top"
default_mdp = "../mdps/step6.0_minimization.mdp"
min_mdp = "../mdps/step5_mini.mdp"
# -----------------------

intermediate_gro = "noions.gro"
intermediate_top = "noions.top"
resized_gro = "resized.gro"
minimized_gro = "membrane.gro"
resolvated_gro = "resolvated.gro"
resolvated_top = "resolvated.top"
cleaned_gro = "system_solv.gro"
cleaned_top = "topol_updated.top"
gro_wat = "system_solv_wat.gro"
top_wat = "topol_updated_wat.top"
ndx = "solv.ndx"
tpr = "ions.tpr"

REMOVE_RESNAMES = ["WAT", "SOL", "Na", "Cl","K", "Ca", "Mg", "Zn"]

def run(cmd, return_output=False):
    print("[RUN]", cmd)
    if return_output:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return result.stdout
    else:
        subprocess.run(cmd, shell=True, check=True)

def remove_ions_and_water_gro(input_gro, output_gro, remove_resnames=None):
    if remove_resnames is None:
        remove_resnames = ["WAT", "SOL", "Na", "Cl","K", "Ca", "Mg", "Zn"]
    remove_resnames_set = set([x.upper() for x in remove_resnames])
    with open(input_gro) as f:
        lines = f.readlines()
    new_lines = [lines[0], lines[1]]
    natoms = 0
    for l in lines[2:-1]:
        if len(l) < 10:
            continue
        resname = l[5:10].strip()
        if resname not in remove_resnames:
            new_lines.append(l)
            natoms += 1
    new_lines.append(lines[-1])
    new_lines[1] = f"{natoms}\n"
    with open(output_gro, "w") as f:
        f.writelines(new_lines)

def resize_box(input_gro, output_gro):
    z_length = float(run(f"tail -1 {input_gro} | awk '{{print $3}}'", return_output=True).strip())
    run(f'gmx editconf -f {input_gro} -d -.15 -o {output_gro}')
    new_z_length = float(run(f"tail -1 {output_gro} | awk '{{print $3}}'", return_output=True).strip())
    last_line_nr = int(run(f"wc -l {output_gro} | awk '{{print $1}}'", return_output=True).strip())
    run(f'sed -i "{last_line_nr}s/{new_z_length}.*/{z_length}/" {output_gro}')
    run(f'gmx editconf -f {output_gro} -c -o {output_gro}')

def remove_ions_and_water_top(input_top, output_top, remove_resnames=None):
    if remove_resnames is None:
        remove_resnames = ["WAT", "SOL", "Na", "Cl","K", "Ca", "Mg", "Zn"]
    remove_resnames_set = set([x.upper() for x in remove_resnames])
    with open(input_top) as f:
        lines = f.readlines()
    with open(output_top, "w") as f:
        in_molecules = False
        for line in lines:
            stripped = line.strip()
            if stripped.lower().startswith("[ molecules ]"):
                in_molecules = True
                f.write(line)
                continue
            if in_molecules:
                if not stripped or stripped.startswith("["):
                    in_molecules = False
                    f.write(line)
                    continue
                parts = line.split()
                if parts and parts[0].upper() in remove_resnames_set:
                    continue
                f.write(line)
            else:
                f.write(line)

def get_pc_n31_box(gro_file):
    minz = 1e9
    maxz = -1e9
    with open(gro_file) as f:
        for l in f:
            #if len(l) > 44 and l[5:10].strip() == "PC" and l[10:15].strip() == "N31": #PE - N31, PS - N31, PH- -P31, PC - N31, PGR - C32, SPM - N31
            if len(l) > 44 and (l[5:10].strip(),l[10:15].strip()) in [("PE", 'N31'), ("PS", 'N31'),("PC", 'N31'), ("SPM", 'N31', ("PH-", 'P31'), ("PGR", 'C32'))]: #PE - N31, PS - N31, PH- -P31, PC - N31, PGR - C32, SPM - N31
                z = float(l[36:44])
                minz = min(minz, z)
                maxz = max(maxz, z)
    print(f"[INFO] PC N31 box min: {minz:.3f}, max: {maxz:.3f}")
    return minz, maxz

def remove_sol_in_zrange(input_gro, output_gro, minz, maxz, pad=0.0):
    minz -= pad
    maxz += pad
    removed = 0
    with open(input_gro) as f:
        lines = f.readlines()
    new_lines = [lines[0], lines[1]]
    natoms = 0
    i = 2
    while i < len(lines)-1:
        l = lines[i]
        if len(l) < 44:
            new_lines.append(l)
            i += 1
            continue
        resname = l[5:10].strip()
        resid = l[0:5]
        if resname == "SOL":
            j = i
            remove = False
            atoms = []
            while j < len(lines)-1 and lines[j][0:5] == resid and lines[j][5:10].strip() == "SOL":
                z = float(lines[j][36:44])
                if minz <= z <= maxz:
                    remove = True
                atoms.append(lines[j])
                j += 1
            if not remove:
                new_lines.extend(atoms)
                natoms += len(atoms)
            else:
                removed += 1
            i = j
        else:
            new_lines.append(l)
            natoms += 1
            i += 1
    new_lines.append(lines[-1])
    new_lines[1] = f"{natoms}\n"
    with open(output_gro, "w") as f:
        f.writelines(new_lines)
    print(f"[INFO] Wrote {output_gro} (removed {removed} SOL molecules in z=[{minz:.2f}, {maxz:.2f}])")

def count_sol_molecules(gro_file):
    count = 0
    last_resid = None
    with open(gro_file) as f:
        for line in f:
            if len(line) > 10 and line[5:10].strip() == "SOL":
                resid = line[0:5]
                if resid != last_resid:
                    count += 1
                    last_resid = resid
    return count

def update_sol_count_top(input_top, output_top, new_sol_count):
    with open(input_top) as f:
        lines = f.readlines()
    with open(output_top, "w") as f:
        in_molecules = False
        for line in lines:
            if line.strip().lower().startswith("[ molecules ]"):
                in_molecules = True
                f.write(line)
                continue
            if in_molecules:
                if line.strip() == "" or line.strip().startswith("["):
                    in_molecules = False
                    f.write(line)
                    continue
                parts = line.split()
                if parts and parts[0] == "SOL":
                    f.write(f"SOL     {new_sol_count}\n")
                else:
                    f.write(line)
            else:
                f.write(line)

def rename_sol_to_wat_gro(input_gro, output_gro):
    with open(input_gro) as f:
        lines = f.readlines()
    with open(output_gro, "w") as f:
        for l in lines:
            if len(l) > 10 and l[5:10].strip() == "SOL":
                l = l[:5] + "WAT" + l[8:]
            f.write(l)

def rename_sol_to_wat_top(input_top, output_top):
    with open(input_top) as f:
        lines = f.readlines()
    with open(output_top, "w") as f:
        for line in lines:
            parts = line.split()
            if parts and parts[0] == "SOL":
                line = line.replace("SOL", "WAT", 1)
            f.write(line)

def make_ndx_for_wat(grofile, ndxfile):
    cmd = f"gmx make_ndx -f {grofile} -o {ndxfile}"
    process = subprocess.Popen(cmd.split(), stdin=subprocess.PIPE, text=True)
    process.communicate("r WAT\nq\n")
    print(f"[INFO] Created {ndxfile} with group WAT")

def run_softcore_minimization(input_gro, input_top, min_mdp, output_tpr):
    run(f'export GMX_MAXCONSTRWARN=-1; gmx grompp -f {min_mdp} -r {input_gro} -c {input_gro} -p {input_top} -o {output_tpr} -maxwarn 3; gmx mdrun -deffnm softcore -v')
    run(f'gmx grompp -f {default_mdp} -r softcore.gro -c softcore.gro -p {input_top} -o membrane.tpr -maxwarn 3')
    run(f'gmx mdrun -deffnm membrane -v')
    run(f'unset GMX_MAXCONSTRWARN')

def prepare_once(base="system"):
    #remove_ions_and_water_gro(input_gro, intermediate_gro)
    #remove_ions_and_water_top(input_top, intermediate_top)
    run(f'cp {input_gro} {intermediate_gro}')
    run(f'cp {input_top} {intermediate_top}')
    resize_box(intermediate_gro, resized_gro)
    run_softcore_minimization(resized_gro, intermediate_top, min_mdp, "softcore.tpr")
    print("[INFO] Preparation done. membrane.gro and noions.top ready.")

def fix_topol():
    protein = ['CYS','ASP','SER','GLN','LYS','ILE','PRO','THR','PHE','ASN','GLY','HIS','LEU','ARG','TRP','ALA','VAL','GLU','TYR','MET', 'GLH','NME','HIE','ACE']
    lipid = ['PA', 'PH-', 'AR', 'PS', 'ST', 'MY', 'PC', 'DHA', 'SA', 'SPM', 'LAL', 'PGR', 'OL', 'PE']
    nucleic = ['A','C','G', 'U', 'T']
    lipid_residues = {"DAPC": ["AR", "PC", "AR"],
    "DLPC": ["LAL", "PC", "LAL"],
    "DLPG": ["LAL", "PGR", "LAL"],
    "DMPC": ["MY", "PC", "MY"],
    "DMPG": ["MY", "PGR", "MY"],
    "DOPC": ["OL", "PC", "OL"],
    "DOPG": ["OL", "PGR", "OL"],
    "DOPS": ["OL", "PS", "OL"],
    "DPPC": ["PA", "PC", "PA"],
    "DPPG": ["PA", "PGR", "PA"],
    "DSPC": ["ST", "PC", "ST"],
    "DSPG": ["ST", "PGR", "ST"],
    "POPA": ["PA", "PH-", "OL"],
    "POPC": ["PA", "PC", "OL"],
    "POPE": ["PA", "PE", "OL"],
    "POPG": ["PA", "PGR", "OL"],
    "POPS": ["PA", "PS", "OL"],
    "PSM": ["PA", "SPM", "SA"],
    "SDPC": ["ST", "PC", "DHA"],
    "SSM": ["ST", "SPM", "SA"]}
    run(f'mv ff.itp toppar/')
    run(f'sed -i "s/ff/toppar\/ff/" topol.top')
    with open('topol.top') as f:
        lines = f.readlines()
    systems = []
    for line in lines:
        l = line.split()
        if len(l) == 0:
            break
        if l[0] == '#include':
            if l[1].split('/')[1][:6] == 'system':
                print(l[1].split('/')[1][:-1])
                systems.append(l[1][1:-1])
    system_types = {}
    number_protein = 1
    number_nucleic = 1
    for idx, system in enumerate(systems):
        #print(idx, system)
        with open(system) as f:
            lines = f.readlines()
            atoms = False
            residues = []
            for line in lines:
                # print(line[:6] )
                if atoms:
                    if line.startswith('\n'):
                        atoms = False
                        residues = sorted(list(set(residues)))
                        for l, r in lipid_residues.items():
                            #print(l,r, residues, sorted(set(r)))
                            if residues == sorted(set(r)):
                                system_types[idx+1] = l
                        break
                    res = line.split()[3]
                    if res in protein:
                        system_types[idx+1] = f'protein{number_protein}'
                        number_protein += 1
                        break
                    elif res in nucleic:
                        system_types[idx+1] = f'nucleic{number_nucleic}'
                        number_nucleic += 1
                        break
                    else:
                        residues.append(res)
                if line[:6] == ';   nr':
                    atoms = True
    for system, system_type in system_types.items():
        print(system, system_type)
        run(f'mv toppar/system{system}.itp toppar/{system_type}.itp')
        run(fr'sed -i "s/system{system}\./{system_type}\./g" toppar/{system_type}.itp')
        run(fr'sed -i "s/system{system}\./{system_type}\./g" topol.top')
        run(f'sed -i "s/system{system} /{system_type} /g" toppar/{system_type}.itp')
        run(f'sed -i "s/system{system} /{system_type} /g" topol.top')



def resolvate_only(base="system", mdp="mdps/step6.0_minimization.mdp", water='OPC', ions="Na+,Cl-", conc="0.15"):
    water_type = {'OPC': 'tip4p', 'TIP3P': 'spc216', 'TIP4PEW': 'tip4p'}
    run(f"gmx solvate -cp membrane.gro -cs {water_type[water.upper()]}.gro -o {resolvated_gro} -p noions.top") ## change water model to user specified
    minz, maxz = get_pc_n31_box(resolvated_gro)
    remove_sol_in_zrange(resolvated_gro, cleaned_gro, minz, maxz)
    num_sol = count_sol_molecules(cleaned_gro)
    update_sol_count_top("noions.top", cleaned_top, num_sol)
    rename_sol_to_wat_gro(cleaned_gro, gro_wat)
    rename_sol_to_wat_top(cleaned_top, top_wat)
    make_ndx_for_wat(gro_wat, ndx)
    pion, nion = ions.split(',')
    run(f"gmx grompp -f {mdp} -r {gro_wat} -c {gro_wat} -p {top_wat} -o {tpr} -maxwarn 3")
    run(f'echo WAT | gmx genion -s {tpr} -p {top_wat} -pname {pion} -nname {nion} -neutral -conc {conc} -o {base}.gro -n {ndx}') ## Change to ions as per user preference
    run(f"gmx grompp -f {mdp} -r {base}.gro -c {base}.gro -p {top_wat} -o pbc.tpr -maxwarn 3")
    run(f'echo 0 | gmx trjconv -f {base}.gro -s pbc.tpr -o {base}.gro -pbc whole')
    run(f'cp {top_wat} topol.top')
    fix_topol()
    print(f"[DONE] See {base}.gro and {base}.top for output.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicas", type=int, default=1, help="Number of replicas to run (default: 1)")
    parser.add_argument("--base", type=str, default="system", help="Base name for final output (default: system)")
    parser.add_argument('--water', default='OPC', help='Which water model to use (default: OPC) (Options: OPC, TIP3P, TIP4PEW)')
    parser.add_argument('--ions', default='Na+,Cl-', help='Which ions to use for solvation')
    parser.add_argument('--conc', default='0.15', help='Which ion concentration to build')
    parser.add_argument('--Z', default='10', help='Which Z height the box needs')
    args = parser.parse_args()

    if args.replicas == 1:
        prepare_once(base=args.base)
        resolvate_only(base=args.base, mdp=default_mdp, water=args.water, ions=args.ions, conc=args.conc)
    else:
        prepare_once(base=args.base)
        for i in range(1, args.replicas + 1):
            rdir = f"r{i}"
            os.makedirs(rdir, exist_ok=True)
            shutil.copy("membrane.gro", os.path.join(rdir, "membrane.gro"))
            shutil.copy("noions.top", os.path.join(rdir, "noions.top"))

            # Copy mdps
            shutil.copytree("../mdps", os.path.join(rdir, "mdps"), dirs_exist_ok=True)

            # Copy toppar (force field includes)
            shutil.copytree("toppar", os.path.join(rdir, "toppar"), dirs_exist_ok=True)

            # Copy ff.itp if present
            if os.path.exists("ff.itp"):
                shutil.copy("ff.itp", os.path.join(rdir, "ff.itp"))

            print(f"\n[INFO] Running replica {i} in {rdir}")
            os.chdir(rdir)
            try:
                local_mdp = "mdps/step6.0_minimization.mdp"
                resolvate_only(base=args.base, mdp=local_mdp, water=args.water, ions=args.ions, conc=args.conc)
            finally:
                os.chdir("..")
