#!/usr/bin/env python3

import subprocess
import sys

# ---- USER SETTINGS ----
input_gro = "system.gro"
input_top = "topol.top"
mdp = "../mdps/step6.0_minimization.mdp"
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

REMOVE_RESNAMES = ["WAT", "SOL", "Na+", "Cl-"]

def run(cmd, return_output=False):
    print("[RUN]", cmd)
    if return_output:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        # result.stdout is the output, result.stderr is error output
        return result.stdout
    else:
        subprocess.run(cmd, shell=True, check=True)

def remove_ions_and_water_gro(input_gro, output_gro, remove_resnames=REMOVE_RESNAMES):
    with open(input_gro) as f:
        lines = f.readlines()
    new_lines = [lines[0], lines[1]]  # Title and atom count (will fix atom count at end)
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
    print(f"[INFO] Wrote {output_gro} (no ions/water), atoms: {natoms}")

def resize_box(input_gro, output_gro):
    z_length = float(run(f"tail -1 {input_gro} | awk '{{print $3}}'", return_output=True).strip())
    run(f'gmx editconf -f {input_gro} -d -.15 -o {output_gro}')
    new_z_length = float(run(f"tail -1 {output_gro} | awk '{{print $3}}'", return_output=True).strip())
    last_line_nr = int(run(f"wc -l {output_gro} | awk '{{print $1}}'", return_output=True).strip())
    print(z_length, new_z_length)
    run(f' sed -i "{last_line_nr}s/{new_z_length}.*/{z_length}/" {output_gro}')


def remove_ions_and_water_top(input_top, output_top, remove_resnames=None):
    if remove_resnames is None:
        remove_resnames = ["WAT", "SOL", "Na+", "Cl-"]
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
                # If line is empty, comment, or next section header, finish section
                if not stripped or stripped.startswith("["):
                    in_molecules = False
                    f.write(line)
                    continue
                parts = line.split()
                if parts and parts[0].upper() in remove_resnames_set:
                    continue  # skip this molecule line entirely
                f.write(line)
            else:
                f.write(line)
    print(f"[INFO] Wrote {output_top} (removed molecules: {', '.join(remove_resnames)})")


def get_pc_n31_box(gro_file):
    minz = 1e9
    maxz = -1e9
    with open(gro_file) as f:
        for l in f:
            if len(l) > 44 and l[5:10].strip() == "PC" and l[10:15].strip() == "N31":
            # if len(l) > 44 and l[5:10].strip() == "CHL" and l[10:15].strip() == "O1":
                z = float(l[36:44])
                minz = min(minz, z)
                maxz = max(maxz, z)
    print(f"[INFO] PC N31 box min: {minz}, max: {maxz}")
    return minz, maxz

def in_box(x, y, z, minz, maxz):
    return all(mn <= v <= mx for v, mn, mx in zip([x, y, z], minz, maxz))

def remove_sol_in_zrange(input_gro, output_gro, minz, maxz, pad=0.0):
    """
    Remove entire SOL (water) residues if ANY atom of the residue has a z-coordinate within [minz-pad, maxz+pad].
    """
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
                atom_l = lines[j]
                z = float(atom_l[36:44])
                if minz <= z <= maxz:
                    remove = True
                atoms.append(atom_l)
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
                #if line.strip() == "" or line.strip().startswith(";") or line.strip().startswith("["):
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
    print(f"[INFO] Updated SOL count in {output_top}")

def rename_sol_to_wat_gro(input_gro, output_gro):
    with open(input_gro) as f:
        lines = f.readlines()
    new_lines = []
    for l in lines:
        if len(l) > 10 and l[5:10].strip() == "SOL":
            l = l[:5] + "WAT" + l[8:]
        new_lines.append(l)
    with open(output_gro, "w") as f:
        f.writelines(new_lines)
    print(f"[INFO] Wrote {output_gro} with SOL→WAT")

def rename_sol_to_wat_top(input_top, output_top):
    with open(input_top) as f:
        lines = f.readlines()
    with open(output_top, "w") as f:
        for line in lines:
            parts = line.split()
            if parts and parts[0] == "SOL":
                line = line.replace("SOL", "WAT", 1)
            f.write(line)
    print(f"[INFO] Wrote {output_top} with SOL→WAT")

def make_ndx_for_wat(grofile, ndxfile):
    # Use GROMACS tool
    cmd = f"gmx make_ndx -f {grofile} -o {ndxfile}"
    process = subprocess.Popen(cmd.split(), stdin=subprocess.PIPE, text=True)
    process.communicate("r WAT\nq\n")
    print(f"[INFO] Created {ndxfile} with group WAT")

def run_softcore_minimization(input_gro, input_top, min_mdp, output_tpr):
    run(f'export GMX_MAXCONSTRWARN=-1; gmx grompp -f {min_mdp} -r {input_gro} -c {input_gro} -p {input_top} -o {output_tpr} -maxwarn 3; gmx mdrun -deffnm softcore -v')
    run(f'gmx grompp -f {mdp} -r softcore.gro -c softcore.gro -p {input_top} -o membrane.tpr -maxwarn 3')
    run(f'gmx mdrun -deffnm membrane -v')
    run(f'unset GMX_MAXCONSTRWARN')




def main(base="system"):
    final_gro = f"{base}.gro"
    final_top = f"{base}.top"
    # 1. Remove all ions and water
    remove_ions_and_water_gro(input_gro, intermediate_gro)
    remove_ions_and_water_top(input_top, intermediate_top)
    # 2. Resize box
    resize_box(intermediate_gro, resized_gro)
    # 3. Soft-core minimization of membrane/protein complex
    run_softcore_minimization(resized_gro, intermediate_top, min_mdp, "softcore.tpr")
    # 4. Re-solvate (add water molecules)
    run(f"gmx solvate -cp {minimized_gro} -cs tip4p.gro -o {resolvated_gro} -p {intermediate_top}")
    # 5. Remove SOLs in PA C116 box
    minz, maxz = get_pc_n31_box(resolvated_gro)
    remove_sol_in_zrange(resolvated_gro, cleaned_gro, minz, maxz)
    # 6. Count SOLs and update topology
    num_sol = count_sol_molecules(cleaned_gro)
    update_sol_count_top(intermediate_top, cleaned_top, num_sol)
    print(f"[INFO] Final number of SOL: {num_sol}")
    # 7. Rename SOL to WAT in gro and top
    rename_sol_to_wat_gro(cleaned_gro, gro_wat)
    rename_sol_to_wat_top(cleaned_top, top_wat)
    # 8. Make index for WAT
    make_ndx_for_wat(gro_wat, ndx)
    # 9. Create tpr and run genion
    run(f"gmx grompp -f {mdp} -r {gro_wat} -c {gro_wat} -p {top_wat} -o {tpr} -maxwarn 3")
    run(f'echo WAT | gmx genion -s {tpr} -p {top_wat} -pname Na+ -nname Cl- -neutral -conc 0.15 -o {final_gro} -n {ndx}')
    run(f'cp {top_wat} {final_top}')
    print(f"[DONE] See {final_gro} and {final_top} for output.")

if __name__ == "__main__":
    main(sys.argv[1])
