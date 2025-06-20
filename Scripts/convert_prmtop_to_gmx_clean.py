import sys
import os
import re
from parmed import load_file
from collections import defaultdict

def merge_itps(toppar_dir: str, output_ff_itp: str = "ff.itp") -> str:
    sections = [
        "defaults", "atomtypes", "bondtypes", "pairtypes",
        "angletypes", "dihedraltypes", "cmaptypes"
    ]
    section_data = {section: [] for section in sections}
    seen_lines = {section: set() for section in sections}
    for filename in sorted(os.listdir(toppar_dir)):
        if filename.endswith(".itp"):
            with open(os.path.join(toppar_dir, filename)) as f:
                lines = f.readlines()
            current_section = None
            for line in lines:
                stripped = line.strip()
                if stripped.startswith("[") and "]" in stripped:
                    header = re.findall(r"\[\s*(\w+)", stripped)
                    if header and header[0].lower() in section_data:
                        current_section = header[0].lower()
                        continue
                    else:
                        current_section = None
                if current_section and not stripped.startswith(";"):
                    key = line.strip()
                    if key and key not in seen_lines[current_section]:
                        section_data[current_section].append(line)
                        seen_lines[current_section].add(key)
    ff_itp_lines = []
    for section in sections:
        if section_data[section]:
            ff_itp_lines.append(f"[ {section} ]\n")
            ff_itp_lines.extend(section_data[section])
            ff_itp_lines.append("\n")
    output_path = os.path.join(toppar_dir, output_ff_itp)
    with open(output_path, "w") as f:
        f.writelines(ff_itp_lines)
    return output_path

def main(prmtop_file, inpcrd_file, toppar_dir):
    base = os.path.splitext(os.path.basename(prmtop_file))[0]
    amber = load_file(prmtop_file, inpcrd_file)
    amber.save(f"{base}.gro", format="gro", overwrite=True)
    amber.save(f"{base}.top", format="gromacs", overwrite=True)

    ff_path = merge_itps(toppar_dir)

    molecules = amber.split()
    mol_counts = defaultdict(int)
    itp_files = []
    posre_written = False

    for idx, (mol, _) in enumerate(molecules, 1):
        name = f"system{idx}"
        itp_name = f"{name}.itp"
        mol_counts[name] += 1

        mol.save(itp_name, format="gromacs", overwrite=True)
        with open(itp_name, "r") as f_in:
            lines = f_in.readlines()
        filtered = []
        inside_skip_section = False
        skip_sections = [
            '[ atomtypes ]', '[ bondtypes ]', '[ pairtypes ]', '[ angletypes ]',
            '[ dihedraltypes ]', '[ cmaptypes ]', '[ defaults ]'
        ]
        for line in lines:
            lower = line.strip().lower()
            if lower.startswith("[ system ]"):
                break
            if any(lower.startswith(sec) for sec in skip_sections):
                inside_skip_section = True
                continue
            if inside_skip_section and lower.startswith("[") and "]" in lower:
                inside_skip_section = False
            if not inside_skip_section:
                filtered.append(line)

        # Prepend the proper [ moleculetype ] section
        prepend = [
            "[ moleculetype ]\n",
            "; Name            nrexcl\n",
            f"{name:<16}3\n",
            "\n"
        ]
        filtered = prepend + filtered

        with open(itp_name, "w") as f_out:
            f_out.writelines(filtered)

        print(f"Wrote {itp_name}")
        itp_files.append(itp_name)

        if not posre_written:
            posre_lines = ["[ position_restraints ]\n; atom  type  fx  fy  fz\n"]
            for atom in mol.atoms:
                if atom.atomic_number > 1:
                    posre_lines.append(f"{atom.idx + 1:6d}    1   1000 1000 1000\n")
            with open("posre.itp", "w") as f:
                f.writelines(posre_lines)
            print(f"Wrote posre.itp for {name}")
            posre_written = True

    with open("topol.top", "w") as f:
        f.write("; GROMACS topology\n\n")
        f.write(f'#include "{os.path.basename(ff_path)}"\n')
        for itp in itp_files:
            f.write(f'#include "{itp}"\n')
        f.write('#include "posre.itp"\n')
        f.write("\n[ system ]\n; Name\nGeneric title\n\n")
        f.write("[ molecules ]\n; Compound       #mols\n")
        for name, count in mol_counts.items():
            f.write(f"{name:<16} {count}\n")

    print("\n✅ Conversion complete.")
    print(f"- {base}.gro")
    print(f"- {base}.top")
    print("- ff.itp")
    print("- topol.top")
    print("- posre.itp")
    print("- *.itp per unique molecule (system1, system2, ...) in coordinate order")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python convert_prmtop_to_gmx.py <file.prmtop> <file.inpcrd> <toppar_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
