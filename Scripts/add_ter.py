import sys

def add_ter_on_chain_change(infile, outfile):
    with open(infile) as fin:
        lines = fin.readlines()

    out_lines = []
    last_chain = None
    last_resnum = None

    for i, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")):
            chain = line[21]
            resnum = line[22:26]
            if last_chain is None:
                last_chain = chain
            elif chain != last_chain:
                # Insert TER before this line
                ter_line = f"TER   {line[6:11]}      {last_resname} {last_chain}{last_resnum}\n"
                out_lines.append(ter_line)
                last_chain = chain
            last_resname = line[17:20]
            last_resnum = resnum
        out_lines.append(line)

    # Add TER at the end if the last line was an atom
    if lines and lines[-1].startswith(("ATOM", "HETATM")):
        ter_line = f"TER   {lines[-1][6:11]}      {last_resname} {last_chain}{last_resnum}\n"
        out_lines.append(ter_line)

    with open(outfile, "w") as fout:
        fout.writelines(out_lines)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python add_ter_on_chain_change.py input.pdb output.pdb")
        sys.exit(1)
    add_ter_on_chain_change(sys.argv[1], sys.argv[2])

