import sys
import argparse

alt = {
    'ASP': 'ASH',
    'GLU': 'GLH',
    'LYS': 'LYN',
    'HIS': 'HIP'
}

def parse_pka(pka_file, ph):
    acid = {}
    base = {}
    his = {}
    with open(pka_file) as f:
        lines = f.readlines()
    c1 = False
    for line in lines:
        if c1:
            if line[0] == '-':
                break
            l = line.split()
            if l[0] in ['ASP', 'GLU'] and float(l[3]) - ph > 0:
                acid[f'{l[0]}-{l[1]}-{l[2]}'] = float(l[3]) - ph
            if l[0] in ['LYS'] and float(l[3]) - ph < 0:
                base[f'{l[0]}-{l[1]}-{l[2]}'] = float(l[3]) - ph
            if l[0] == 'HIS' and float(l[3]) - ph > 0:
                his[f'{l[0]}-{l[1]}-{l[2]}'] = float(l[3]) - ph
        if line == '       Group      pKa  model-pKa   ligand atom-type\n':
            c1 = True
    selected = []
    for k in list(acid.keys()) + list(base.keys()) + list(his.keys()):
        l = k.split('-')
        resname, resid, chain = l[0], l[1], l[2]
        selected.append((resname, chain, resid))
    return selected

def parse_list(listfile):
    selected = []
    with open(listfile) as f:
        for line in f:
            if not line.strip() or line.strip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 3:
                print(f"Warning: skipping malformed line: {line.strip()}", file=sys.stderr)
                continue
            resname, chain, resid = parts
            selected.append((resname, chain, resid))
    return selected

def write_sed_cmds(selected, outfile):
    with open(outfile, "w") as f:
        f.write("#!/bin/bash\n\n")
        for resname, chain, resid in selected:
            if resname not in alt:
                print(f"Warning: {resname} not in supported list, skipping", file=sys.stderr)
                continue
            f.write(f'sed -i "s/{resname} {chain}{int(resid):>4}/{alt[resname]} {chain}{int(resid):>4}/g" $1\n')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pka_or_list", help="pKa file or residue list file")
    # parser.add_argument("pdb_file", help="Target PDB file (for sed script)")
    parser.add_argument("--ph", type=float, default=7.4, help="pH value (default: 7.4)")
    parser.add_argument("--list", action="store_true", help="Treat input as residue list, not pKa file")
    parser.add_argument("-o", "--outfile", default="mutate_residues.sh", help="Output bash file (default: mutate_residues.sh)")
    args = parser.parse_args()
    if args.list:
        selected = parse_list(args.pka_or_list)
    else:
        selected = parse_pka(args.pka_or_list, args.ph)
    write_sed_cmds(selected, args.outfile)
    print(f"Wrote sed commands to {args.outfile}")

if __name__ == "__main__":
    main()
