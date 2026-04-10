#!/usr/bin/env python3

import argparse
import os
import re
import shlex
import subprocess
import sys
from collections import defaultdict, deque
from pathlib import Path


LOGFILE = None


def run(cmd, check=True, shell=None, capture_output=True, input_text=None):
    """Run a shell command or Python script with optional logging."""
    global LOGFILE

    if isinstance(cmd, str):
        if shell is None:
            if re.search(r"[><|;&]", cmd):
                shell = True
            else:
                shell = False
        cmd_list = cmd if shell else shlex.split(cmd)
    else:
        cmd_list = list(map(str, cmd))
        if shell is None:
            shell = False

    cmd_str = cmd if isinstance(cmd, str) else " ".join(map(str, cmd_list))
    print(f"Running: {cmd_str}")

    result = subprocess.run(
        cmd_list,
        check=check,
        shell=shell,
        text=True,
        capture_output=capture_output,
        input=input_text,
    )

    if capture_output:
        if result.stdout:
            print(result.stdout, end="")
        if result.stderr:
            print(result.stderr, end="", file=sys.stderr)

    if LOGFILE:
        with open(LOGFILE, "a") as f:
            f.write(f"\n$ {cmd_str}\n")
            if input_text:
                f.write("[stdin]\n")
                f.write(input_text)
                if not input_text.endswith("\n"):
                    f.write("\n")
            if capture_output and result.stdout:
                f.write(result.stdout)
            if capture_output and result.stderr:
                f.write(result.stderr)

    return result


SECTION_RE = re.compile(r"^\s*\[\s*(\w+)\s*\]")


def read_topology_sections(top_file):
    sections = defaultdict(list)
    current = None

    with open(top_file) as fh:
        for line in fh:
            m = SECTION_RE.match(line)
            if m:
                current = m.group(1).lower()
                sections[current].append(line)
                continue
            if current is not None:
                sections[current].append(line)

    return sections


def parse_atoms_section(atom_lines):
    """
    Parse [ atoms ] into:
      atom_to_resnr: topology atom index -> resnr
    """
    atom_to_resnr = {}
    for line in atom_lines:
        s = line.strip()
        if not s or s.startswith(";") or s.startswith("["):
            continue
        parts = s.split()
        if len(parts) < 3:
            continue
        if not parts[0].isdigit():
            continue

        atomnr = int(parts[0])
        resnr = int(parts[2])
        atom_to_resnr[atomnr] = resnr

    return atom_to_resnr


def parse_bonds_section(bond_lines):
    bonds = []
    for line in bond_lines:
        s = line.strip()
        if not s or s.startswith(";") or s.startswith("["):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        if not (parts[0].isdigit() and parts[1].isdigit()):
            continue

        ai = int(parts[0])
        aj = int(parts[1])
        bonds.append((ai, aj))

    return bonds


def read_pdb_atom_chain_order(pdb_file):
    """
    Read ATOM/HETATM records from PDB in file order.
    Returns:
      atom_to_chain: topology-like atom index (1-based by order) -> chain ID
      chain_order: order of first appearance of chains in file
    """
    atom_to_chain = {}
    chain_order = []
    seen_chains = set()

    idx = 0
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                idx += 1
                chain = line[21].strip() or "-"
                atom_to_chain[idx] = chain
                if chain not in seen_chains:
                    seen_chains.add(chain)
                    chain_order.append(chain)

    return atom_to_chain, chain_order


def connected_components(atom_numbers, bonds):
    adj = defaultdict(set)
    for ai, aj in bonds:
        adj[ai].add(aj)
        adj[aj].add(ai)

    visited = set()
    comps = []

    for atom in sorted(atom_numbers):
        if atom in visited:
            continue

        q = deque([atom])
        visited.add(atom)
        comp = []

        while q:
            u = q.popleft()
            comp.append(u)
            for v in adj[u]:
                if v not in visited:
                    visited.add(v)
                    q.append(v)

        comps.append(sorted(comp))

    return comps


def summarise_components(components, atom_to_chain, atom_to_resnr):
    summaries = []

    for i, comp in enumerate(components, start=1):
        chains = []
        seen = set()
        for a in comp:
            c = atom_to_chain.get(a, "?")
            if c not in seen:
                seen.add(c)
                chains.append(c)

        resnrs = sorted({atom_to_resnr[a] for a in comp if a in atom_to_resnr})

        summaries.append({
            "index": i,
            "natoms": len(comp),
            "chains": chains,
            "atom_min": min(comp),
            "atom_max": max(comp),
            "res_min": min(resnrs) if resnrs else None,
            "res_max": max(resnrs) if resnrs else None,
            "atoms": comp,
        })

    return summaries


def union_find_groups_from_components(component_summaries, chain_order):
    """
    Convert component summaries into chain groups in original chain order.
    """
    parent = {c: c for c in chain_order}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for comp in component_summaries:
        chains = [c for c in comp["chains"] if c in parent]
        if len(chains) >= 2:
            root = chains[0]
            for c in chains[1:]:
                union(root, c)

    groups = defaultdict(list)
    for c in chain_order:
        groups[find(c)].append(c)

    ordered_groups = []
    seen = set()
    for c in chain_order:
        root = find(c)
        if root not in seen:
            seen.add(root)
            ordered_groups.append(groups[root])

    return ordered_groups


def merge_answers_from_groups(chain_order, chain_groups):
    """
    For chain order A B C D E F and groups [ABCD, EF],
    generate answers y y y n y
    """
    chain_to_group = {}
    for gi, group in enumerate(chain_groups):
        for c in group:
            chain_to_group[c] = gi

    answers = []
    for i in range(len(chain_order) - 1):
        c1 = chain_order[i]
        c2 = chain_order[i + 1]
        answers.append("y" if chain_to_group[c1] == chain_to_group[c2] else "n")

    return answers

def build_pdb2gmx_cmd(
    gmx,
    input_pdb,
    output_pdb,
    output_top,
    ff,
    water,
    merge_mode,
    ignh=False,
    chainsep=None,
    extra_args=None,
):
    cmd = [
        gmx, "pdb2gmx",
        "-f", str(input_pdb),
        "-o", str(output_pdb),
        "-p", str(output_top),
        "-ff", ff,
        "-water", water,
        "-merge", merge_mode,
    ]

    if ignh:
        cmd.append("-ignh")
    if chainsep:
        cmd.extend(["-chainsep", chainsep])
    if extra_args:
        cmd.extend(extra_args)

    return cmd


def main():
    parser = argparse.ArgumentParser(
        description="Run pdb2gmx twice: first with merge all to detect connectivity, then rerun with merge interactive."
    )
    parser.add_argument("-f", "--input", required=True, help="Input PDB file")
    parser.add_argument("-ff", "--forcefield", default='amber19sb_gromacs_lipid21', help="Force field name for pdb2gmx")
    parser.add_argument("-water", "--water", default='opc', help="Water model for pdb2gmx")
    parser.add_argument("--gmx", default="gmx", help="GROMACS executable (default: gmx)")
    parser.add_argument("--ignh", action="store_true", help="Pass -ignh to pdb2gmx")
    parser.add_argument("--chainsep", default=None, help="Pass -chainsep value to pdb2gmx")
    parser.add_argument("--log", default=None, help="Optional logfile")
    parser.add_argument("--prefix", default="pdb2gmx", help="Prefix for intermediate/final output files")
    parser.add_argument(
        "--extra",
        nargs=argparse.REMAINDER,
        help="Extra arguments passed through to pdb2gmx after '--extra'",
    )

    args = parser.parse_args()

    global LOGFILE
    LOGFILE = args.log

    input_pdb = Path(args.input).resolve()
    prefix = args.prefix

    pass1_pdb = Path(f"{prefix}_pass1_allmerge.pdb")
    pass1_top = Path(f"{prefix}_pass1_allmerge.top")

    final_pdb = Path(f"{prefix}.pdb")
    final_top = Path(f"{prefix}.top")

    extra_args = args.extra if args.extra else []

    # -------------------------
    # First pass: merge all
    # -------------------------
    cmd1 = build_pdb2gmx_cmd(
        gmx=args.gmx,
        input_pdb=input_pdb,
        output_pdb=pass1_pdb,
        output_top=pass1_top,
        ff=args.forcefield,
        water=args.water,
        merge_mode="all",
        ignh=args.ignh,
        chainsep=args.chainsep,
        extra_args=extra_args,
    )
    run(cmd1, check=True, capture_output=True)

    # -------------------------
    # Parse connectivity
    # -------------------------
    atom_to_chain, chain_order = read_pdb_atom_chain_order(pass1_pdb)
    sections = read_topology_sections(pass1_top)

    if "atoms" not in sections:
        raise RuntimeError(f"No [ atoms ] section found in {pass1_top}")
    if "bonds" not in sections:
        raise RuntimeError(f"No [ bonds ] section found in {pass1_top}")

    atom_to_resnr = parse_atoms_section(sections["atoms"])
    bonds = parse_bonds_section(sections["bonds"])
    topo_atoms = sorted(atom_to_resnr.keys())

    if len(atom_to_chain) != len(topo_atoms):
        print(
            f"Warning: {pass1_pdb} has {len(atom_to_chain)} ATOM/HETATM records but "
            f"{pass1_top} [ atoms ] has {len(topo_atoms)} atoms.\n"
            "Connectivity mapping assumes atom order is identical.",
            file=sys.stderr,
        )

    components = connected_components(topo_atoms, bonds)
    summaries = summarise_components(components, atom_to_chain, atom_to_resnr)
    chain_groups = union_find_groups_from_components(summaries, chain_order)

    print("\nDetected connected molecule(s):\n")
    for i, group in enumerate(chain_groups, start=1):
        matching = [s for s in summaries if set(s["chains"]) == set(group)]
        if matching:
            s = matching[0]
            print(f"Molecule {i}:")
            print(f"  chains:   {', '.join(group)}")
            print(f"  atoms:    {s['atom_min']} - {s['atom_max']} ({s['natoms']} atoms)")
            if s["res_min"] is not None:
                print(f"  residues: {s['res_min']} - {s['res_max']}")
            print()
        else:
            print(f"Molecule {i}: chains {', '.join(group)}\n")

    print("Chain groups:", " | ".join("".join(g) for g in chain_groups))

    # -------------------------
    # Generate merge answers
    # -------------------------
    merge_answers = merge_answers_from_groups(chain_order, chain_groups)
    merge_input = "\n".join(merge_answers) + "\n"

    print("\nChain order seen by pdb2gmx:")
    print("  " + " ".join(chain_order))
    print("Automatic merge answers:")
    print("  " + " ".join(merge_answers))
    print()

    # -------------------------
    # Second pass: interactive merge
    # -------------------------
    cmd2 = build_pdb2gmx_cmd(
        gmx=args.gmx,
        input_pdb=input_pdb,
        output_pdb=final_pdb,
        output_top=final_top,
        ff=args.forcefield,
        water=args.water,
        merge_mode="interactive",
        ignh=args.ignh,
        chainsep=args.chainsep,
        extra_args=extra_args,
    )
    run(cmd2, check=True, capture_output=True, input_text=merge_input)

    print("\nDone.")
    print(f"Final structure: {final_pdb}")
    print(f"Final topology:  {final_top}")


if __name__ == "__main__":
    main()