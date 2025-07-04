#!/usr/bin/env python3

import subprocess
import sys

def run(cmd, return_output=False):
    print("[RUN]", cmd)
    if return_output:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        # result.stdout is the output, result.stderr is error output
        return result.stdout
    else:
        subprocess.run(cmd, shell=True, check=True)

def make_ndx(input_gro):
    run(f'echo q |  gmx make_ndx -f {input_gro}.gro -o tmp.ndx > /dev/null')
    groups = run(f'grep "\[" tmp.ndx  | wc -l', return_output=True)
    print(groups)
    groups = int(groups.strip()) - 1
    print(groups)
    run(f"printf '1 \n name {groups + 1} SOLU \n \"PA\" | \"PC\" | \"OL\" \n name {groups + 2} MEMB \n {groups + 1} | {groups + 2} \n ! {groups + 3} \n name {groups + 4} SOLV \n q \n' | gmx make_ndx -f {input_gro}.gro -o index.ndx")

    
def main(base="system"):
    input_gro = f"{base}"
    final_gro = f"{base}"
    # 1. Make index file
    make_ndx(input_gro)
    print(f"[DONE] See {final_gro} for output.")

if __name__ == "__main__":
    main(sys.argv[1])
