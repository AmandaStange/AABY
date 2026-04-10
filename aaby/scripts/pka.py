from propka.run import single
from pathlib import Path
import shutil
import sys

if len(sys.argv) < 2:
    print("Usage: python pka.py <pdb_file> [pH]")
    sys.exit(1)

pdb_file = Path(sys.argv[1]).resolve()

ph = 7.4
if len(sys.argv) > 2:
    try:
        ph = float(sys.argv[2])
    except ValueError:
        print("Invalid pH value. Using default pH 7.4.")

# Expected output filename produced by PROPKA
generated_pka = Path.cwd() / f"{pdb_file.stem}.pka"

mol = single(
    str(pdb_file),
    optargs=[f"--pH={ph}", "-q"],
)

# Move the generated file next to the input PDB
target_pka = pdb_file.with_suffix(".pka")

if generated_pka.exists():
    shutil.move(str(generated_pka), str(target_pka))
else:
    raise FileNotFoundError(f"Expected PROPKA output not found: {generated_pka}")