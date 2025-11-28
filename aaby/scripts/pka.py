from propka.run import single
import sys

if len(sys.argv) < 2:
    print("Usage: python pka.py <pdb_file>")
    sys.exit(1)

pdb_file = sys.argv[1]
ph = 7.4
if len(sys.argv) > 2:
    try:
        ph = float(sys.argv[2])
    except ValueError:
        print("Invalid pH value. Using default pH 7.4.")


mol = single(
    pdb_file,
    optargs=[f"--pH={ph}", '-q'],
)