from .cli import parse_args
from .workflow import run_pipeline

def main():
    args = parse_args()
    run_pipeline(args)

if __name__ == "__main__":
    main()
