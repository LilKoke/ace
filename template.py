import argparse
from pathlib import Path


def hoge(path: Path):
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str)
    args = parser.parse_args()
    hoge(Path(args.path))
