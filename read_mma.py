# TODO read the mma_bool files and return a list of boolean values following the README.md file

import json
from typing import List

def read_mma_bool(path: str) -> List[bool]:
    with open(path, 'r') as file:
        data = json.load(file)
    return data

if __name__ == "__main__":
    path = "data/armadillo-rollers/mma_bool/5vf_mma_bool.json"
    print(read_mma_bool(path))


if __name__ == "__main__":
    import sys
    print(read_mma_bool(sys.argv[1]))