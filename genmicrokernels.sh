#!/usr/bin/env bash

set -e

source $CODE_DIR/merge_git_repos/sfem/venv/bin/activate

python3 numerr.py > snumerr.hpp
python3 narrowphase.py > snumtol.hpp
python3 tuv.py > stuv.hpp
