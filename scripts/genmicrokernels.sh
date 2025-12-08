#!/usr/bin/env bash

source $CODE_DIR/merge_git_repos/sfem/venv/bin/activate

set -e

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR=$SCRIPTPATH/..

python3 $PROJECT_DIR/python/numerr.py 		> $PROJECT_DIR/src/snumerr.hpp
python3 $PROJECT_DIR/python/narrowphase.py 	> $PROJECT_DIR/src/snumtol.hpp
python3 $PROJECT_DIR/python/tuv.py 			> $PROJECT_DIR/src/stuv.hpp
