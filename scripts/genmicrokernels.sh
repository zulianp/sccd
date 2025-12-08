#!/usr/bin/env bash

set -e

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR=$SCRIPTPATH/..
source $PROJECT_DIR/data/venv/bin/activate

python3 $PROJECT_DIR/python/numerr.py 		> $PROJECT_DIR/src/snumerr.hpp
python3 $PROJECT_DIR/python/narrowphase.py 	> $PROJECT_DIR/src/snumtol.hpp
python3 $PROJECT_DIR/python/tuv.py 			> $PROJECT_DIR/src/stuv.hpp
