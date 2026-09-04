#!/usr/bin/env bash
#
# Render wip/COMPARE.md to HTML with this directory's stylesheet.
#
# The document used to sit beside this script; it moved to wip/ when the
# maintainer notes were separated from the user documentation, and this script
# kept looking for it here.

set -e

SCRIPT_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR="$( cd -- "${SCRIPT_DIR}/../../.." >/dev/null 2>&1 ; pwd -P )"

markdown_py "${PROJECT_DIR}/wip/COMPARE.md" \
    -x codehilite -x tables -x fenced_code > "${SCRIPT_DIR}/COMPARE.html"

echo "wrote ${SCRIPT_DIR}/COMPARE.html"
