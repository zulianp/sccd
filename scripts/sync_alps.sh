#!/usr/bin/env bash
#
# Sync this working tree to the CSCS Alps cluster, for building CUDA.
#
#     scripts/sync_alps.sh [destination]
#
# The destination defaults to $SCCD_ALPS_DIR, then to
# alps:/capstor/scratch/cscs/zulianp/sccd. `alps` is an ssh alias. See
# wip/ALPS.md for building and running once the source is there.

set -euo pipefail

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR="$( cd -- "${SCRIPTPATH}/.." >/dev/null 2>&1 ; pwd -P )"

DEST="${1:-${SCCD_ALPS_DIR:-alps:/capstor/scratch/cscs/zulianp/sccd}}"

# What gets sent is what git tracks, and the filter says so directly: rsync reads
# the same .gitignore files git does, rather than carrying a hand-maintained
# exclude list that drifts from them. The four --include rules ahead of it
# restore the tracked files whose .gitignore entry is a negation, since rsync's
# per-directory merge reads `!` as a literal pattern rather than as one.
#
# Measured against this tree: 151 files, 2.3 MB, and every tracked file present.
# The hand-maintained list this replaced sent 836 MB, 770 MB of which was
# benchmark/alps -- run output that had come *from* the cluster.
#
# --delete is the point, not a flourish. Sources are globbed by CMake
# (file(GLOB_RECURSE ...)), so a file left behind from an earlier sync -- a
# header that has since been renamed, a demoted kernel -- is compiled into the
# next build, and the SCCD_PUBLIC_HEADERS guard turns a stale public header into
# a configure-time FATAL_ERROR. Without it the remote tree only ever grows.
#
# Anchor every pattern with a leading slash or a trailing one. An unanchored
# pattern matches at every level: the version of this script that this replaced
# excluded 'api', which silently dropped src/api/ -- the library's only compiled
# translation unit -- and '*git', which matches .git but would also match a
# directory named 'digit'.
#
# data/ is excluded with everything else git ignores. Fetch it on the cluster
# with benchmark/scripts/download_datasets.sh rather than pushing it over the link.
exec rsync -av --delete \
    --include '/docs/*.md' \
    --exclude '/benchmark/oracle/violations-*.csv' \
    --include '/benchmark/oracle/*.csv' \
    --include '/benchmark/assessment/*.csv' \
    --include '/spikes/research/**/*.csv' \
    --include '**/figures/README.md' \
    --filter='dir-merge,- .gitignore' \
    --exclude '/.git/' \
    "${PROJECT_DIR}/" "${DEST}/"
