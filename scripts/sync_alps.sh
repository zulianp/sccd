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

# What gets sent is what git tracks: rsync reads the same .gitignore files git
# does, so there is no exclude list here to keep in step with them. The four
# --include rules ahead of the merge restore the tracked files whose .gitignore
# entry is a negation, since rsync's per-directory merge reads `!` as a literal
# pattern rather than as one. That comes to 151 files and 2.3 MB.
#
# --delete matters. Sources are globbed by CMake (file(GLOB_RECURSE ...)), so a
# file left behind from an earlier sync -- a renamed header, a demoted kernel --
# is compiled into the next build, and the SCCD_PUBLIC_HEADERS guard turns a
# stale public header into a configure-time FATAL_ERROR. Excluded paths are not
# deleted, so the remote build trees survive a re-sync.
#
# Anchor every pattern with a leading or trailing slash. An unanchored pattern
# matches at every level: `api` would drop src/api/, the library's only compiled
# translation unit, and `*git` matches .git but also a directory named `digit`.
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
