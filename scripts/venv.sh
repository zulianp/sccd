# Sourced by the analysis scripts. Not executable on its own.
#
# The analysis and plotting scripts need SymPy, NumPy, pandas and Matplotlib,
# none of which the library itself requires. Point SCCD_VENV at a virtual
# environment holding them, or leave it unset and have them on the system
# python3. python/README.md says how to create one.
#
# Nothing here fails when no environment is found: the scripts then run against
# whatever python3 is on PATH, and fail on the import if it is not enough.

for _sccd_venv in "${SCCD_VENV:-}" "${PROJECT_DIR}/.venv" "${PROJECT_DIR}/data/venv"; do
    if [ -n "${_sccd_venv}" ] && [ -f "${_sccd_venv}/bin/activate" ]; then
        # shellcheck source=/dev/null
        source "${_sccd_venv}/bin/activate"
        break
    fi
done
unset _sccd_venv
