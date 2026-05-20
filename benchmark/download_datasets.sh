#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCCD_DATA_DIR:-"${SCRIPT_DIR}/../data"}"
CACHE_DIR="${SCCD_DOWNLOAD_CACHE_DIR:-"${DATA_DIR}/.downloads"}"

is_enabled() {
    case "${1:-0}" in
        1|ON|on|On|TRUE|true|True|YES|yes|Yes) return 0 ;;
        *) return 1 ;;
    esac
}

fetch() {
    local url="$1"
    local output="$2"

    if command -v curl >/dev/null 2>&1; then
        if [[ -s "${output}" ]]; then
            set +e
            curl -L --fail --retry 3 --continue-at - --output "${output}" "${url}"
            local status=$?
            set -e

            if [[ "${status}" -eq 33 ]]; then
                echo "server does not support resumed downloads; restarting ${output}" >&2
                rm -f "${output}"
                curl -L --fail --retry 3 --output "${output}" "${url}"
            else
                return "${status}"
            fi
        else
            curl -L --fail --retry 3 --output "${output}" "${url}"
        fi
    elif command -v wget >/dev/null 2>&1; then
        wget -c -O "${output}" "${url}"
    else
        echo "error: curl or wget is required to download datasets" >&2
        return 1
    fi
}

download_archive() {
    local name="$1"
    local url="$2"
    local required_path="$3"
    local archive="${CACHE_DIR}/${name}.tar.gz"

    if [[ -e "${required_path}" ]]; then
        echo "dataset archive already present: ${required_path}"
        return 0
    fi

    mkdir -p "${CACHE_DIR}" "${DATA_DIR}"
    if [[ -f "${archive}" ]] && tar -tzf "${archive}" >/dev/null 2>&1; then
        echo "using cached dataset archive: ${archive}"
    else
        fetch "${url}" "${archive}"
    fi
    tar -xzf "${archive}" -C "${DATA_DIR}"
}

download_dataset() {
    local name="$1"
    local url="$2"

    download_archive \
        "${name}" \
        "${url}" \
        "${DATA_DIR}/${name}"
}

if is_enabled "${SCCD_ENABLE_ARMADILLO_ROLLERS:-0}"; then
    download_dataset \
        "armadillo-rollers" \
        "https://archive.nyu.edu/bitstream/2451/74508/3/armadillo-rollers.tar.gz"
fi

if is_enabled "${SCCD_ENABLE_CLOTH_BALL:-0}"; then
    download_dataset \
        "cloth-ball" \
        "https://archive.nyu.edu/bitstream/2451/74508/4/cloth-ball.tar.gz"
fi

if is_enabled "${SCCD_ENABLE_CLOTH_FUNNEL:-0}"; then
    download_dataset \
        "cloth-funnel" \
        "https://archive.nyu.edu/bitstream/2451/74508/5/cloth-funnel.tar.gz"
fi

if is_enabled "${SCCD_ENABLE_N_BODY_SIMULATION:-0}"; then
    download_dataset \
        "n-body-simulation" \
        "https://archive.nyu.edu/bitstream/2451/74508/6/n-body-simulation.tar.gz"
fi

if is_enabled "${SCCD_ENABLE_PUFFER_BALL:-0}"; then
    download_archive \
        "puffer-ball-boxes+queries+mma_bool+roots" \
        "https://archive.nyu.edu/bitstream/2451/74508/8/puffer-ball-boxes%2bqueries%2bmma_bool%2broots.tar.gz" \
        "${DATA_DIR}/puffer-ball/boxes"
    download_archive \
        "puffer-ball-frames" \
        "https://archive.nyu.edu/bitstream/2451/74508/9/puffer-ball-frames.tar.gz" \
        "${DATA_DIR}/puffer-ball/frames"
fi

if is_enabled "${SCCD_ENABLE_ROD_TWIST:-0}"; then
    download_archive \
        "rod-twist-boxes+queries+mma_bool+roots" \
        "https://archive.nyu.edu/bitstream/2451/74508/7/rod-twist-boxes%2bqueries%2bmma_bool%2broots.tar.gz" \
        "${DATA_DIR}/rod-twist/boxes"
    download_archive \
        "rod-twist-frames-0-999" \
        "https://archive.nyu.edu/bitstream/2451/74508/10/rod-twist-frames-0-999.tar.gz" \
        "${DATA_DIR}/rod-twist/frames/0.ply"
    download_archive \
        "rod-twist-frames-1000-1999" \
        "https://archive.nyu.edu/bitstream/2451/74508/11/rod-twist-frames-1000-1999.tar.gz" \
        "${DATA_DIR}/rod-twist/frames/1000.ply"
    download_archive \
        "rod-twist-frames-2000-2999" \
        "https://archive.nyu.edu/bitstream/2451/74508/12/rod-twist-frames-2000-2999.tar.gz" \
        "${DATA_DIR}/rod-twist/frames/2000.ply"
    download_archive \
        "rod-twist-frames-3000-4000" \
        "https://archive.nyu.edu/bitstream/2451/74508/13/rod-twist-frames-3000-4000.tar.gz" \
        "${DATA_DIR}/rod-twist/frames/3000.ply"
fi
