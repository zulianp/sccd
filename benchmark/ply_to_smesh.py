#!/usr/bin/env python3
"""Convert a binary PLY triangle surface into the raw-array folder smesh reads.

smesh's mesh reader takes a directory, not a file: coordinates as ``x.float64``,
``y.float64``, ``z.float64`` and connectivity as ``i0.int32`` .. ``iN.int32``,
one file per element corner. The CCD benchmark frames under ``data/<scene>/frames``
are binary PLY, so they need converting before ``sccd_refine_scaling`` can use
them.

Two frames of the same scene share a connectivity, which is what lets the
scaling benchmark refine both consistently: it refines each mesh independently
and relies on refinement being a deterministic function of the topology.

    ply_to_smesh.py <in.ply> <out_dir>
"""

import struct
import sys
from pathlib import Path


# PLY spells the same types several ways; map them all to struct codes.
_SCALAR = {
    "float": "f", "float32": "f", "double": "d", "float64": "d",
    "char": "b", "int8": "b", "uchar": "B", "uint8": "B",
    "short": "h", "int16": "h", "ushort": "H", "uint16": "H",
    "int": "i", "int32": "i", "uint": "I", "uint32": "I",
}


def read_ply(path):
    """Return (n_vertices, coords[3][n], faces[3][m]).

    Handles ascii, binary_little_endian and binary_big_endian, because the CCD
    benchmark scenes are not consistent: cloth-funnel is little-endian binary and
    cloth-ball is big-endian binary.
    """
    with open(path, "rb") as handle:
        header = []
        while True:
            line = handle.readline()
            if not line:
                raise ValueError(f"{path}: truncated header")
            header.append(line.strip().decode("ascii", "replace"))
            if header[-1] == "end_header":
                break

        fmt_line = next((h for h in header if h.startswith("format")), "")
        if "binary_little_endian" in fmt_line:
            byte_order = "<"
        elif "binary_big_endian" in fmt_line:
            byte_order = ">"
        elif "ascii" in fmt_line:
            byte_order = None
        else:
            raise ValueError(f"{path}: unrecognised format line {fmt_line!r}")

        n_vertices = n_faces = 0
        vertex_props = []
        face_list = None
        element = None
        for line in header:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "element":
                element = parts[1]
                if element == "vertex":
                    n_vertices = int(parts[2])
                elif element == "face":
                    n_faces = int(parts[2])
            elif parts[0] == "property" and element == "vertex" and parts[1] != "list":
                vertex_props.append(parts[1])
            elif parts[0] == "property" and element == "face" and parts[1] == "list":
                face_list = (parts[2], parts[3])

        coords = [[0.0] * n_vertices for _ in range(3)]
        faces = [[], [], []]

        def fan(idx):
            for k in range(1, len(idx) - 1):
                faces[0].append(idx[0])
                faces[1].append(idx[k])
                faces[2].append(idx[k + 1])

        if byte_order is None:
            tokens = handle.read().split()
            cursor = 0
            stride = len(vertex_props)
            for i in range(n_vertices):
                coords[0][i] = float(tokens[cursor])
                coords[1][i] = float(tokens[cursor + 1])
                coords[2][i] = float(tokens[cursor + 2])
                cursor += stride
            for _ in range(n_faces):
                count = int(tokens[cursor])
                cursor += 1
                fan([int(t) for t in tokens[cursor:cursor + count]])
                cursor += count
        else:
            try:
                codes = "".join(_SCALAR[p] for p in vertex_props)
            except KeyError as exc:
                raise ValueError(f"{path}: unsupported vertex property {exc}") from exc

            reader = struct.Struct(byte_order + codes)
            blob = handle.read(reader.size * n_vertices)
            for i in range(n_vertices):
                values = reader.unpack_from(blob, i * reader.size)
                coords[0][i], coords[1][i], coords[2][i] = values[0], values[1], values[2]

            count_code = _SCALAR[face_list[0]] if face_list else "B"
            index_code = _SCALAR[face_list[1]] if face_list else "i"
            count_size = struct.calcsize(byte_order + count_code)
            index_size = struct.calcsize(byte_order + index_code)
            rest = handle.read()
            offset = 0
            for _ in range(n_faces):
                (count,) = struct.unpack_from(byte_order + count_code, rest, offset)
                offset += count_size
                idx = struct.unpack_from(byte_order + index_code * count, rest, offset)
                offset += index_size * count
                fan(idx)

    return n_vertices, coords, faces


def write_smesh_folder(out_dir, coords, faces):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, values in zip(("x", "y", "z"), coords):
        with open(out / f"{name}.float64", "wb") as handle:
            handle.write(struct.pack(f"<{len(values)}d", *values))
    for k, values in enumerate(faces):
        with open(out / f"i{k}.int32", "wb") as handle:
            handle.write(struct.pack(f"<{len(values)}i", *values))


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 1
    n_vertices, coords, faces = read_ply(sys.argv[1])
    write_smesh_folder(sys.argv[2], coords, faces)
    print(f"{sys.argv[1]}: {n_vertices} vertices, {len(faces[0])} triangles -> {sys.argv[2]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
