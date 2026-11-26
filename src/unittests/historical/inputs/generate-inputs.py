#!/usr/bin/env python3
"""Generate the small, deterministic trajectories used by the workflow tests.

The generated files are committed; this script documents how they were made
and is not run by `make check`. Run it from any directory:

    python3 src/unittests/historical/inputs/generate-inputs.py

It writes, next to this script:
  water4.xyz             4 flexible waters, 8 frames (IR and Raman)
  water4-sfg.lammpstrj   4 waters in two slabs with dipole columns (SFG)
  water4-vel.lammpstrj   4 waters with velocities (power spectrum)

All numbers are rounded to 6 decimals and files use LF line endings with a
final newline, as the readers require.
"""

import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
NFRAMES = 8
R_OH = 0.9572           # Angstrom
HALF_ANGLE = math.radians(104.52) / 2


def rotation(a, b, c):
    """Rotation matrix for successive rotations about x, y and z."""
    ca, sa = math.cos(a), math.sin(a)
    cb, sb = math.cos(b), math.sin(b)
    cc, sc = math.cos(c), math.sin(c)
    rx = [[1, 0, 0], [0, ca, -sa], [0, sa, ca]]
    ry = [[cb, 0, sb], [0, 1, 0], [-sb, 0, cb]]
    rz = [[cc, -sc, 0], [sc, cc, 0], [0, 0, 1]]

    def mul(p, q):
        return [[sum(p[i][k] * q[k][j] for k in range(3)) for j in range(3)]
                for i in range(3)]
    return mul(rz, mul(ry, rx))


def apply(m, v):
    return [sum(m[i][k] * v[k] for k in range(3)) for i in range(3)]


def water(center, frame, index):
    """O, H, H positions of one flexible, rotating, drifting water."""
    t = frame + 0.37 * index
    r1 = R_OH + 0.020 * math.sin(1.3 * t + index)
    r2 = R_OH + 0.020 * math.sin(1.7 * t + 2 * index + 0.5)
    half = HALF_ANGLE + 0.030 * math.sin(0.9 * t + index)
    rot = rotation(0.3 * index + 0.15 * math.sin(0.5 * t),
                   0.2 * index + 0.10 * math.cos(0.7 * t),
                   0.5 * index + 0.12 * math.sin(0.4 * t + 1.0))
    shift = [center[0] + 0.05 * math.sin(0.6 * t),
             center[1] + 0.05 * math.cos(0.8 * t),
             center[2] + 0.04 * math.sin(0.5 * t + 0.3)]
    local = [[0.0, 0.0, 0.0],
             [r1 * math.sin(half), r1 * math.cos(half), 0.0],
             [-r2 * math.sin(half), r2 * math.cos(half), 0.0]]
    return [[shift[i] + apply(rot, p)[i] for i in range(3)] for p in local]


def fmt(x):
    value = round(x, 6)
    if value == 0:
        value = 0.0  # avoid "-0.000000"
    return f"{value:.6f}"


def write(name, text):
    with open(os.path.join(HERE, name), "w", newline="\n") as handle:
        handle.write(text)


# Bulk-like cluster: molecules 1-3 are hydrogen-bond neighbours (~2.9 A),
# molecule 4 sits ~5.9 A from molecule 1, inside the 2B switching region.
CLUSTER = [(10.0, 10.0, 10.0), (12.85, 10.3, 10.2),
           (10.4, 12.9, 9.6), (15.6, 11.5, 11.0)]

# Two slabs centred on z = 0: a neighbouring pair on each side.
SLABS = [(0.0, 0.0, 8.0), (2.9, 0.4, 8.3),
         (0.5, 0.2, -8.0), (-2.4, 1.5, -8.2)]


def xyz_trajectory():
    lines = []
    for frame in range(NFRAMES):
        lines.append(str(3 * len(CLUSTER)))
        lines.append(f"frame {frame}")
        for index, center in enumerate(CLUSTER):
            for label, pos in zip("OHH", water(center, frame, index)):
                lines.append(" ".join([label] + [fmt(c) for c in pos]))
    return "\n".join(lines) + "\n"


def sfg_trajectory():
    lines = []
    for frame in range(NFRAMES):
        lines += ["ITEM: TIMESTEP", str(250 * frame),
                  "ITEM: NUMBER OF ATOMS", str(3 * len(SLABS)),
                  "ITEM: BOX BOUNDS pp pp pp",
                  "-10.0 10.0", "-10.0 10.0", "-20.0 20.0",
                  "ITEM: ATOMS id type x y z v_dipx_perm v_dipy_perm "
                  "v_dipz_perm v_dipx_ind v_dipy_ind v_dipz_ind"]
        atom = 0
        for index, center in enumerate(SLABS):
            o, h1, h2 = water(center, frame, index)
            # Permanent dipole along the H-O-H bisector, ~1.85 D; a small
            # induced dipole that varies in time. Only the O row carries it.
            bis = [(h1[i] + h2[i]) / 2 - o[i] for i in range(3)]
            norm = math.sqrt(sum(b * b for b in bis))
            perm = [1.85 * b / norm for b in bis]
            ind = [0.10 * math.sin(0.8 * frame + index + k) for k in range(3)]
            for kind, pos in zip((1, 2, 2), (o, h1, h2)):
                atom += 1
                dip = perm + ind if kind == 1 else [0.0] * 6
                lines.append(" ".join([str(atom), str(kind)]
                                      + [fmt(c) for c in pos]
                                      + [fmt(d) for d in dip]))
    return "\n".join(lines) + "\n"


def velocity_trajectory():
    lines = []
    for frame in range(NFRAMES):
        lines += ["ITEM: TIMESTEP", str(250 * frame),
                  "ITEM: NUMBER OF ATOMS", str(3 * len(CLUSTER)),
                  "ITEM: BOX BOUNDS pp pp pp",
                  "0.0 20.0", "0.0 20.0", "0.0 20.0",
                  "ITEM: ATOMS id type x y z vx vy vz"]
        atom = 0
        for index, center in enumerate(CLUSTER):
            for kind, pos in zip((1, 2, 2), water(center, frame, index)):
                atom += 1
                vel = [0.002 + 0.01 * math.sin(1.1 * frame + 0.7 * atom + k)
                       for k in range(3)]
                lines.append(" ".join([str(atom), str(kind)]
                                      + [fmt(c) for c in pos]
                                      + [fmt(v) for v in vel]))
    return "\n".join(lines) + "\n"


write("water4.xyz", xyz_trajectory())
write("water4-sfg.lammpstrj", sfg_trajectory())
write("water4-vel.lammpstrj", velocity_trajectory())
