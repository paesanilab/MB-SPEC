# Trajectory regression suites

The suites validate IR dipoles and absorption spectra, Raman polarizabilities and
spectra, collective-velocity power spectra, and SFG responses, complex spectra,
and Fresnel-corrected intensities.

| Suite | Numerical references | Analytical | Other | Total |
|---|---:|---:|---:|---:|
| IR | 24 | 0 | 13 | **37** |
| Raman | 18 | 0 | 12 | **30** |
| POWER | 10 | 0 | 5 | **15** |
| SFG | 48 | 6 | 10 | **64** |
| **Total** | **100** | **6** | **40** | **146** |

Numerical checks compare saved reference tables. Analytical checks evaluate SFG
interface weighting and neighbor exclusion. Other checks verify equivalent input
representations, atom exclusion, switching, normalization, and thread agreement.
A named check may contain multiple calculations.

`make check` from the top-level build directory runs these suites, including
the extended cases and the historical workflow check described below. After
building the programs, run these suites on their own from the build directory.

```bash
make -C src/unittests check
```

Fixtures are in `inputs/`. Saved settings and outputs are in `references/`,
and `test-suite.log` and `.build/runs/` in the test build directory contain
diagnostics. `make distclean` removes the generated test runs.

To regenerate references deliberately, run `python3 generate_references.py --replace`
from this source directory. For an out-of-tree build, point `--bindir` at that
build's `src/` directory.

## Historical workflow check

`check-historical` adds four end-to-end cases with eight-frame fixtures and 20
saved outputs from the `master` commit recorded in [`historical/references/PROVENANCE`](historical/references/PROVENANCE).

| Workflow | Coverage |
|---|---|
| IR | Internal `1b+2b` dipoles, mean-subtracted response, and spectrum |
| Raman | Full `1b+2b+nb` polarizabilities and both response modes with Bose–Einstein correction |
| SFG | Three components with `1b+2b` dipoles and polarizabilities, both correlation modes, and Fresnel factors with nonzero redshift and background |
| Power | LAMMPS velocities, response, and spectrum |

The check runs with one thread. It compares finite numeric tables with a
relative tolerance of `1e-8`, an absolute tolerance of `1e-10` times the
reference column scale, and one unit in the last printed digit of the
reference. The IR response uses `1e-8` for the column-scaled term because
rounding in the printed dipoles propagates through mean subtraction and
autocorrelation. Generate the compact trajectories with
`python3 historical/inputs/generate-inputs.py` from this directory.

CI also builds the pinned commit independently with Clang and Intel MKL and
runs the same check against the old executables. After configuring the current
checkout and initializing MKL (`MKLROOT`), run these commands from the repository
root to reproduce that verification locally. The historical build requires GSL.

```bash
git worktree add --detach /tmp/mbspec-master c516736c3946166e0ee81b2198854a67da660c99
sh src/unittests/build-historical.sh /tmp/mbspec-master
make -C src/unittests check-historical
src/unittests/check-historical /tmp/mbspec-master/src
```

Use Clang for the historical build. The original C++11 LAMMPS parser depends
on operand evaluation order when assigning column indices. Clang evaluates
them in the order used to generate the references. The current code's C++17
build specifies that order.

