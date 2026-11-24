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
the extended cases. After building the programs, use the following command
from the top-level build directory to run these suites on their own.

```bash
make -C src/unittests check
```

Fixtures are in `inputs/`. Saved settings and outputs are in `references/`,
and `test-suite.log` and `.build/runs/` in the test build directory contain
diagnostics. `make distclean` removes the generated test runs.

To regenerate references deliberately, run `python3 generate_references.py --replace`
from this source directory. For an out-of-tree build, point `--bindir` at that
build's `src/` directory. The [small workflow tests](../../tests/README.md) have
their own fixtures and reference-generation procedure.

