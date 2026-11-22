# MB-SPEC tests

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

After configuring and building MB-SPEC, run from `src/unittests` in the build tree:

```bash
make check
```

Fixtures are in `inputs/`; saved settings and values are in `references/`.
`generate_references.py` explicitly regenerates references; test execution never
updates them.

