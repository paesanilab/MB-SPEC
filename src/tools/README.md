# Shared library and water property models

`src/tools` provides the internal `libtools.a` library used by the IR, Raman,
SFG, and power spectrum programs. It contains trajectory readers, configuration
handling, geometry and numerical utilities, and the water dipole and
polarizability models. It does not provide a standalone command-line program.

## File overview

In the tables below, `.cpp/.h` denotes an implementation and its public interface.

### Configuration, trajectories, and data structures

| Files | Role |
|---|---|
| `system.cpp/.h` | JSON configuration, defaults, option accessors, and simulation-cell conversion. |
| `read-lammpstrj.cpp/.h` | LAMMPS and XYZ frame readers; coordinate and dipole fields are identified from LAMMPS column headers. |
| `struct-lammpstrj-frame.h` | Frame storage for coordinates, atom information, dipoles, cell matrices, and times, including allocation and copying. |
| `struct-sfg-pols.h` | Molecular property components, interface weights, and correlation partner lists used by SFG. |
| `io-tools.cpp/.h` | Input-error reporting, inverse-cell construction, and atom-label utilities. |
| `json/json.h` | Bundled nlohmann JSON library, version 3.7.0; its MIT license notice is included in the header. |

### Geometry and numerical utilities

| Files | Role |
|---|---|
| `pbc-tools.cpp/.h` | Periodic imaging, minimum-image distances, monomer imaging, and slab centering. |
| `orient.cpp/.h` | Vector operations, angles, rotations, Euler transformations, centering, and inertia calculations. |
| `mat_tools.cpp/.h` | Matrix inversion, eigenvalues/eigenvectors, multiplication, norms, and convergence checks. Inversion and diagonalization use LAPACKE. |
| `spectral_smoothing.h` | Gaussian spectral smoothing with a three-standard-deviation kernel and local boundary renormalization. |
| `constants.h` | Physical constants, atomic masses, and unit conversions used by library routines. |
| `macros.h` | Shared compiler and atom-label macros. |

### Water dipole and polarizability surfaces

| Files | Role |
|---|---|
| `ltp2011.cpp/.h` | LTP2011 one-body water dipole surface and laboratory-frame dipole evaluation. |
| `ps.cpp/.h` | Partridge–Schwenke (PS) water monomer potential and dipole routines; supplies geometry-dependent charges to the TTM4 models and optional permanent dipoles to SFG. |
| `dip_poly2b5.cpp/.h` | Fitted coefficients and evaluation of the explicit two-body dipole polynomial. |
| `poly2b5.cpp/.h`, `poly2b5-{ha1,ha2,hb1,hb2,oa,ob}.cpp` | Polynomial basis evaluation used by `dip_poly2b5`, with a wrapper and six separately compiled output blocks. |
| `avila.cpp/.h` | Avila one-body water polarizability and transformation into the laboratory frame. |
| `alpha-el-poly2b3.cpp/.h` | Fitted coefficients and evaluation of the explicit two-body polarizability correction. |
| `poly-l-2b3.cpp/.h` | Polynomial basis evaluation used by `alpha-el-poly2b3`. |
| `lp-site.cpp/.h` | Construction of two auxiliary lone-pair sites from an O–H–H monomer geometry. |

### TTM4 electrostatics and response

| Files | Role |
|---|---|
| `ttm4-dip.cpp/.h` | Finite-system induced-dipole evaluation without energy or gradient output; used for the isolated-dimer term replaced in the full IR and SFG dipole decompositions. |
| `ttm4-es.cpp/.h` | Finite-system electrostatic and induced energies, gradients, dipoles, and polarizability response. |
| `ttm4-periodic-es.cpp/.h` | Cell-dependent TTM4 dipole and polarizability response; used for Raman nonadditive tensors. |
| `ttm4-smear.cpp/.h` | Damped electrostatic interaction functions used by the TTM4 models. |
| `gammq.cpp/.h` | Gamma-function routines used by the electrostatic damping functions. |
| `alpha-ttm4.cpp/.h` | Wrappers for total finite-system TTM4 polarizability and subtraction of isolated-monomer tensors. |

## Model dependencies and use

The main model paths are:

- **IR:** `ltp2011` supplies the internal one-body dipole. `dip_poly2b5` calls
  `poly2b5` for explicit two-body dipoles. Full external decompositions use
  `ttm4-dip` to replace the corresponding short-range induced-dimer term.
- **Raman:** `avila` supplies the one-body tensor. `alpha-el-poly2b3` uses
  `lp-site` and `poly-l-2b3` for two-body corrections. NB-containing modes use
  `ttm4-periodic-es` for the total-system tensor. The caller subtracts the
  isolated-monomer tensors evaluated with `alpha_tot` in `alpha-ttm4`.
  The separate `alpha_no_1b` wrapper performs finite-system monomer subtraction
  but is not called by the current spectroscopy programs.
- **SFG:** selected dipole components are combined with `avila` tensors and
  optional explicit two-body corrections. Permanent dipoles come from the
  trajectory by default; `use_ps_dipole=true` replaces them with internal PS
  monomer dipoles. SFG currently has no internal LTP2011 option. Properties
  and partner lists are stored in
  `struct-sfg-pols.h`.
- **POWER:** shared configuration and smoothing are used. The native LAMMPS
  velocity reader is implemented in `src/power/response.cpp`, rather than in
  the shared coordinate/dipole reader.

The TTM4 models use `ttm4-smear`, which in turn uses `gammq`. TTM4 model classes
retain mutable work storage; use a separate instance for each concurrent
worker. Dipole and polarizability response accessors require the associated
model evaluation to have initialized that storage.
The `dipind` accessors sum induced dipoles for the supplied system; their
output is a dimer contribution only when the evaluated system is a dimer.

## Units and calling conventions

- Water model coordinates are Cartesian Å, ordered as contiguous O–H–H
  triplets. A monomer occupies nine values; a dimer is passed as two monomer
  coordinate arrays. Callers must provide consistently imaged geometries.
- Dipole vectors use the order `x, y, z` and are interpreted in Debye by the
  spectroscopy workflows. Polarizability tensors use Å³. Electrostatic
  routines also use internal scaled quantities; consult their interfaces and
  conversion code before using energies or gradients directly.
- Nine-component polarizability tensors use row-major order:
  `xx, xy, xz, yx, yy, yz, zx, zy, zz`.
  `avila::alpha6` instead returns `xx, yy, zz, xy, xz, yz`. The Raman property
  output uses `xx, xy, xz, yy, yz, zz`; these six-component layouts differ.
- `dip_poly2b5::eval` and `alpha_el_poly2b3::eval` return polynomial values
  without the distance switch. Spectroscopy callers apply switching,
  periodic imaging, pair summation, and any induced-dipole replacement.
- Frame times are assigned by callers in ps. LAMMPS coordinates are found by
  header, with one complete triplet of `x y z`, `xu yu zu`, `xs ys zs`, or
  `xsu ysu zsu`. Scaled coordinates are converted with the cell and origin.
  When dipoles are required, the shared reader expects the permanent and
  induced `v_dip{x,y,z}_perm` and `v_dip{x,y,z}_ind` components.
- The shared XYZ reader reads atom labels and three numeric components;
  their meaning depends on the calling workflow. Coordinate workflows use
  positions, while POWER supplies velocities. The cell is supplied separately.

## Scientific and numerical limitations

The built-in property surfaces are water-specific. External dipole definitions
are determined by the trajectory producer; PS, LTP2011, and externally supplied
permanent dipoles must not be treated as interchangeable without justification.
The published MB-μ construction uses LTP2011 for its explicit one-body term
[4]. The optional PS permanent term in SFG is a different model choice;
PS-derived charges inside TTM4 induction serve a separate role.

The internal cell-dependent TTM4 treatment uses only the diagonal cell lengths
for minimum imaging and therefore assumes an orthorhombic cell, including
cubic cells. Triclinic support in the shared trajectory reader does not extend
to this model. It has no Ewald summation. Raman
NB-containing modes are not validated for periodic bulk calculations. The word
`periodic` in the filename does not establish a complete periodic electrostatic
treatment.

Gaussian smoothing excludes zero-frequency input and sets the zero-frequency
output to zero. This changes the POWER zero-frequency value when smoothing is
enabled.

## Build and generated files

`Makefile.am` defines `libtools.a`, its sources, and additional distribution
files. The library is built as part of the repository build. Matrix utilities
require LP64 LAPACKE with 32-bit integers. The default `--with-lapack=generic`
uses `lapacke.h` and tries `-llapacke`, then `-llapacke -llapack -lblas`.
`--with-lapack=mkl` selects `mkl_lapacke.h` and uses supplied flags or the
compiler environment, without guessing MKL paths, libraries, or threading
runtimes. Configuration checks compilation and linkage of `LAPACKE_dgetrf`,
`LAPACKE_dgetri`, and `LAPACKE_dsyev` and stops if they are unavailable.

Set `LAPACKE_CFLAGS` and `LAPACKE_LIBS` for a custom installation or library
combination. Explicit `LAPACKE_LIBS` are used unchanged, without fallback.
OpenBLAS must include LAPACKE support or be paired with a separate LAPACKE
library. Headers and libraries must use the same LP64 interface; ILP64 is
unsupported. The old `MKL_CFLAGS` and `MKL_LIBS` variables are rejected:
rename them to the LAPACKE variables and select `--with-lapack=mkl` for MKL.
The legacy `--enable-mkl`, `--disable-mkl`, and `--with-mkl-*` options remain
removed. See the root README for generic, custom-library, and MKL examples.
Build configuration supplies the selected LAPACKE flags to the spectroscopy
programs using this library. `make check` runs numerical matrix tests against
the selected implementation; the test executable is not installed.

`Makefile.in`, `Makefile`, `.deps/`, `*.o`, and `libtools.a` are generated build
files.

The `poly2b5.cpp` wrapper evaluates the six output blocks in their original
order. Their explicit polynomial assignments are kept in
`poly2b5-{ha1,ha2,hb1,hb2,oa,ob}.cpp`, with private declarations in
`poly2b5-internal.h`. Separate translation units reduce optimization work per
compiler invocation and allow `make -j4` to compile the blocks concurrently,
without lowering the optimization level or changing the polynomial expressions.
The public `cpt_poly2b5` interface is unchanged.

`make check` also runs polynomial regression tests against reference values
captured from the original implementation. Binary-fraction inputs allow exact
comparison of all 9,408 outputs. The tests cover zero inputs, two asymmetric
input pairs, output bounds, and overlapping output arrays.

The large `poly-l-2b3.cpp` file also contains explicit polynomial basis
expressions. Fitted coefficients remain in `dip_poly2b5.cpp` and
`alpha-el-poly2b3.cpp`.

## Scientific references

1. H. Partridge and D. W. Schwenke, *The determination of an accurate isotope
   dependent potential energy surface for water from extensive ab initio
   calculations and experimental data*, J. Chem. Phys. **106**, 4618–4639
   (1997). [Paper](https://history.arc.nasa.gov/hist_pdfs/awards/hjaa_1999.pdf).
2. L. Lodi, J. Tennyson, and O. L. Polyansky, *A global, high accuracy ab initio
   dipole moment surface for the electronic ground state of the water
   molecule*, J. Chem. Phys. **135**, 034113 (2011).
   [DOI: 10.1063/1.3604934](https://doi.org/10.1063/1.3604934).
3. G. Avila, *Ab initio dipole polarizability surfaces of water molecule:
   Static and dynamic at 514.5 nm*, J. Chem. Phys. **122**, 144310 (2005).
   [DOI: 10.1063/1.1867437](https://doi.org/10.1063/1.1867437).
4. G. R. Medders and F. Paesani, *Infrared and Raman Spectroscopy of Liquid
   Water through “First-Principles” Many-Body Molecular Dynamics*, J. Chem.
   Theory Comput. **11**, 1145–1154 (2015). Defines the MB-μ and MB-α
   constructions, including their TTM4-Ind contributions.
   [DOI: 10.1021/ct501131j](https://doi.org/10.1021/ct501131j).
5. C. J. Burnham, D. J. Anick, P. K. Mankoo, and G. F. Reiter, *The vibrational
   proton potential in bulk liquid water and ice*, J. Chem. Phys. **128**,
   154519 (2008). Introduces TTM4-F.
   [DOI: 10.1063/1.2895750](https://doi.org/10.1063/1.2895750).
