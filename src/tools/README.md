# Shared library and water property models

`src/tools` provides the internal `libtools.a` library used by the IR, Raman,
SFG, and power spectrum programs. It contains trajectory readers, configuration
handling, geometry and numerical utilities, and the water dipole and
polarizability models. 

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
| `poly2b5.cpp/.h` | Polynomial basis evaluation used by `dip_poly2b5`. |
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
The `dipind` accessors sum induced dipoles for the supplied system.

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
  output uses `xx, xy, xz, yy, yz, zz`.
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

- The built-in property surfaces are water-specific. External dipole definitions
are determined by the trajectory producer; PS, LTP2011, and externally supplied
permanent dipoles must not be treated as interchangeable without justification.
The published MB-μ construction uses LTP2011 for its explicit one-body term
[4]. The optional PS permanent term in SFG is a different model choice, and the
PS-derived charges inside TTM4 induction serve a separate role.

- The internal cell-dependent TTM4 treatment uses only the diagonal cell lengths
for minimum imaging and therefore assumes an orthorhombic cell, including
cubic cells. Triclinic support in the shared trajectory reader does not extend
to this model. It has no Ewald summation. Raman NB-containing modes are not 
validated for periodic bulk calculations. 

## Build and generated files

- `Makefile.am` defines `libtools.a`, its sources, and additional distribution
files. The library is built as part of the repository build. 

- The current matrix utility implementation requires Intel MKL's `mkl_lapacke.h` and LAPACKE
libraries. Configuration checks compilation and linkage of `LAPACKE_dgetrf`,
`LAPACKE_dgetri`, and `LAPACKE_dsyev` and stops if they are unavailable.
`--disable-mkl` is unsupported. Set `--with-mkl-cflags` and
`--with-mkl-libs` when the installation differs from the default `/opt/intel`
paths. Build configuration supplies `MKL_CFLAGS` and `MKL_LIBS` to the 
spectroscopy programs using this library.
A generic LAPACKE installation alone does not provide the named MKL header.

- `Makefile.in`, `Makefile`, `.deps/`, `*.o`, and `libtools.a` are generated build
files.

- The large `poly2b5.cpp` and `poly-l-2b3.cpp` files contain explicit polynomial
basis expressions; their fitted coefficients reside in the corresponding
`dip_poly2b5.cpp` and `alpha-el-poly2b3.cpp` model files.

# TTM4 dipole moment and polarizability

TTM4 represents water electrostatics using geometry-dependent permanent charges
and mutually induced dipoles. The dipole moment describes the system's charge
separation; the polarizability describes how that dipole changes under an
applied electric field.

MB-SPEC implements **TTM4-Ind**, the TTM4-based induction model described by 
[Medders and Paesani (2015)](https://doi.org/10.1021/ct501131j), rather than the full **TTM4-F** 
[water potential](https://doi.org/10.1063/1.2895750). 

## Dipole calculation

Each water molecule has four electrostatic sites: O, H, H, and a virtual M site
near the oxygen. MB-SPEC obtains geometry-dependent charges from the
Partridge–Schwenke dipole surface and redistributes them onto the H and M sites.

The permanent charges produce an electric field at the polarizable sites.
Each site develops an induced dipole:

$$
\boldsymbol{\mu}_{i}^{\mathrm{ind}}
=a_i\left[
\mathbf{E}_{i}^{\mathrm{perm}}
+\sum_{j\ne i}\mathbf{T}_{ij}^{\mathrm{damped}} {\boldsymbol{\mu}_{j}^{\mathrm{ind}}}
\right].
$$

Here, $a_i$ is the site polarizability, $\mathbf E_i^{\mathrm{perm}}$ is the
field from the permanent charges, and $\mathbf T_{ij}^{\mathrm{damped}}$
describes the damped dipole interaction.

The dipoles are coupled: each induced dipole contributes to the field inducing
the others. MB-SPEC solves these simultaneous linear equations using a
Cholesky decomposition.

For the polarizable O and H sites, the damped interaction tensor is

$$
\mathbf T_{ij}^{\mathrm{damped}}
=\frac{3\left[1-\left(1+\frac43x_{ij}\right)e^{-x_{ij}}\right]}
{r_{ij}^{5}}\,
\mathbf r_{ij}\mathbf r_{ij}^{\mathsf T}
-\frac{1-e^{-x_{ij}}}{r_{ij}^{3}}\mathbf I,
$$

with

$$
x_{ij}=d_{ij}
\left[\frac{r_{ij}}{(a_i a_j)^{1/6}}\right]^4.
$$

Here, $\mathbf r_{ij}=\mathbf r_i-\mathbf r_j$,
$r_{ij}=|\mathbf r_{ij}|$, $\mathbf I$ is the identity tensor, and $d_{ij}$
is the dimensionless damping parameter. The implementation uses
$d_{ij}=0.055$ for intermolecular pairs and intramolecular H–H pairs, and
$d_{ij}=0.626$ for intramolecular O–H pairs. The M site has zero polarizability
and carries no induced dipole.

At large separation, the exponential terms vanish and the tensor reduces to
the usual point-dipole interaction. At short separation, damping suppresses
the singular behavior. These expressions use the model's electrostatic units;
with distances in Å and site polarizabilities in Å³, the tensor has units of
Å⁻³. See `smear2` in [ttm4-smear.cpp](tools/ttm4-smear.cpp) and the tensor
construction in [ttm4-dip.cpp](tools/ttm4-dip.cpp).

The total dipole is

$$
\mathbf M^{\mathrm{TTM4}}
=\underbrace{\sum_i q_i\mathbf r_i}_{\text{permanent}}
+\underbrace{\sum_i\boldsymbol\mu_i^{\mathrm{ind}}}_{\text{induced}}.
$$

`ttm4_es::dip` returns permanent plus
induced dipoles, whereas `dipind` and the separate `ttm4_dip` routine return
only the induced contribution. Dipole outputs are in Debye. See
[ttm4-es.cpp](tools/ttm4-es.cpp) and [ttm4-dip.cpp](tools/ttm4-dip.cpp).

## Polarizability calculation

At fixed atomic positions, a small uniform external field $\mathbf F$ changes
the induced dipoles throughout the system. The system polarizability tensor
measures this response:

$$
A_{ab}^{\mathrm{TTM4}}
=\left.\frac{\partial M_a}{\partial F_b}\right|_{\mathbf F=0}.
$$

When polarizability is expressed in SI units, this derivative has units of
C·m²/V. The Å³ convention used by MB-SPEC expresses polarizability as a
volume, obtained by dividing the SI polarizability by $4\pi\epsilon_0$ and
converting m³ to Å³. The coupled-dipole equations here are understood in the
model's consistent electrostatic units.

Differentiating the induced-dipole equations gives

$$
\frac{\partial\boldsymbol\mu_i^{\mathrm{ind}}}{\partial F_b}
=a_i\left[
\mathbf e_b+
\sum_{j\ne i}\mathbf T_{ij}^{\mathrm{damped}}
\frac{\partial\boldsymbol\mu_j^{\mathrm{ind}}}{\partial F_b}
\right],
$$

where $\mathbf e_b$ is the unit vector along Cartesian direction $b$.
Polarizability therefore follows from the same coupled dipole response.
MB-SPEC evaluates it analytically from the polarization matrix. It does not
perform separate finite-field calculations. See `ttm4_es::ptensor` in
[ttm4-es.cpp](tools/ttm4-es.cpp).

The site parameters are

$$
a_{\mathrm O}=1.310\ \text{Å}^3,\qquad
a_{\mathrm H}=0.294\ \text{Å}^3,\qquad
a_{\mathrm M}=0.
$$

Although these site parameters are isotropic and constant, their interactions
produce a geometry-dependent molecular polarizability tensor,
reported in Å³.

## Use in MB-SPEC

The total TTM4 polarizability already contains isolated-molecule
contributions. Raman's nonadditive contribution therefore subtracts them:

$$
\mathbf A_{\mathrm{nonadd}}^{\mathrm{TTM4}}
=\mathbf A_{\mathrm{system}}^{\mathrm{TTM4}}
-\sum_m\mathbf A_m^{\mathrm{TTM4}}.
$$

The full implemented Raman combination is

$$
\mathbf A^{\mathrm{MB\text{-}\alpha}}
=\sum_m\mathbf A_m^{\mathrm{Avila}}
+\mathbf A_{\mathrm{2B,corr}}
+\mathbf A_{\mathrm{nonadd}}^{\mathrm{TTM4}}.
$$

This subtraction prevents counting both the TTM4 and Avila monomer
polarizabilities. The nonadditive term includes two-body induction as well as
higher-body effects. See [polarizability.cpp](raman/polarizability.cpp) and
[alpha-ttm4.cpp](tools/alpha-ttm4.cpp).

The explicit short-range polarizability correction is added on top of
induction, whereas the MB-mu dipole construction replaces the corresponding
short-range induced-dimer contribution. See
[Medders and Paesani, Eqs. (7), (10), and (11)](https://doi.org/10.1021/ct501131j).

The internal cell-dependent implementation lacks Ewald electrostatics and is
not validated for periodic bulk calculations.
## References

1. Partridge and Schwenke, J. Chem. Phys. 106, 4618 (1997)
2. Lodi, Tennyson, and Polyansky, J. Chem. Phys. 135, 034113 (2011)
3. Avila, J. Chem. Phys. 122, 144310 (2005)
4. Medders and Paesani, J. Chem. Theory Comput. 11, 1145 (2015)
5. Burnham, Anick, Mankoo, and Reiter, J. Chem. Phys. 128, 154519 (2008).
