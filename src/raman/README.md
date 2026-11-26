# Raman spectra with MB-SPEC

The Raman workflow calculates correlations of the total polarizability tensor
and their reduced Raman spectra. It consists of three executables:

1. `mbspec-raman-polarizability` constructs the total water polarizability tensor for each trajectory frame.
2. `mbspec-raman-response` calculates its isotropic or depolarized autocorrelation.
3. `mbspec-raman-spectrum` Fourier transforms the response and applies the Raman frequency factor.

Copy [the Raman template](../../json_template/raman.json) to your working directory
as `raman.json`, set the trajectory parameters, and run the stages in order.

```bash
OMP_NUM_THREADS=4 mbspec-raman-polarizability trajectory.lammpstrj raman.json > polarizability.dat
mbspec-raman-response polarizability.dat raman.json > response.dat
mbspec-raman-spectrum response.dat raman.json > spectrum.dat
```

## Units

| Quantity | Unit / convention |
|---|---|
| Coordinates | Å |
| Time | ps |
| Polarizabilities | Å³ |
| Polarizability correlations | Å⁶ |
| Temperature | K |
| Transform frequency $\nu$, angular frequency $\omega$ | ps⁻¹ and rad/ps |
| Wavenumber $\widetilde\nu$, smoothing FWHM | cm⁻¹ |
| Raman intensity | Arbitrary units |


## Outputs

| File | Columns, in order |
|---|---|
| `polarizability.dat` | Time [ps], αxx, αxy, αxz, αyy, αyz, αzz [Å³] |
| `response.dat` | Lag time [ps], selected polarizability correlation [Å⁶] |
| `spectrum.dat` | Wavenumber [cm⁻¹], intensity [arbitrary units] |
| `molecular_polarizabilities.dat` | Molecule index, time [ps], αxx, αxy, αxz, αyx, αyy, αyz, αzx, αzy, αzz [Å³] |


## JSON options

```jsonc
{
  "RAMAN": {
    "box": [0.00, 19.66, 0.00, 19.66, 0.00, 19.66], // Cell bounds in Å; XYZ imaging; spectrum requires positive volume
    "imcon": 2,                                 // Cell convention: 1 cubic, 2 orthorhombic (default), 3 triclinic
    "dump_frequency": 0.004,                    // Frame interval in ps
    "skip_initial_atoms": 0,                    // Omit leading atom records; retain O–H–H triplets; default 0
    "trajectory_duration_ps": 0.0,              // Duration read in ps; 0 reads all frames (default)
    "int_max": 10.0,                            // Exclusive correlation lag limit in ps
    "temperature": 300.0,                       // Temperature in K; sets reduced Raman frequency factor
    "spectrum_transform_time_ps": 5.0,          // Response duration transformed, in ps
    "tau": 0.0,                                 // Exponential damping time in ps; 0 disables it (default)
    "smoothing_fwhm_cm1": 15.0,                 // Gaussian spectral FWHM in cm⁻¹; omitted or 0 disables it
    "subtract_mean_polarizability": true,       // Remove mean tensor before correlation; default false
    "write_molecular_properties": false,        // Write molecular tensors; only with 1b; default false
    "raman_response": "depolarized",            // Choose "isotropic" or "depolarized"; default isotropic
    "bose_einstein_correction": false,          // Output I_BE instead of I when true; default false
    "polarizability_mb": "1b",                  // Decomposition; see Polarizability below
    "cutoff_1": 5.5,                            // Inner O–O switching radius in Å; default 5.5
    "cutoff_2": 7.5                             // Outer O–O switching radius in Å; default 7.5
  }
}
```

For `box`, six values are `[xmin,xmax,ymin,ymax,zmin,zmax]`. 
Triclinic cells use nine matrix elements `[xx,xy,xz,yx,yy,yz,zx,zy,zz]`.

## Trajectory input

`mbspec-raman-polarizability` accepts LAMMPS or an XYZ trajectory. For LAMMPS, the
reader identifies coordinates from the `ITEM: ATOMS` header and accepts exactly
one complete triplet of `x y z`, `xu yu zu`, `xs ys zs`, or `xsu ysu zsu`.
XYZ frames contain an atom count, one comment line,
and atom-label/x/y/z rows; 
Water atoms must be ordered O--H--H in every frame.

The internal polarizability model is specific to water. For other molecular
systems, supply a total symmetric polarizability tensor for each frame directly
to `mbspec-raman-response`, then run `mbspec-raman-spectrum`. The external file must contain exactly seven
numeric columns per row, without a header:

```text
time_ps alpha_xx alpha_xy alpha_xz alpha_yy alpha_yz alpha_zz
```


## Polarizability

Raman scattering probes changes in the polarizability of a system. The
polarizability tensor relates an applied electric field to the resulting
induced dipole. MB-SPEC uses the total system tensor,
$\boldsymbol\alpha(t)$, expressed in the laboratory Cartesian frame.

For water, the many-body expansion separates isolated-molecule contributions
from changes caused by interactions between molecules:

$$
\boldsymbol\alpha_{\mathrm{tot}}(1,\ldots,N)=
\sum_i\boldsymbol\alpha_{\mathrm{1B}}(i)
+\sum_{i\lt j}\boldsymbol\alpha_{\mathrm{2B}}(i,j)+\cdots.
$$

The MB-alpha surface uses the Avila one-body polarizability, a short-range
polynomial correction for water dimers, and TTM4-Ind induction.
Its two-body term is

$$
\boldsymbol\alpha_{\mathrm{2B}}(i,j)=
s(r_{ij})\boldsymbol\alpha_{\mathrm{2B,poly}}(i,j)
+\boldsymbol\alpha_{\mathrm{2B,ind}}(i,j).
$$

The polynomial is fitted as a correction to the TTM4-Ind polarizability.
Consequently, the full induced two-body term is retained when the polynomial
is added.

### Available water decompositions

To describe the implemented modes, define

$$
A=\sum_i\boldsymbol\alpha_{\mathrm{1B}}^{\mathrm{Avila}}(i),\qquad
T=\boldsymbol\alpha^{\mathrm{TTM4}}(1,\ldots,N),
$$

$$
T_1=\sum_i\boldsymbol\alpha^{\mathrm{TTM4}}_{\mathrm{1B}}(i),\qquad
P_\alpha=\sum_{i\lt j}s(r_{ij})\boldsymbol\alpha_{\mathrm{2B,poly}}(i,j).
$$

For the truncated `1b+2b` mode, define the switched induction increment

$$
J_\alpha=\sum_{i\lt j}s(r_{ij})
\left[\boldsymbol\alpha^{\mathrm{TTM4}}(i,j)
-\boldsymbol\alpha^{\mathrm{TTM4}}(i)
-\boldsymbol\alpha^{\mathrm{TTM4}}(j)\right].
$$

Each dimer and monomer tensor in this expression is evaluated in isolation.
Subtracting the monomer tensors removes the TTM one-body contribution.

Here, $T$ is the total-system TTM4 tensor and $T_1$ is the sum of
isolated-monomer TTM4 tensors. Their difference, $T-T_1$, is the nonadditive
induction contribution, including two-body and higher-body terms.

| `polarizability_mb` | Total tensor per frame | Contributions |
|---|---|---|
| `1b` | $A$ | Avila one-body |
| `nb` | $T-T_1$ | TTM4 nonadditive |
| `1b+nb` | $A+T-T_1$ | Avila one-body + nonadditive |
| `1b+2b` | $A+P_\alpha+J_\alpha$ | Avila one-body + smoothly truncated complete pair contribution |
| `1b+2b+nb` | $A+P_\alpha+T-T_1$ | Avila one-body + switched 2B polynomial + nonadditive |

The full combination replaces the TTM4 monomer tensors with Avila tensors
and adds the polynomial correction. The `1b+2b` mode includes the polynomial
and isolated-pair induction, switching both to zero at the outer cutoff. It
excludes higher-body induction. This deliberately truncated approximation
differs from the published MB-alpha two-body term in Eq. (2), which retains
long-range induction [Medders and Paesani (2015), Section 2, Eq. (11)](https://doi.org/10.1021/ct501131j).

### Two-body switching

Pair contributions are summed over unique water dimers. The switching
function depends on the O–O distance $r$, with inner radius
$R_i$ (`cutoff_1`) and outer radius $R_o$ (`cutoff_2`):

$$
s(r)=\begin{cases}
1,&r\le R_i,\\
1+x^2(2x-3),&R_i\lt r\lt R_o,\\
0,&r\ge R_o,
\end{cases}
\qquad x=\frac{r-R_i}{R_o-R_i}.
$$

The default radii are 5.5 and 7.5 Å. In `1b+2b`, the same switch multiplies
both the polynomial and the TTM pair increment: both have full weight below
5.5 Å, taper smoothly between 5.5 and 7.5 Å, and vanish at and beyond 7.5 Å.
Only pairs inside the outer cutoff require dimer evaluations; isolated TTM
monomer tensors are calculated once per frame. In `1b+2b+nb`, only the
polynomial is switched; the existing full-system induction is retained.

## Polarizability autocorrelation

The isotropic part of the total tensor describes the average response to an
electric field along the three Cartesian directions. The remaining,
traceless tensor describes its directional anisotropy:

$$
\alpha_{\mathrm{iso}}(t)=\frac13\mathrm{Tr}\boldsymbol\alpha(t),\qquad
\boldsymbol\beta(t)=\boldsymbol\alpha(t)-\alpha_{\mathrm{iso}}(t)\mathbf I,
$$

where $\mathbf I$ is the identity tensor. These quantities give the two
correlations supported by `mbspec-raman-response`:

$$
C_{\mathrm{iso}}(t)=
\left\langle\alpha_{\mathrm{iso}}(0)\alpha_{\mathrm{iso}}(t)\right\rangle,
\qquad
C_{\mathrm{dep}}(t)=
\left\langle\mathrm{Tr}[\boldsymbol\beta(0)\boldsymbol\beta(t)]\right\rangle.
$$

The isotropic correlation tracks
the scalar part of the polarizability; the depolarized correlation tracks
its anisotropic part, which contributes to depolarized Raman scattering.

For $N$ equally spaced frames, `mbspec-raman-response` estimates the equilibrium average
using every available pair separated by $k\Delta t$:

$$
C_{\mathrm{iso}}(k\Delta t)=\frac1{N-k}\sum_{i=0}^{N-k-1}
\alpha_{\mathrm{iso}}(t_i)\alpha_{\mathrm{iso}}(t_{i+k}),
$$

$$
C_{\mathrm{dep}}(k\Delta t)=\frac1{N-k}\sum_{i=0}^{N-k-1}
\mathrm{Tr}[\boldsymbol\beta(t_i)\boldsymbol\beta(t_{i+k})].
$$

Each lag is normalized by its own number of time origins, $N-k$.
For a symmetric tensor, the trace in Eq. (9) sums the three diagonal
products and twice each of the three independent off-diagonal products.
If `subtract_mean_polarizability` is enabled, the correlations use
fluctuations about the trajectory-average tensor:

$$
\langle\boldsymbol\alpha\rangle_{\mathrm{traj}}=
\frac1N\sum_{i=0}^{N-1}\boldsymbol\alpha(t_i),\qquad
\delta\boldsymbol\alpha(t_i)=
\boldsymbol\alpha(t_i)-\langle\boldsymbol\alpha\rangle_{\mathrm{traj}}.
$$

## Fourier transform and Raman spectrum

The Fourier transform resolves the time dependence of the selected
polarizability correlation ($C(t)$) into frequency components. 

$$
\widetilde C(\omega)=\int_{-\infty}^{\infty}
C(t)e^{-i\omega t}\,dt.
$$

For an equilibrium autocorrelation, $C(-t)=C(t)$, so positive lag times
suffice for a cosine transform:

$$
\widetilde C(\omega)=2\int_0^\infty C(t)\cos(\omega t)\,dt.
$$

A reduced Raman spectrum also includes a frequency-dependent
factor involving temperature and the incident laser frequency. MB-SPEC
applies the following form to both response types, defining the reduced
intensity $I(\widetilde\nu)$ before Gaussian smoothing:

$$
I(\widetilde\nu)=
\mathcal A\,\frac{\widetilde\nu}{(\widetilde\nu_I-\widetilde\nu)^4}
\tanh\left(\frac{\hbar\,2\pi c\widetilde\nu}{2k_BT}\right)
S(\widetilde\nu).
$$

Here, $\widetilde\nu$ is the Raman-shift wavenumber in cm⁻¹,
$\widetilde\nu_I=20500\;\mathrm{cm}^{-1}$ is the fixed incident laser wavenumber, 
$c$ is the speed of light in cm/s, 
$T$ is the temperature, 
$k_B$ is the Boltzmann constant, 
$\hbar$ is the reduced Planck constant,
$\mathcal A$ is a frequency-independent constant setting the arbitrary intensity scale, and
$S$ is the finite cosine transform defined below.

Set `bose_einstein_correction` to `true` to output the Bose–Einstein-corrected
spectrum $I_{\mathrm{BE}}(\widetilde\nu)$, following
[Iuchi, Morita, and Kato (2002), Eq. (26)](https://doi.org/10.1021/jp013773y):

$$
I_{\mathrm{BE}}(\widetilde\nu)=
\frac{I(\widetilde\nu)}{2\pi c\widetilde\nu},
$$

The correction applies to either response type before Gaussian smoothing. 
The output is reported in arbitrary units, with no volume prefactor or absolute cross-section
normalization. 

### Spectral processing

Exponential damping suppresses long-time contributions and broadens spectral features. 
It acts on the correlation before Fourier transformation.
For $C_k=C(k\Delta t)$, the samples used by the transform are

$$
D_k=C_k
\begin{cases}
\exp(-k\Delta t/\tau),&\tau>0,\\
1,&\tau=0.
\end{cases}
$$


For transform duration $T_{\mathrm{transform}}$, the sample count is

$$
N_T=\left\lfloor\frac{T_{\mathrm{transform}}}{\Delta t}+10^{-10}\right\rfloor.
$$

Here, $T_{\mathrm{transform}}$ is the JSON parameter `spectrum_transform_time_ps`.
The small dimensionless offset prevents rounding just below an integer from
dropping a sample. MB-SPEC evaluates

$$
S(\widetilde\nu_n)=2\Delta t\left[
\frac12D_0+\sum_{k=1}^{N_T-2}D_k
\cos(2\pi c_{\mathrm{cm/ps}}\widetilde\nu_n k\Delta t)\right],
$$

on the grid

$$
\Delta\nu=\frac1{(2N_T-1)\Delta t},\qquad
\nu_n=n\Delta\nu,\qquad\omega_n=2\pi\nu_n,\qquad
\widetilde\nu_n=\nu_n/c_{\mathrm{cm/ps}},
\quad n=0,\ldots,N_T-1.
$$

$S$ has units Å⁶·ps.
The last selected sample, $D_{N_T-1}$, is excluded from the sum; the
implementation assumes the correlation has decayed to zero at the upper
integration limit.

Gaussian smoothing acts on the selected intensity, $I$ or $I_{\mathrm{BE}}$.
A positive
`smoothing_fwhm_cm1` sets its full width at half maximum (FWHM) in cm⁻¹;
omission or zero disables it. For FWHM $w$ and wavenumber spacing
$\Delta\widetilde\nu$, the kernel width and length are

$$
\sigma_{\mathrm{grid}}=\frac{w}{2\sqrt{2\ln2}\,\Delta\widetilde\nu},
\qquad K=2\left\lceil3\sigma_{\mathrm{grid}}\right\rceil+1,
$$

$$
g_j=\frac{\exp[-(j-(K-1)/2)^2/(2\sigma_{\mathrm{grid}}^2)]}
{\sum_{m=0}^{K-1}\exp[-(m-(K-1)/2)^2/(2\sigma_{\mathrm{grid}}^2)]}.
$$

The Gaussian is truncated at three standard deviations.

## Scientific limitations

**Periodic induction.** 
The `nb`, `1b+nb`, and `1b+2b+nb` modes evaluate $T$ using an internal TTM4 
response that uses only orthorhombic cell lengths for imaging. These modes 
lack Ewald treatment and have not been validated for periodic bulk simulations. 
`mbspec-raman-polarizability` accepts triclinic cells (`imcon=3`) only for the internal `1b` and `1b+2b` modes.

