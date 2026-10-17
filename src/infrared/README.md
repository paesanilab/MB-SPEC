# Infrared spectra with MB-SPEC

The infrared workflow calculates the total dipole time correlation and its absorption spectrum. It consists of three executables:

1. dipole   : construct total dipole for each trajectory frame;
2. response : calculate the total dipole autocorrelation function; and
3. spectrum : Fourier transform the response and apply the IR normalization.

Run the three stages in order:

```bash
OMP_NUM_THREADS=4 ./dipole trajectory.lammpstrj mb-spec.json > dipole.dat
./response dipole.dat mb-spec.json > response.dat
./spectrum response.dat mb-spec.json > spectrum.dat
```

## Units

| Quantity | Unit / convention |
|---|---|
| Coordinates | Å |
| Time | ps; converted to seconds |
| Dipole moments | Debye; converted to C·m |
| Dipole autocorrelation | (C·m)² |
| Temperature | K |
| Frequency $\nu$, angular frequency $\omega$ | Hz and rad/s |
| Wavenumber $\widetilde\nu$, smoothing FWHM | cm⁻¹ |
| IR intensity $n(\omega)\alpha(\omega)$ | cm⁻¹ |

## Outputs

| File | Columns, in order |
|---|---|
| `dipole.dat` | Time [ps], μx, μy, μz [Debye] |
| `response.dat` | Lag time [ps], total-dipole autocorrelation [(C·m)²] |
| `spectrum.dat` | Wavenumber [cm⁻¹], n(ω)α(ω) [cm⁻¹] |
| `molecular_dipoles.dat` (optional) | Molecule index, time [ps], μx, μy, μz [Debye] |


## JSON options

Put all settings inside `INFRARED`.
**remove the `//` comments before saving it as an input JSON file.**

```jsonc
{
  "INFRARED": {
    "box": [0.0,19.65,0.0,19.65,0.0,19.65], // Required cell bounds in Å; sets IR normalization volume
    "imcon": 2,                             // Cell convention: 1 cubic, 2 orthorhombic (default), 3 triclinic
    "dump_frequency": 0.004,                // Required frame interval in ps
    "skip_initial_atoms": 0,                // Omit this many leading atom records per frame; default 0
    "trajectory_duration_ps": 0.0,          // Duration read in ps; 0 reads all frames (default)
    "int_max": 10.0,                        // Required exclusive correlation lag limit in ps
    "temperature": 300.0,                   // Required temperature in K; sets IR intensity prefactor
    "spectrum_transform_time_ps": 10.0,     // Required response duration in ps
    "tau": 0.0,                             // Exponential damping time in ps; 0 disables it (default)
    "smoothing_fwhm_cm1": 15.0,             // Gaussian spectral FWHM in cm⁻¹; 0 disables it
    "subtract_mean_dipole": false,          // Remove the mean total dipole before correlation; default false
    "write_molecular_properties": false,    // Write molecular_dipoles.dat; only DMS=true with 1b; default false
    "DMS": false,                           // false: external dipoles (default); true: internal LTP2011 water surface
    "dipole_mb": "1b",                      // Required dipole source
    "cutoff_1": 5.5,                        // Inner O–O switching radius in Å for 2B terms; default 5.5
    "cutoff_2": 7.5                         // Outer O–O switching radius in Å for 2B terms; default 7.5
  }
}
```

For `box`, six values are `[xmin,xmax,ymin,ymax,zmin,zmax]`; triclinic cells
use nine matrix elements `[xx,xy,xz,yx,yy,yz,zx,zy,zz]`.

## Trajectory input

For `DMS=false`, input must be a standard LAMMPS trajectory. 
The reader identifies fields from the `ITEM: ATOMS` header. It requires
one complete coordinate triplet: `x y z`, `xu yu zu`, `xs ys zs`, or
`xsu ysu zsu`. External dipole modes additionally require
`v_dipx_perm v_dipy_perm v_dipz_perm` and
`v_dipx_ind v_dipy_ind v_dipz_ind`.
Other fields, including `id`, `mol`, `type`, `q`, and dipole magnitudes, are
optional. 

For `DMS=true`, either the same LAMMPS format or a simple XYZ trajectory may be
used. Each XYZ frame consists of an atom count, one comment line, and atom
label/x/y/z lines. The JSON box supplies the cell for XYZ
input. Water atoms must be ordered O--H--H.


## Dipole moments

With `DMS=false`, `dipole` reads externally supplied permanent and induced dipole moments:

| dipole_mb | Total dipole per frame | Contributions |
|---|---|---|
| `1b` | $\mathbf M_{\mathrm{perm}}$ | Supplied permanent dipoles |
| `nb` | $\mathbf M_{\mathrm{ind}}$ | Supplied induced dipoles |
| `1b+nb` | $\mathbf M_{\mathrm{perm}}+\mathbf M_{\mathrm{ind}}$ | Permanent + induced |
| `1b+2b` | $\mathbf M_{\mathrm{perm}}+P$ | Permanent + switched 2B polynomial |
| `1b+2b+nb` | $\mathbf M_{\mathrm{perm}}+\mathbf M_{\mathrm{ind}}+P-I$ | Permanent + induced; replace switched induced-dimer term with polynomial |


With `DMS=true`, `dipole` evaluates the internal LTP2011 water surface:

| dipole_mb | Total dipole per frame | Contributions |
|---|---|---|
| `1b` | $\sum_i\boldsymbol\mu_i^{\mathrm{LTP2011}}$ | LTP2011 one-body dipoles |
| `1b+2b` | $\sum_i\boldsymbol\mu_i^{\mathrm{LTP2011}}+P$ | LTP2011 one-body + switched 2B polynomial |

The permanent and induced dipole moments are summed over atoms as:

$$
\mathbf M_{\mathrm{perm}}=\sum_a\boldsymbol\mu_a^{\mathrm{perm}},
\qquad\mathbf M_{\mathrm{ind}}=\sum_a\boldsymbol\mu_a^{\mathrm{ind}}.
\tag{1}
$$

The switched polynomial $P$ and induced replacement $I$ are defined for unique water dimers as:

$$
P\equiv\Delta\mathbf M_{\mathrm{2B,poly}}
=\sum_{i\lt j}s(r_{ij})\boldsymbol\mu_{ij}^{\mathrm{2B,poly}},
\tag{2}
$$

$$
I\equiv\Delta\mathbf M_{\mathrm{2B,ind}}
=\sum_{i\lt j}s(r_{ij})\boldsymbol\mu_{ij}^{\mathrm{2B,ind}}.
\tag{3}
$$

The short-range polynomial dipole is constructed from effective charges on
the six atoms of each water dimer:

$$
\boldsymbol\mu_{ij}^{\mathrm{2B,poly}}
=\sum_{a=1}^{6}q_a(\xi_0,\ldots,\xi_{14})\,\mathbf r_a.
$$

The effective charges are fifth-order polynomials in 15 auxiliary variables formed from exponentiated interatomic distances. 
The basis contains 1,826 terms for hydrogen charges and 1,052 for oxygen charges, with symmetry under permutations of equivalent atoms. 
The fit used 42,495 dimer configurations with counterpoise-corrected RI-MP2/aug-cc-pVTZ two-body dipoles. 
`dip_poly2b5::eval` returns the unswitched polynomial dipole and the IR caller applies $s(r_{ij})$. 

The explicit two-body switching function is

$$
s(r)=\begin{cases}
1,&r\le R_i,\\
1+x^2(2x-3),&R_i\lt r\lt R_o,\\
0,&r\ge R_o,
\end{cases}
\qquad x=\frac{r-R_i}{R_o-R_i}.
\tag{4}
$$

The full two-body term is

$$
\boldsymbol\mu_{ij}^{\mathrm{2B}}
=s(r_{ij})\boldsymbol\mu_{ij}^{\mathrm{2B,poly}}
+[1-s(r_{ij})]\boldsymbol\mu_{ij}^{\mathrm{2B,ind}}.
\tag{5}
$$

With higher-body induction,

$$
\mathbf M_{\mathrm{MB-\mu}}=
\sum_i\boldsymbol\mu_i^{\mathrm{1B}}
+\sum_{i\lt j}\left[s(r_{ij})\boldsymbol\mu_{ij}^{\mathrm{2B,poly}}
+(1-s(r_{ij}))\boldsymbol\mu_{ij}^{\mathrm{2B,ind}}\right]
+\mathbf M_{\mathrm{higher}}^{\mathrm{ind}}.
\tag{6}
$$

The full `1b+2b+nb` route implements

$$
\mathbf M_{\mathrm{MB-\mu}}=\mathbf M_{\mathrm{perm}}+\mathbf M_{\mathrm{ind}}
+\sum_{i\lt j}s(r_{ij})\left[\boldsymbol\mu_{ij}^{\mathrm{2B,poly}}
-\boldsymbol\mu_{ij}^{\mathrm{2B,ind}}\right].
\tag{7}
$$



## Dipole autocorrelation

Infrared absorption is calculated from the time-dependent total dipole moment
of the system, $\mathbf M(t)$. The selected dipole model or external input
provides this quantity; the subsequent correlation and spectrum calculations
follow the same procedure for all dipole sources.

The equilibrium total-dipole autocorrelation function is

$$
C(t)=\left\langle\mathbf M(0)\cdot\mathbf M(t)\right\rangle.
\tag{8}
$$

The brackets denote an equilibrium average, and the dot product sums the
correlations of the three Cartesian dipole components. The correlation
describes how the total dipole remains correlated over a time separation $t$.
Its Fourier transform resolves this temporal behavior into frequency components.

For a trajectory containing $N$ equally spaced frames, `response` estimates
this average using every available pair of frames separated by $k\Delta t$:

$$
C(k\Delta t)=\frac1{N-k}\sum_{i=0}^{N-k-1}
\mathbf M(t_i)\cdot\mathbf M(t_{i+k}).
\tag{9}
$$

Each lag is normalized by its own number of available time origins, $N-k$.
Longer lags therefore contain fewer contributions to the average.

If `subtract_mean_dipole` is enabled, the correlation is calculated from
fluctuations about the trajectory mean, replacing $\mathbf M$ in Eq. (9) by
$\delta\mathbf M$:

$$
\overline{\mathbf M}=\frac1N\sum_{i=0}^{N-1}\mathbf M(t_i),
\qquad\delta\mathbf M(t_i)=\mathbf M(t_i)-\overline{\mathbf M}.
\tag{10}
$$

## Fourier transform and IR spectrum

The frequency-domain dipole correlation is defined by

$$
\widetilde C(\omega)=\int_{-\infty}^{\infty}
C(t)e^{-i\omega t}\,dt.
\tag{11}
$$

For an equilibrium autocorrelation function, $C(-t)=C(t)$, so the transform
can be evaluated from positive lag times using a cosine transform:

$$
\widetilde C(\omega)=2\int_0^\infty C(t)\cos(\omega t)\,dt.
\tag{12}
$$

MB-SPEC converts the transformed correlation into the infrared
absorption quantity using

$$
n(\omega)\alpha(\omega)=
\frac{\beta\omega^2}{6\epsilon_0cV}\widetilde C(\omega),
\qquad\beta=\frac1{k_BT}.
\tag{13}
$$

where 
$n(\omega)$ is the frequency dependent refractive index, 
$\alpha(\omega)$ is the Beer-Lambert absorption coefficient, 
$T$ is the temperature, 
$V$ is the system volume, 
$\epsilon_0$ is the vacuum permittivity, 
and $c$ is the speed of light. 
Equation (13) can also be interpreted using the Kubo-transformed dipole
correlation of quantum linear-response theory, with the classical correlation
used as an approximation. 
MB-SPEC applies no additional frequency-dependent quantum correction factor.

### Spectral processing

A finite trajectory provides the correlation at discrete lag times over a
limited duration. `spectrum_transform_time_ps` selects the number of samples
used to approximate the continuous transform. For $N_T$ samples, with
$\Delta t$ expressed in seconds, MB-SPEC evaluates

$$
\widetilde C_n=\Delta t\left[C_0+2\sum_{k=1}^{N_T-1}
C_k\cos\left(\frac{2\pi nk}{2N_T-1}\right)\right],
\qquad C_k=C(k\Delta t).
\tag{14}
$$

The corresponding frequency grid is

$$
\Delta\nu=\frac1{(2N_T-1)\Delta t},\qquad
\nu_n=n\Delta\nu,\qquad\omega_n=2\pi\nu_n,
\qquad n=0,\ldots,N_T-1.
\tag{15}
$$

Two optional processing steps affect the resulting spectrum. 
1. Exponential damping `tau` acts on the correlation before Fourier transform, suppresses long-time contributions, and broadens spectral features. `tau = 0` disables damping.

$$
C_k\leftarrow C_k\exp\left(-\frac{k\Delta t}{\tau}\right).
\tag{17}
$$

2. Gaussian smoothing acts on the completed spectrum. 
A positive `smoothing_fwhm_cm1` applies Gaussian spectral smoothing with 
the specified full width at half maximum (FWHM) in cm⁻¹. Zero disables smoothing.
For FWHM $w$ and wavenumber spacing
$\Delta\widetilde\nu$, the kernel width $\sigma_{\mathrm{grid}}$, length $K$, and normalized kernel $g_j$ are

$$
\sigma_{\mathrm{grid}}=\frac{w}{2\sqrt{2\ln2}\,\Delta\widetilde\nu},
\qquad K=2\left\lceil3\sigma_{\mathrm{grid}}\right\rceil+1.
\tag{21}
$$

$$
g_j=\frac{\exp\left[-(j-(K-1)/2)^2/(2\sigma_{\mathrm{grid}}^2)\right]}
{\sum_{m=0}^{K-1}\exp\left[-(m-(K-1)/2)^2/(2\sigma_{\mathrm{grid}}^2)\right]}.
\tag{22}
$$
