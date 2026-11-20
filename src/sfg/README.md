# Sum-frequency generation spectra with MB-SPEC

The SFG workflow calculates a truncated dipole-polarizability cross-correlation function for water.
It requires a water slab trajectory with the permanent and induced dipole columns required by the reader.

Run:

```bash
OMP_NUM_THREADS=4 ./response trajectory.lammpstrj sfg.json > response.dat
./spectrum response.dat sfg.json > spectrum.dat
```

Optionally apply Fresnel factors after calculating the component spectra:

```bash
./fresnel fresnel.json
```


## Units

| Quantity | Unit / convention |
|---|---|
| Cartesian coordinates | Å |
| Time | ps |
| Permanent and induced dipoles | Debye |
| Polarizabilities | Å³ |
| Molecular response | D·Å³ |
| Temperature | K|
| Frequency $\nu$, angular frequency $\omega$ | ps⁻¹ and rad/ps  |
| Wavenumber $\widetilde\nu$, smoothing FWHM, Fresnel redshift | cm⁻¹ |
| Reduced complex signal | D·Å³ × ps × cm⁻¹ |
| Fresnel wavelengths and optical frequencies | µm and GHz |
| Fresnel incidence angles (in JSON) | Degrees from the surface normal |
| Refractive index and extinction coefficient | Dimensionless |

## Outputs

| File | Columns, in order |
|---|---|
| `response.dat` | Lag time [ps], response [D·Å³], unnormalized sum [D·Å³], time-origin count |
| `spectrum.dat` | Wavenumber [cm⁻¹], imaginary reduced signal, real reduced signal, squared magnitude |
| Fresnel polarization files | Shifted wavenumber [cm⁻¹], squared complex amplitude |


## JSON options

```jsonc
{
  "SFG": {
    "box": [0.0,24.79,0.0,24.79,0.0,100.0],      // Reference cell in Å; unused by current SFG stages
    "imcon": 2,                                  // Cell convention: 1 cubic, 2 orthorhombic (default), 3 triclinic
    "dump_frequency": 0.004,                     // Frame interval in ps
    "skip_initial_atoms": 0,                     // Omit leading records; retain O–H–H triplets; default 0
    "trajectory_duration_ps": 0.0,               // Duration read in ps; 0 reads all frames (default)
    "temperature": 300.0,                        // Positive K; validated but unused in reduced transform
    "int_max": 10.0,                             // Maximum response lag in ps
    "dipole_mb": "1b",                           // Dipole decomposition; see Molecular properties below
    "polarizability_mb": "1b",                   // 1b or 1b+2b; no NB polarizability mode
    "cutoff_1": 5.5,                             // Inner O–O switching radius in Å; default 5.5
    "cutoff_2": 7.5,                             // Outer O–O switching radius in Å; default 7.5
    "ijk": "xxz",                                // Tensor component: zzz, zxx, xzx, or xxz
    "zc1": 6.0,                                  // Inner interface-weight radius in Å
    "zc2": 8.0,                                  // Outer interface-weight radius in Å
    "surface": "both",                           // positive, negative, or both (default); signed lower-interface weight
    "center_slab": true,                         // Center mass-weighted slab COM at origin; default false
    "shift_slab_by_half_cell": false,            // Shift/reimage by half z cell before centering; default false
    "time_dependent_switch": false,              // Same-molecule weights: g(0)g(t)² if true, g(0)³ if false (default)
    "use_ps_dipole": false,                      // Replace external permanent dipoles with PS; default false
    "include_autocorrelation": true,             // Include same-molecule terms; default true
    "include_cross_correlation": false,          // Include selected neighbor terms; default true
    "cross_correlation_cutoff": 4.0,             // Positive O–O neighbor cutoff in Å; default 4.0
    "cross_correlation_sampling_stride": 250,    // Positive neighbor-sampling frame stride; default 250
    "subtract_initial_response": true,           // Boolean: subtract C(0) before transform if true
    "spectrum_transform_time_ps": 4.0,           // Response duration transformed, in ps
    "tau": 0.0,                                  // Exponential damping time in ps; 0 disables it (default)
    "smoothing_fwhm_cm1": 10.0,                  // Gaussian FWHM in cm⁻¹; omitted/0 disables; smooth complex components
    "spectrum_max_wavenumber_cm1": 4500.0,       // Grid limit in cm⁻¹; points above Nyquist are aliased
    "spectrum_num_wavenumbers": 4500             // Integer ≥ 2; grid point count, including both endpoints
  }
}
```

## Trajectory input

`response` accepts one LAMMPS dump trajectory.
The reader identifies fields from the `ITEM: ATOMS` header.
It requires exactly one coordinate triplet of `x y z`, `xu yu zu`, `xs ys zs`, or `xsu ysu zsu`;
scaled coordinates are converted with the LAMMPS box matrix and origin.
It also requires `v_dipx_perm v_dipy_perm v_dipz_perm` and `v_dipx_ind v_dipy_ind v_dipz_ind`.
Frame times are frame index times `dump_frequency`; LAMMPS timestep records
do not supply the time axis.

## Molecular properties

Before explicit pair corrections,

$$
\mu_{a,k}=I_{\mathrm{1B}}\mu_{a,k}^{\mathrm{perm}}
+I_{\mathrm{NB}}\mu_{a,k}^{\mathrm{ind}},
$$

$$
\alpha_{a,ij}=\alpha_{a,ij}^{\mathrm{Avila,1B}}.
$$

The permanent term comes from the trajectory or internal PS surface when
`use_ps_dipole=true`; there is no internal LTP2011 option. The indicators
$I_{\mathrm{1B}}$ and $I_{\mathrm{NB}}$ select the base contributions:

| `dipole_mb` | $I_{\mathrm{1B}}$ | $I_{\mathrm{NB}}$ | Unswitched pair term added through allocation below |
|---|---|---|---|
| `1b` | 1 | 0 | None |
| `nb` | 0 | 1 | None |
| `1b+nb` | 1 | 1 | None |
| `1b+2b` | 1 | 0 | $\boldsymbol\mu_{ab}^{\mathrm{2B,poly}}$ |
| `1b+2b+nb` | 1 | 1 | $\boldsymbol\mu_{ab}^{\mathrm{2B,poly}}-\boldsymbol\mu_{ab}^{\mathrm{2B,ind}}$ |

| `polarizability_mb` | Base tensor | Unswitched pair term added through allocation below |
|---|---|---|
| `1b` | Avila one-body | None |
| `1b+2b` | Avila one-body | $\boldsymbol\alpha_{ab}^{\mathrm{2B,poly}}+\boldsymbol\alpha_{ab}^{\mathrm{2B,ind}}$ |

The polarizability induction increment is evaluated for an isolated dimer:

$$
\boldsymbol\alpha_{ab}^{\mathrm{2B,ind}}
=\boldsymbol\alpha^{\mathrm{TTM4}}(a,b)
-\boldsymbol\alpha^{\mathrm{TTM4}}(a)
-\boldsymbol\alpha^{\mathrm{TTM4}}(b).
$$

The monomer subtraction avoids counting one-body polarizabilities twice.
Isolated TTM monomer tensors are cached once per frame. The polynomial is an
additive correction to this TTM pair increment.
There is no internal nonadditive-polarizability mode. Both external dipole
triplets remain required for every decomposition.

For an eligible pair $a,b$, the switched correction allocated to molecule $a$ is

$$
\Delta X_{a(ab)}=s(r_{ab})\frac{|g_a|}{|g_a|+|g_b|}\Delta X_{ab}^{\mathrm{2B}},
\qquad X\in\{\mu_k,\mu_{k'},\alpha_{ij},\alpha_{i'j'}\}.
$$

Allocation to $b$ is analogous; pairs with zero denominator are skipped. For dipole `1b+2b+nb`, the pair term is the MB-mu polynomial minus the TTM4 induced-dimer contribution. For polarizability `1b+2b`, it is the sum of the MB-alpha polynomial correction and the TTM pair increment. Both use the same switch and molecular allocation. For eligible pairs, the two allocated contributions sum to the switched pair tensor before the correlation applies additional interface weights. The allocation is an existing molecular-response convention; conservation of the pair tensor alone does not establish a unique interfacial response.

### Two-body switching

The explicit two-body switching function is

$$
s(r)=\begin{cases}
1,&r\le R_i,\\
1+x^2(2x-3),&R_i\lt r\lt R_o,\\
0,&r\ge R_o,
\end{cases}
\qquad x=\frac{r-R_i}{R_o-R_i}.
$$

The switching function vanishes at the outer radius. For polarizability
`1b+2b`, both polynomial and pair induction have full weight below the inner
radius (default 5.5 Å), taper smoothly to zero at the outer radius (default
7.5 Å), and are omitted at and beyond that radius. 
It excludes higher-body induction.

## Dipole-polarizability correlation

### Interface weights and tensor selection

Molecular oxygen coordinates are $\mathbf r_a$ with laboratory coordinate $z_a$. 
For $z_1=\mathrm{zc1}$, $z_2=\mathrm{zc2}$, and $0\le z_1\lt z_2$,

$$
g(z)=\begin{cases}
0,&|z|\lt z_1,\\
\mathrm{sgn}(z)\cos^2\left[\frac{\pi(|z|-z_2)}{2(z_1-z_2)}\right],&z_1\le|z|\le z_2,\\
\mathrm{sgn}(z),&|z|\gt z_2.
\end{cases}
$$

For dipole origin $m$ and later polarizability frame $n$,

$$
B_{ab}^{mn}=\frac12\left[\mu_{a,k}^{m}\alpha_{b,ij}^{n}
+\mu_{a,k'}^{m}\alpha_{b,i'j'}^{n}\right].
$$

| Component | $(k,ij)$ | $(k',i'j')$ |
|---|---|---|
| `zzz` | $(z,zz)$ | $(z,zz)$ |
| `zxx` | $(x,zx)$ | $(y,zy)$ |
| `xzx` | $(x,xz)$ | $(y,yz)$ |
| `xxz` | $(z,xx)$ | $(z,yy)$ |

For example, `xxz` uses $[\mu_z\alpha_{xx}+\mu_z\alpha_{yy}]/2$.

### Truncated cross-correlation partners

For each origin $m$ and central molecule $a$ with $g_a^m\ne0$, candidates $b\ne a$ first require $r_{ab}(m)\lt 2R$, where $R$ is the `cross_correlation_cutoff`. A partner is retained if

$$
\overline r_{ab,m}=\frac1Q\sum_{\substack{q=0,h,2h,\ldots\\m+q\le L_m}}
|\mathbf r_a(m)-\mathbf r_b(m+q)|\lt R,
\qquad L_m=\min(N-1,m+N_{\mathrm{lag}}).
$$

Here, $h$ is the `cross_correlation_sampling_stride` and $Q$ counts sampled terms. The central oxygen is held at time $m$. Distances are minimum-imaged using the cell at the origin frame. Holding the central oxygen fixed defines an asymmetric, time-dependent partner selection and truncates the full system correlation.

### Accumulation and normalization

$$
N_{\mathrm{lag}}=\left\lfloor\frac{\mathrm{int\_max}}{\Delta t}\right\rfloor.
$$

For $\ell=0,\ldots,N_{\mathrm{lag}}-1$, $n=m+\ell$, central set $C_m$, and partner set $P_{ma}$,

$$
S_\ell^{\mathrm{AC}}=\sum_{m=0}^{N-\ell-1}\sum_{a\in C_m}
 w_a^{mn}B_{aa}^{mn},
$$

```math
w_a^{mn}=\begin{cases}
(g_a^m)^3,&\mathrm{time\_dependent\_switch=false},\\
g_a^m(g_a^n)^2,&\mathrm{time\_dependent\_switch=true},
\end{cases}
```

$$
S_\ell^{\mathrm{CC}}=\sum_{m=0}^{N-\ell-1}\sum_{a\in C_m}\sum_{b\in P_{ma}}
 g_a^m(g_b^n)^2B_{ab}^{mn},
$$

$$
C_\ell=\frac{S_\ell^{\mathrm{AC}}+S_\ell^{\mathrm{CC}}}{N-\ell}.
$$

`include_autocorrelation` and `include_cross_correlation` select the respective sums. The time-dependent switch affects only the same-molecule term. Normalization is by origin count, not molecule or pair count. The response also writes the unnormalized numerator and the origin count.

## Fourier transform and SFG spectrum

### Damping

Damping suppresses long-time contributions and broadens spectral features.
`tau = 0` disables damping. 
A positive `tau` multiplies a correlation by

$$
C_k\leftarrow C_k\exp\left(-\frac{k\Delta t}{\tau}\right).
$$

### Reduced complex transform

Only the first two response columns are read. For $t_n=n\Delta t$,

$$
\widehat C_n=[C_n-I_{\mathrm{shift}}C_0]e^{-t_n/\tau},
$$

where $I_{\mathrm{shift}}=1$ for `subtract_initial_response=true`, otherwise zero.
The exponential is omitted when $\tau=0$. 
Values read at $t\gt \mathrm{int\_max}$ are zeroed before this operation.

$$
N_T=\left\lfloor\frac{T}{\Delta t}+10^{-10}\right\rfloor,
$$

$$
A_c(\nu)=\frac12\widehat C_0+\sum_{n=1}^{N_T-2}\widehat C_n\cos(2\pi\nu t_n),
$$

$$
A_s(\nu)=\sum_{n=1}^{N_T-2}\widehat C_n\sin(2\pi\nu t_n).
$$

where $T = $ `spectrum_transform_time_ps`.
The final otherwise eligible sample is excluded. For $K =$ `spectrum_num_wavenumbers`,

$$
\widetilde\nu_k=\frac{k\widetilde\nu_{\max}}{K-1},\qquad
\nu_k=c\widetilde\nu_k,
$$

$$
S_{\mathrm{im}}(\widetilde\nu_k)=\Delta t\,\widetilde\nu_k A_c(\nu_k),
$$

$$
S_{\mathrm{re}}(\widetilde\nu_k)=-\Delta t\,\widetilde\nu_k A_s(\nu_k),
$$

$$
M(\widetilde\nu_k)=S_{\mathrm{im}}^2+S_{\mathrm{re}}^2.
$$

When spectral smoothing is enabled, it is applied separately to the real and imaginary components before the squared magnitude is calculated. 
The fourth output column equals the sum of squares of the reported components. 
Fresnel post-processing reads the complex components and ignores the fourth column.

There is no quantum correction, temperature prefactor, or absolute susceptibility conversion. Increasing the requested grid density does not extend the independent physical frequency range:

$$
\widetilde\nu_{\mathrm{Nyquist}}=\frac1{2c\Delta t}.
$$


### Gaussian smoothing

For Gaussian spectral smoothing with full width at half maximum (FWHM) $w$ and grid spacing $\Delta\widetilde\nu$, the kernel width and length are

$$
\sigma_{\mathrm{grid}}=\frac{w}{2\sqrt{2\ln2}\,\Delta\widetilde\nu},
\qquad K=2\left\lceil3\sigma_{\mathrm{grid}}\right\rceil+1,
$$

$$
g_j=\frac{\exp\left[-(j-(K-1)/2)^2/(2\sigma_{\mathrm{grid}}^2)\right]}
{\sum_{m=0}^{K-1}\exp\left[-(m-(K-1)/2)^2/(2\sigma_{\mathrm{grid}}^2)\right]}.
$$

The Gaussian is truncated at three standard deviations. Only positive-frequency points contribute to the convolution, and the available weights are renormalized at each boundary. The shared smoothing routine sets the zero-frequency output to zero. The unsmoothed SFG signal is already zero at zero frequency. Smoothing a completed spectrum is distinct from damping its time correlation.

### Relation to the literature

The resonant SFG expression in Medders and Paesani [3, Eq. 1] is

$$
\mathrm{Im}\chi_{ijk}^{(2,R)}(\omega)=
\mathrm{Im}\left[\frac{i\omega}{k_BT}\int_0^\infty dt\,
 e^{-i\omega t}\langle\alpha_{ij}(t)\mu_k(0)\rangle\right].
$$

The indices $i$, $j$, and $k$ refer to the sum-frequency, visible, and infrared polarizations, respectively. MB-SPEC does not implement this absolute prefactor or an unrestricted total dipole–polarizability correlation. `SFG.temperature` is required to be positive by the spectrum executable but does not enter its numerical transform.

## Fresnel post-processing

This optional optical step reads reduced complex component spectra and produces
squared amplitudes for SSP, PSS, SPS, and PPP spectra.

### Fresnel JSON options

```jsonc
{
  "xxz_file": "xxz.dat",                // xxz component: wavenumber, imaginary, real
  "zzz_file": "zzz.dat",                // zzz component; supplies frequency grid
  "xzx_file": "xzx.dat",                // xzx component; zxx is assumed equal to xzx
  "visible_incidence_angle_deg": 67.0,  // visible incidence angle in degrees from normal
  "infrared_incidence_angle_deg": 63.0, // IR incidence angle in degrees from normal
  "visible_wavelength_um": 0.8,         // visible wavelength in µm
  "nonresonant_background": 0.0,        // real background added to every component; 0 for none
  "redshift_cm1": 0.0,                  // shift subtracted from wavenumbers in cm⁻¹; 0 for none
  "water_refractive_index_n_file": "fresnel_factors/refractive_index_water_n.txt", // [µm], real index table
  "water_refractive_index_k_file": "fresnel_factors/refractive_index_water_k.txt", // [µm], extinction coefficient
  "ssp_output_file": "sfg_ssp.dat",     // SSP output filename
  "pss_output_file": "sfg_pss.dat",     // PSS output filename
  "sps_output_file": "sfg_sps.dat",     // SPS output filename
  "ppp_output_file": "sfg_ppp.dat"      // PPP output filename
}
```

### Input convention and frequency conversion

For $q\in\{xxz,zzz,xzx\}$,

$$
\widetilde\chi_q(\widetilde\nu)=R_q(\widetilde\nu)+iI_q(\widetilde\nu).
$$

Inputs are `wavenumber`, `imaginary`, `real`; an optional fourth column is ignored. Frequency is taken from `zzz`. 

$$
\widetilde\nu=\widetilde\nu_0-\Delta\widetilde\nu,
\qquad\Delta\widetilde\nu=\mathrm{redshift\_cm1}.
$$

Frequencies below are ordinary frequencies in GHz, with wavelengths in µm:

$$
\nu_{\mathrm{vis}}=\frac{c}{1000\lambda_{\mathrm{vis}}},
$$

$$
\nu_{\mathrm{IR}}=\widetilde\nu\frac{100c}{10^9},\qquad
\lambda_{\mathrm{IR}}=\frac{10^4}{\widetilde\nu},
$$

$$
\nu_{\mathrm{SFG}}=\nu_{\mathrm{vis}}+\nu_{\mathrm{IR}},\qquad
\lambda_{\mathrm{SFG}}=\frac{c}{1000\nu_{\mathrm{SFG}}}.
$$

The same real constant background $B$ is added to each component:

$$
X_q=[R_q+B]+iI_q,\qquad X_{zxx}=X_{xzx}.
$$

The code does not read a separate `zxx` spectrum. 

### Optical model and local-field factors

The incident medium is air, $n_1=1$. Water's complex refractive index is interpolated linearly in wavelength from the configured $n$ and $k$ tables:

$$
n_2(\nu)=n_{\mathrm{water}}(\lambda)+ik_{\mathrm{water}}(\lambda).
$$

Using principal complex roots,

$$
n_i(\nu)=\left[\frac{n_2(\nu)^4+5n_2(\nu)^2}{4n_2(\nu)^2+2}\right]^{1/2},
$$

$$
\gamma=\sin^{-1}\left(\frac{n_1\sin\theta}{n_2}\right).
$$

Angles are measured from the surface normal, with JSON incidence angles converted from degrees to radians.

$$
L_{xx}(\nu,\theta)=\frac{2n_1\cos\gamma}{n_1\cos\gamma+n_2\cos\theta},
$$

$$
L_{yy}(\nu,\theta)=\frac{2n_1\cos\theta}{n_1\cos\theta+n_2\cos\gamma},
$$

$$
L_{zz}(\nu,\theta)=\frac{2n_2\cos\theta\,(n_1/n_i)^2}
{n_1\cos\gamma+n_2\cos\theta}.
$$

These complex factors change both magnitude and phase. For visible incidence $\beta_1$ and infrared incidence $\beta_2$,

$$
\beta=\sin^{-1}\left[\frac{\nu_{\mathrm{vis}}\sin\beta_1+
\nu_{\mathrm{IR}}\sin\beta_2}{\nu_{\mathrm{SFG}}}\right].
$$

$L_{ab}^{(1)}=L_{ab}(\nu_{\mathrm{vis}},\beta_1)$,
$L_{ab}^{(2)}=L_{ab}(\nu_{\mathrm{IR}},\beta_2)$, and
$L_{ab}^{(3)}=L_{ab}(\nu_{\mathrm{SFG}},\beta)$.

### Polarization amplitudes

$$
A_{\mathrm{SSP}}=L_{zz}^{(2)}L_{yy}^{(1)}L_{yy}^{(3)}\sin\beta_2\,X_{xxz},
$$

$$
A_{\mathrm{PSS}}=L_{yy}^{(2)}L_{yy}^{(1)}L_{zz}^{(3)}\sin\beta\,X_{zxx},
$$

$$
A_{\mathrm{SPS}}=L_{yy}^{(2)}L_{zz}^{(1)}L_{yy}^{(3)}\sin\beta_1\,X_{xzx},
$$

$$
\begin{aligned}
A_{\mathrm{PPP}}={}&-L_{zz}^{(2)}L_{xx}^{(1)}L_{xx}^{(3)}
\cos\beta\cos\beta_1\sin\beta_2\,X_{xxz} \\
&-L_{xx}^{(2)}L_{zz}^{(1)}L_{xx}^{(3)}
\cos\beta\sin\beta_1\cos\beta_2\,X_{xzx} \\
&+L_{xx}^{(2)}L_{xx}^{(1)}L_{zz}^{(3)}
\sin\beta\cos\beta_1\cos\beta_2\,X_{zxx} \\
&+L_{zz}^{(2)}L_{zz}^{(1)}L_{zz}^{(3)}
\sin\beta\sin\beta_1\sin\beta_2\,X_{zzz}.
\end{aligned}
$$

$$
I_{\mathrm{pol}}(\widetilde\nu)=|A_{\mathrm{pol}}(\widetilde\nu)|^2.
$$
