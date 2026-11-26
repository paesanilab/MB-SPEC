# Power spectra with MB-SPEC

POWER computes the autocorrelation of the total atomic velocity and transforms
it into a normalized spectrum. This collective observable includes cross terms
between atoms; it differs from a sum of per-atom velocity autocorrelations.

Copy [the power template](../../json_template/power.json) to your working directory
as `power.json`. Set the frame interval and skip count for your trajectory before
running the two stages.

```bash
mbspec-power-response trajectory.lammpstrj power.json > response.dat
mbspec-power-spectrum response.dat power.json > spectrum.dat
```

## Units

| Quantity | Unit / convention |
|---|---|
| Velocities | Input units are retained |
| Velocity autocorrelation | Squared input velocity units |
| Time | ps; transform converts time to seconds |
| XYZ cell lengths | Å |
| Frequency $\nu$, angular frequency $\omega$ | Hz and rad/s |
| Wavenumber $\widetilde\nu$, smoothing FWHM | cm⁻¹; $\nu=c\widetilde\nu$ with $c$ in cm/s |
| Intensity | Implemented normalized line shape; not established as an absolute density of states |

Uniform velocity rescaling cancels from the spectrum's normalization by $C(0)$.
It does not cancel from the response file.

## Outputs

| File | Columns, in order |
|---|---|
| `response.dat` | Lag time [ps], collective velocity autocorrelation [squared input velocity units] |
| `spectrum.dat` | Wavenumber [cm⁻¹], normalized intensity |


## JSON options

Put all settings inside `POWER`.
**Remove the `//` comments before saving this JSONC example as an input file.**

```jsonc
{
  "POWER": {
    "box": [0.00, 19.66, 0.00, 19.66, 0.00, 19.66], // Required for XYZ only, in Å; unused by the velocity observable
    "imcon": 2,                                  // XYZ cell convention: 1 cubic, 2 orthorhombic (default), 3 triclinic
    "dump_frequency": 0.004,                     // Required frame interval in ps
    "trajectory_duration_ps": 0.0,               // Duration read in ps; 0 reads all frames (default)
    "int_max": 10.0,                             // Required exclusive correlation lag limit in ps
    "skip_initial_atoms": 1,                     // Omit leading atom records; default 0; example omits the first atom
    "subtract_mean_velocity": false,             // Remove mean total velocity before correlation; default false
    "spectrum_transform_time_ps": 10.0,          // Required response duration transformed, in ps
    "tau": 0.0,                                  // Exponential damping time in ps; 0 disables it (default)
    "smoothing_fwhm_cm1": 15.0                   // Gaussian spectral FWHM in cm⁻¹; omitted or 0 disables it
  }
}
```

For XYZ, `box` uses six bounds `[xmin,xmax,ymin,ymax,zmin,zmax]` or nine matrix
elements `[xx,xy,xz,yx,yy,yz,zx,zy,zz]` for a triclinic cell. LAMMPS input does
not use the JSON box. `imcon` does not introduce periodic corrections to the
velocity correlation.

## Trajectory input

The power response reads LAMMPS atom-column names and requires native `vx`,
`vy`, and `vz` velocity columns. Other atom fields are ignored. The LAMMPS box
and timestep records must be present.

An XYZ velocity trajectory is selected only when its filename ends in lowercase
`.xyz`; an uppercase `.XYZ` suffix is not recognized as XYZ input. Each frame
has the standard atom-count and comment lines, followed by one row per atom:

```text
label vx vy vz
```

Frame times are frame index times `dump_frequency`. LAMMPS timestep records
and XYZ comments do not supply the time axis.

## Velocity autocorrelation

$$
\mathbf V(t_i)=\sum_{a\;\mathrm{included}}\mathbf v_a(t_i),
$$

$$
C(k\Delta t)=\frac1{N-k}\sum_{i=0}^{N-k-1}
\mathbf V(t_i)\cdot\mathbf V(t_{i+k}).
$$

This unweighted collective correlation includes cross terms between atoms. It differs from $\sum_a\langle\mathbf v_a(0)\cdot\mathbf v_a(t)\rangle$. With `subtract_mean_velocity`, the trajectory mean of $\mathbf V$ is removed first.

## Fourier transform and power spectrum

### Damping

A positive `tau` multiplies a correlation by

$$
C_k\leftarrow C_k\exp\left(-\frac{k\Delta t}{\tau}\right).
$$

`tau = 0` disables damping. Damping suppresses long-time contributions and broadens spectral features.

### Transform and normalization

The spectrum reader retains rows with $t\lt \mathrm{int\_max}$ and transforms the first $N_T=\lfloor T_{\mathrm{transform}}/\Delta t_{\mathrm{ps}}+10^{-10}\rfloor$ samples. Here $T_{\mathrm{transform}}$ is the JSON parameter `spectrum_transform_time_ps`. With $\Delta t=10^{-12}\Delta t_{\mathrm{ps}}$ in seconds and $c=2.99792458\times10^{10}$ cm/s, the cosine sum is

$$
D(f_n)=C_0+2\sum_{k=1}^{N_T-1}C_k
\cos\left(\frac{2\pi nk}{2N_T-1}\right),\qquad
f_n=\frac{n}{(2N_T-1)\Delta t},
$$

$$
I_{\mathrm{code}}(\widetilde\nu_n)=\frac{2cD(c\widetilde\nu_n)}{C(0)}.
$$

The discrete sum has no $\Delta t$ factor. The output is therefore an implemented normalized line shape whose spectral-density normalization requires analytical validation; it is not established as an absolute vibrational density of states. The implementation requires $C(0)>0$ and has no temperature, volume, dielectric, Planck-constant, or quantum-correction prefactor. 

The dimensionless $10^{-10}$ in the sample count prevents floating-point
rounding just below an integer from dropping one sample.

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

The Gaussian is truncated at three standard deviations. 
