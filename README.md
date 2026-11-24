# MB-SPEC

MB-SPEC calculates vibrational spectra from molecular dynamics trajectories, including those generated with MB-pol and other MB-nrg models. It supports infrared (IR), Raman, sum-frequency generation (SFG), and velocity power spectra. 
For water, MB-SPEC includes built-in many-body dipole moment (MB-$`\mu`$) and polarizability (MB-$`\alpha`$) surfaces that can be applied to compatible trajectories generated with any potential energy surface.

## Build

MB-SPEC requires a C++17 or newer compiler with OpenMP support and Intel MKL
with its LAPACKE headers and LP64 libraries. A C compiler is also used by
Automake's archiver detection. GSL is not required. Bootstrapping from a Git
checkout requires Autoconf 2.70+ and Automake 1.16+. The GNU Autoconf Archive
C++ standard check is included in `m4/`; no separate Archive installation is
needed. Its upstream revision is recorded in `m4/README.md`. Release tarballs
include the generated `configure` script.

Supply `MKL_CFLAGS` and `MKL_LIBS` for your installation and threading layer.
Configuration does not assume installation paths, architecture flags, or an
OpenMP runtime. After initializing your MKL environment, a Linux build using
the sequential LP64 libraries can be configured from the repository root as:

```bash
./autogen.sh                 # Git checkouts only
./configure CXX=c++ CXXFLAGS="-O2" \
    MKL_CFLAGS="-I${MKLROOT}/include" \
    MKL_LIBS="-L${MKLROOT}/lib -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
make -j4
```

Configuration enables C++17 using compiler-specific flags and rejects modes
older than C++17. To select a newer mode, put its flag in `CXXFLAGS`, for example
`CXXFLAGS="-O2 -std=c++20"`. User `CXXFLAGS` are preserved. OpenMP flags are
detected automatically; disabling OpenMP is unsupported because the programs
call its runtime directly.

The legacy `--enable-mkl`, `--disable-mkl`, and `--with-mkl-*` options have been
removed. MKL remains required; use `MKL_CFLAGS` and `MKL_LIBS` to configure it.

`compile.sh` is a convenience script using `icpx` by default. Executables are
created in `src/infrared`, `src/raman`, `src/sfg`, and `src/power`.

## Workflows

Use `OMP_NUM_THREADS` to set parallelism for stages that support OpenMP. 
The trajectory may be LAMMPS or `.xyz` format.
All variables for `INFRARED`, `RAMAN`, `SFG`, and `POWER` are specified in the JSON file.

## IR

```bash
./dipole trajectory.lammpstrj ir.json > dipole.dat
./response dipole.dat ir.json > response.dat
./spectrum response.dat ir.json > spectrum.dat
```

Needs `INFRARED` JSON and coordinates; external dipole mode also needs dipoles.

## Raman

```bash
./polarizability trajectory.lammpstrj raman.json > polarizability.dat
./response polarizability.dat raman.json > response.dat
./spectrum response.dat raman.json > spectrum.dat
```

Needs `RAMAN` JSON and coordinates.

## SFG

```bash
./response trajectory.lammpstrj sfg.json > response.dat
./spectrum response.dat sfg.json > spectrum.dat
./fresnel fresnel.json
```

Needs `SFG` JSON, coordinates, and dipoles.

The fresnel factor calculation also requires the three SFG tensor component files (xxz.dat, zzz.dat, and xzx.dat) 
and the refractive index files in `src/sfg/fresnel_factors/`. 

## Power

```bash
./response trajectory.lammpstrj power.json > response.dat
./spectrum response.dat power.json > spectrum.dat
```

Needs `POWER` JSON and velocities. 

## References

1. Medders and Paesani, J. Chem. Theory Comput. 9, 4844 (2013)
2. Medders and Paesani, J. Chem. Theory Comput. 11, 1145 (2015)
3. Medders and Paesani, J. Chem. Phys. 142, 212411 (2015)
4. Medders and Paesani, J. Am. Chem. Soc. 138, 3912 (2016)
