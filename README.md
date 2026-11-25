# MB-SPEC

MB-SPEC calculates vibrational spectra from molecular dynamics trajectories, including those generated with MB-pol and other MB-nrg models. It supports infrared (IR), Raman, sum-frequency generation (SFG), and velocity power spectra. 
For water, MB-SPEC includes built-in many-body dipole moment (MB-$`\mu`$) and polarizability (MB-$`\alpha`$) surfaces that can be applied to compatible trajectories generated with any potential energy surface.

## Build

MB-SPEC requires a C++17 or newer compiler with OpenMP support and a LAPACKE
implementation with LP64 (32-bit integer) headers and libraries. The default
backend is generic LAPACKE, with Intel MKL available as an option.
`mbspec-ir-dipole` is the only installed program that links LAPACKE.
Bootstrapping from a Git checkout requires Autoconf 2.70+ and Automake 1.16+.
The required macros from the GNU Autoconf Archive and pkgconf are bundled in
[`m4/`](m4/README.md). Release tarballs include the generated `configure` script.

For a generic LAPACKE installation on the compiler's normal search paths, run
these commands from the source directory.

```bash
autoreconf -fi  # Git checkouts only
./configure
make -j4
make check
make install
```

`./configure` detects a C++ compiler and uses Autoconf's default flags, typically
`-g -O2` with GCC. Set `CXX` to choose a compiler and `CXXFLAGS` to supply your
own flags. For example, `./configure CXX=g++ CXXFLAGS="-O2"` selects GCC and
sets optimization level 2.

The default `--with-lapack=generic` uses `<lapacke.h>`. It first tries the
flags from the `lapacke` pkg-config module, when pkg-config is installed (it is
optional), then `-llapacke`, then `-llapacke -llapack -lblas`. Supply `LAPACKE_CFLAGS` and `LAPACKE_LIBS`
for a custom installation or library combination, for example:

```bash
./configure \
    LAPACKE_CFLAGS="-I/path/to/lapack/include" \
    LAPACKE_LIBS="-L/path/to/lapack/lib -llapacke -llapack -lblas"
```

Explicit `LAPACKE_LIBS`, including an empty value, are validated and used
without fallback. OpenBLAS can provide LAPACKE itself or serve as the backend
for a separate LAPACKE library. Its build can omit LAPACKE, so OpenBLAS alone
does not guarantee that the required header and routines are available.

Select Intel MKL explicitly with `--with-lapack=mkl`, which uses
`<mkl_lapacke.h>`. After initializing your MKL environment (so `MKLROOT` is
set), a Linux build using the sequential LP64 libraries can be configured as:

```bash
./configure --with-lapack=mkl
make -j4
make check
make install
```

To build with the Intel oneAPI C++ compiler, initialize the oneAPI environment
and run `./configure CXX=icpx --with-lapack=mkl` before the same `make` steps.

Without `LAPACKE_LIBS`, an MKL build tries, in order:

1. the compiler's own MKL support (for example `CXXFLAGS="-O2 -qmkl=sequential"`
   with `icpx`)
2. oneAPI's `mkl-dynamic-lp64-seq` pkg-config module
3. `-I${MKLROOT}/include` with the sequential LP64 libraries under `MKLROOT`

The MKL libraries must be on the runtime library path (initializing the MKL
environment sets `LD_LIBRARY_PATH`).
For other threading layers or layouts, pass the flags explicitly:

```bash
./configure --with-lapack=mkl \
    LAPACKE_CFLAGS="-I${MKLROOT}/include" \
    LAPACKE_LIBS="-L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl"
```

For either implementation, the headers and libraries must agree on
the LP64 integer interface. Native configuration runs a small numerical probe
to check the library ABI. Make the selected shared libraries available at
configure time, for example by initializing your MKL environment. A cross-build
checks linkage during configuration and requires `make check` on the target
to verify the runtime ABI. `make check` also tests singular-matrix failure.

Configuration enables C++17 using compiler-specific flags and rejects modes
older than C++17. To select a newer mode, put its flag in `CXXFLAGS`, for example
`CXXFLAGS="-O2 -std=c++20"`. User `CXXFLAGS` are preserved. OpenMP flags are
detected automatically; disabling OpenMP is unsupported because the programs
call its runtime directly.

With GCC-compatible compilers (GCC, Clang, `icpx`), `-Wall` is added
automatically. Two options are meant for development:

- `--enable-werror` treats compiler warnings as errors.
- `--enable-sanitizers` builds with AddressSanitizer and
  UndefinedBehaviorSanitizer (`-fsanitize=address,undefined`). Any finding
  aborts the program, so `make check` fails. The generated polynomial library
  is not instrumented (see `src/tools/README.md`).

The old `MKL_CFLAGS` and `MKL_LIBS` variables are ignored with a warning. Rename them to
`LAPACKE_CFLAGS` and `LAPACKE_LIBS` and add `--with-lapack=mkl` to retain an
MKL build. The legacy `--enable-mkl`, `--disable-mkl`, and `--with-mkl-*`
options remain removed. LAPACKE is required and cannot be disabled.

Use `make -jN` to select the number of parallel build jobs. Run `make clean`
before changing compilers or compilation flags, then rerun `./configure`
with the new settings, rebuild, and rerun `make install` to update `./bin/`. Run
`autoreconf -fi` again after editing `configure.ac`, `Makefile.am`, or files in
`m4/`.

By default, `make install` puts the programs in `./bin/` and their data in
`./share/mb-spec/` relative to the directory where you ran `configure`.
Use `./configure --prefix=/absolute/path` to install into another directory.
Run `make uninstall` from the configured build directory to remove the installed files.
Each program is named `mbspec-<workflow>-<stage>`.

| Workflow | Programs |
|----------|----------|
| IR | `mbspec-ir-dipole`, `mbspec-ir-response`, `mbspec-ir-spectrum` |
| Raman | `mbspec-raman-polarizability`, `mbspec-raman-response`, `mbspec-raman-spectrum` |
| SFG | `mbspec-sfg-response`, `mbspec-sfg-spectrum`, `mbspec-sfg-fresnel` |
| Power | `mbspec-power-response`, `mbspec-power-spectrum` |

Without installing, the programs are in `src/infrared`, `src/raman`, `src/sfg`
and `src/power` of the build tree. The JSON templates are installed to
`$prefix/share/mb-spec/json_template/`. The example Fresnel configuration and
the water refractive-index tables are installed to `$prefix/share/mb-spec/sfg/`.
`fresnel.json` refers to the tables by relative path, so copy
`fresnel.json` and `fresnel_factors/` into your working directory.

## Testing

Run `make check` from the top-level build directory. It runs the matrix and
2B polynomial unit tests and the [trajectory regression suites](src/unittests/README.md).
The trajectory suites check numerical references, analytical cases, input formats,
threading, and complete workflows against a pinned `master` baseline.

## Quickstart

Run these commands after `make install` to calculate an IR spectrum from the
bundled eight-frame water trajectory. Set `source_dir` to the source root when
building in a separate directory.

```bash
source_dir=.
example="$source_dir/src/unittests/historical/inputs"
OMP_NUM_THREADS=1 ./bin/mbspec-ir-dipole "$example/water4.xyz" "$example/ir.json" > dipole.dat
./bin/mbspec-ir-response dipole.dat "$example/ir.json" > response.dat
./bin/mbspec-ir-spectrum response.dat "$example/ir.json" > spectrum.dat
```

`dipole.dat` contains eight rows with time in ps and three dipole components in
Debye. `response.dat` contains six rows with lag in ps and the autocorrelation
in (C·m)². `spectrum.dat` contains six rows with wavenumber and IR intensity,
both in cm⁻¹.

## Workflows

Run the commands below from the top-level build directory after `make install`.
Copy the corresponding JSON file from `./share/mb-spec/json_template/` and
adjust it for your trajectory. To run the programs from another directory, add
the installed `bin/` directory to your `PATH`. Use `OMP_NUM_THREADS` to set
parallelism for stages that support OpenMP.

Use `--help` on any executable to see its arguments. A workflow configuration
must contain the section for the executable you are running. Use the matching
template for recognized keys. A single file can contain several workflow
sections. IR, Raman, and Power spectrum inputs require two finite numeric
columns and an evenly spaced time grid.

## IR

```bash
./bin/mbspec-ir-dipole trajectory.lammpstrj ir.json > dipole.dat
./bin/mbspec-ir-response dipole.dat ir.json > response.dat
./bin/mbspec-ir-spectrum response.dat ir.json > spectrum.dat
```

IR takes LAMMPS coordinates and dipole columns in external dipole mode. The
internal water dipole model (`DMS=true`) accepts LAMMPS coordinates or XYZ.
See the [IR guide](src/infrared/README.md).

## Raman

```bash
./bin/mbspec-raman-polarizability trajectory.lammpstrj raman.json > polarizability.dat
./bin/mbspec-raman-response polarizability.dat raman.json > response.dat
./bin/mbspec-raman-spectrum response.dat raman.json > spectrum.dat
```

The water polarizability stage accepts coordinates in LAMMPS or XYZ format.
See the [Raman guide](src/raman/README.md).

## SFG

```bash
./bin/mbspec-sfg-response trajectory.lammpstrj sfg.json > response.dat
./bin/mbspec-sfg-spectrum response.dat sfg.json > spectrum.dat
```

SFG takes a LAMMPS water slab trajectory with coordinate and dipole columns.
The optional Fresnel step uses the `xxz`, `zzz`, and `xzx` component spectra
and the water refractive-index tables. See [Fresnel post-processing](src/sfg/README.md#fresnel-post-processing).

## Power

```bash
./bin/mbspec-power-response trajectory.lammpstrj power.json > response.dat
./bin/mbspec-power-spectrum response.dat power.json > spectrum.dat
```

Power takes LAMMPS dump files with `vx`, `vy`, and `vz` columns or XYZ files
whose atom rows contain velocities. See the [power spectrum guide](src/power/README.md).

## References

1. Medders and Paesani, J. Chem. Theory Comput. 9, 4844 (2013)
2. Medders and Paesani, J. Chem. Theory Comput. 11, 1145 (2015)
3. Medders and Paesani, J. Chem. Phys. 142, 212411 (2015)
4. Medders and Paesani, J. Am. Chem. Soc. 138, 3912 (2016)
