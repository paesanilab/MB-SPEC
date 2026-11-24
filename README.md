# MB-SPEC

MB-SPEC calculates vibrational spectra from molecular dynamics trajectories, including those generated with MB-pol and other MB-nrg models. It supports infrared (IR), Raman, sum-frequency generation (SFG), and velocity power spectra. 
For water, MB-SPEC includes built-in many-body dipole moment (MB-$`\mu`$) and polarizability (MB-$`\alpha`$) surfaces that can be applied to compatible trajectories generated with any potential energy surface.

## Build

MB-SPEC requires a C++17 or newer compiler with OpenMP support and a LAPACKE
implementation with LP64 (32-bit integer) headers and libraries. Generic
LAPACKE is the default; Intel MKL is optional. No C compiler, GSL or other
libraries are required. Only `mbspec-ir-dipole` links LAPACKE. Bootstrapping from a Git
checkout requires Autoconf 2.70+ and Automake 1.16+. The Autoconf macros it
uses (from the GNU Autoconf Archive and pkgconf) are included in `m4/`, so no
separate Autoconf Archive or pkg-config macro installation is
needed. The vendored macros' upstream revisions are recorded in `m4/README.md`. Release tarballs
include the generated `configure` script.

For a generic LAPACKE installation on the compiler's normal search paths:

```bash
./autogen.sh                 # Git checkouts only
./configure CXX=c++ CXXFLAGS="-O2"
make -j4
make check
```

The default `--with-lapack=generic` uses `<lapacke.h>`. It first tries the
flags from the `lapacke` pkg-config module, when pkg-config is installed (it is
optional), then `-llapacke`, then `-llapacke -llapack -lblas`. Supply `LAPACKE_CFLAGS` and `LAPACKE_LIBS`
for a custom installation or library combination, for example:

```bash
./configure CXX=c++ CXXFLAGS="-O2" \
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
./configure --with-lapack=mkl CXX=c++ CXXFLAGS="-O2"
make -j4
make check
```

Without `LAPACKE_LIBS`, an MKL build tries, in order:

1. the compiler's own MKL support (for example `CXXFLAGS="-O2 -qmkl=sequential"`
   with `icpx`)
2. oneAPI's `mkl-dynamic-lp64-seq` pkg-config module
3. `-I${MKLROOT}/include` with the sequential LP64 libraries under `MKLROOT`

The MKL libraries must be on the runtime library path (initializing the MKL
environment sets `LD_LIBRARY_PATH`).
For other threading layers or layouts, pass the flags explicitly:

```bash
./configure --with-lapack=mkl CXX=c++ CXXFLAGS="-O2" \
    LAPACKE_CFLAGS="-I${MKLROOT}/include" \
    LAPACKE_LIBS="-L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl"
```

For either implementation, the headers and libraries must agree on
the LP64 integer interface; ILP64 is unsupported. Configuration checks the
header and compilation/linkage of all three required LAPACKE routines. It
does not run numerical tests. `make check` tests inversion, singular-matrix
failure and symmetric eigenvalues against the selected implementation.

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

`make check` runs unit tests of the matrix routines and the 2B polynomial,
and end-to-end tests of every workflow. The workflow tests compare each
stage's output with reference outputs from the `master` branch; see
`tests/README.md`. The regression suites in `src/unittests/` also check numerical
references, analytical cases, and agreement across input formats and thread counts.

The old `MKL_CFLAGS` and `MKL_LIBS` variables are ignored with a warning. Rename them to
`LAPACKE_CFLAGS` and `LAPACKE_LIBS` and add `--with-lapack=mkl` to retain an
MKL build. The legacy `--enable-mkl`, `--disable-mkl`, and `--with-mkl-*`
options remain removed. LAPACKE is required and cannot be disabled.

Use `make -jN` to select the number of parallel build jobs. Run `make clean`
before changing compilers or compilation flags, then rerun `./configure`
with the new settings and rebuild. Run
`./autogen.sh` again after editing `configure.ac`, `Makefile.am`, or files in
`m4/`.

`make install` installs the programs to `$prefix/bin` (default `/usr/local`;
choose another with `./configure --prefix=DIR`). Each is named
`mbspec-<workflow>-<stage>`:

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

The workflow commands below assume the programs are on your `PATH`.

## Workflows

Use `OMP_NUM_THREADS` to set parallelism for stages that support OpenMP. 
The trajectory may be LAMMPS or `.xyz` format.
All variables for `INFRARED`, `RAMAN`, `SFG`, and `POWER` are specified in the JSON file.

## IR

```bash
mbspec-ir-dipole trajectory.lammpstrj ir.json > dipole.dat
mbspec-ir-response dipole.dat ir.json > response.dat
mbspec-ir-spectrum response.dat ir.json > spectrum.dat
```

Needs `INFRARED` JSON and coordinates; external dipole mode also needs dipoles.

## Raman

```bash
mbspec-raman-polarizability trajectory.lammpstrj raman.json > polarizability.dat
mbspec-raman-response polarizability.dat raman.json > response.dat
mbspec-raman-spectrum response.dat raman.json > spectrum.dat
```

Needs `RAMAN` JSON and coordinates.

## SFG

```bash
mbspec-sfg-response trajectory.lammpstrj sfg.json > response.dat
mbspec-sfg-spectrum response.dat sfg.json > spectrum.dat
mbspec-sfg-fresnel fresnel.json
```

Needs `SFG` JSON, coordinates, and dipoles.

The fresnel factor calculation also requires the three SFG tensor component files (xxz.dat, zzz.dat, and xzx.dat) 
and the refractive index files in `src/sfg/fresnel_factors/`. 

## Power

```bash
mbspec-power-response trajectory.lammpstrj power.json > response.dat
mbspec-power-spectrum response.dat power.json > spectrum.dat
```

Needs `POWER` JSON and velocities. 

## References

1. Medders and Paesani, J. Chem. Theory Comput. 9, 4844 (2013)
2. Medders and Paesani, J. Chem. Theory Comput. 11, 1145 (2015)
3. Medders and Paesani, J. Chem. Phys. 142, 212411 (2015)
4. Medders and Paesani, J. Am. Chem. Soc. 138, 3912 (2016)
