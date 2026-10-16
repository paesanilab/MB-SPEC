# MB-SPEC

MB-SPEC calculates vibrational spectra from molecular dynamics trajectories, including those generated with MB-pol and other MB-nrg models. It supports infrared (IR), Raman, sum-frequency generation (SFG), and velocity power spectra. 
For water, MB-SPEC includes built-in many-body dipole moment (MB-$`\mu`$) and polarizability (MB-$`\alpha`$) surfaces that can be applied to compatible trajectories generated with any potential energy surface.

## Build

MB-SPEC requires Autotools, a C++ compiler with OpenMP support, GSL, and Intel
MKL with its LAPACKE headers and libraries. From the repository root:

```bash
./autogen.sh
./configure CXX=c++ CC=cc CXXFLAGS="-O2 -fopenmp" CFLAGS="-O2 -fopenmp"
make -j4
```

`compile.sh` is a convenience script for Intel compilers. Executables are
created in `src/infrared`, `src/raman`, `src/sfg`, and `src/power`.

## Workflows

Use `OMP_NUM_THREADS` to set parallelism for stages that support OpenMP. 
The trajectory may be LAMMPS or `.xyz` format.
All variables for `INFRARED`, `RAMAN`, `SFG`, and `POWER` are specified in the JSON file.

## IR

```bash
./dipole trajectory.lammpstrj mb-spec.json > dipole.dat
./response dipole.dat mb-spec.json > response.dat
./spectrum response.dat mb-spec.json > spectrum.dat
```

Needs `INFRARED` JSON and coordinates; external dipole mode also needs dipoles.

## Raman

```bash
./polarizability trajectory.lammpstrj mb-spec.json > polarizability.dat
./response polarizability.dat mb-spec.json > response.dat
./spectrum response.dat mb-spec.json > spectrum.dat
```

Needs `RAMAN` JSON and coordinates.

## SFG

```bash
./response trajectory.lammpstrj sfg.json > response.dat
./spectrum response.dat sfg.json > spectrum.dat
./fresnel fresnel.json
```

Needs `SFG` JSON, coordinates, and dipoles.

## Power

```bash
./response trajectory.lammpstrj mb-spec.json > response.dat
./spectrum response.dat mb-spec.json > spectrum.dat
```

Needs `POWER` JSON and velocities. 

## References

1. Medders and Paesani, J. Chem. Theory Comput. 9, 4844 (2013)
2. Medders and Paesani, J. Chem. Theory Comput. 11, 1145 (2015)
3. Medders and Paesani, J. Chem. Phys. 142, 212411 (2015)
4. Medders and Paesani, J. Am. Chem. Soc. 138, 3912 (2016)
