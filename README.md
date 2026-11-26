# MB-SPEC

MB-SPEC calculates vibrational spectra from molecular dynamics trajectories, including those generated with MB-pol and other MB-nrg models. It supports infrared (IR), Raman, sum-frequency generation (SFG), and velocity power spectra. 
For water, MB-SPEC includes built-in many-body dipole moment (MB-$`\mu`$) and polarizability (MB-$`\alpha`$) surfaces that can be applied to compatible trajectories generated with any potential energy surface.

## Build

MB-SPEC requires Autotools, a C++ compiler with OpenMP support, GSL, and Intel
MKL with its LAPACKE headers and libraries. From the repository root, run `compile.sh` or

```bash
autoreconf -fi
./configure CC=icx CXX=icpx CFLAGS="-g -O2 -fopenmp -Wall -std=c99" CXXFLAGS="-g -O2 -fopenmp -Wall"
make -j4
```

All executables are created in `src/infrared`, `src/raman`, `src/sfg`, and `src/power`.

## Workflows

Use `OMP_NUM_THREADS` to set parallelism for stages that support OpenMP. 
The trajectory may be LAMMPS or `.xyz` format.
All variables for `INFRARED`, `RAMAN`, `SFG`, and `POWER` are specified in their respective 
JSON files in `json_template/`

Each module’s source directory contains a README that describes its variables and equations in detail (e.g., `src/infrared/README.md`).
In addition, `src/README-ind-dip-pol.md` describes the TTM4-Ind model and 
how MB-SPEC uses it to calculate induction contributions to the dipole moment and polarizability of water.

## IR

```bash
./src/ir/dipole trajectory.lammpstrj json_template/ir.json > dipole.dat
./src/ir/response dipole.dat json_template/ir.json > response.dat
./src/ir/spectrum response.dat json_template/ir.json > spectrum.dat
```

Needs coordinates (and dipoles if external).

## Raman

```bash
./src/raman/polarizability trajectory.lammpstrj json_template/raman.json > polarizability.dat
./src/raman/response polarizability.dat json_template/raman.json > response.dat
./src/raman/spectrum response.dat json_template/raman.json > spectrum.dat
```

Needs coordinates.

## SFG

```bash
./src/sfg/response trajectory.lammpstrj json_template/sfg.json > response.dat
./src/sfg/spectrum response.dat json_template/sfg.json > spectrum.dat
./src/sfg/fresnel src/sfg/fresnel.json
```

Needs coordinates (and dipoles if external).
The fresnel factor calculation also requires the three SFG tensor component files (xxz.dat, zzz.dat, and xzx.dat) 
and the refractive index files in `src/sfg/fresnel_factors/`. 

## Power

```bash
./src/power/response trajectory.lammpstrj json_template/power.json > response.dat
./src/power/spectrum response.dat json_template/power.json > spectrum.dat
```

Needs velocities. 

## References

1. Medders and Paesani, J. Chem. Theory Comput. 9, 4844 (2013)
2. Medders and Paesani, J. Chem. Theory Comput. 11, 1145 (2015)
3. Medders and Paesani, J. Chem. Phys. 142, 212411 (2015)
4. Medders and Paesani, J. Am. Chem. Soc. 138, 3912 (2016)
5. Straight and Paesani, J. Phys. Chem. B 120, 8359 (2016)
6. Moberg, Straight, Knight, and Paesani, J. Phys. Chem. Lett. 8, 2579 (2017) 
7. Moberg, Straight, and Paesani, J. Phys. Chem. B 122, 4356 (2018)
8. Hunter, Shakib, and Paesani, J. Phys. Chem. B 122, 10754 (2018)
9. Reddy, ..., and Paesani, Chem 4, 1629 (2018)
10. Moberg, Li, Reddy, and Paesani, J. Chem. Phys. 150, 034701 (2019)
