# MB-SPEC

MB-SPEC calculates vibrational spectra from molecular dynamics trajectories, including those generated with MB-pol and other MB-nrg models. It supports infrared (IR), Raman, sum-frequency generation (SFG), and velocity power spectra. 
For water, MB-SPEC includes built-in many-body dipole moment (MB-$`\mu`$) and polarizability (MB-$`\alpha`$) surfaces that can be applied to compatible trajectories generated with any potential energy surface.

`Release status`: MB-SPEC is currently in beta. We plan to improve the build process in an upcoming release.

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
In addition, `src/tools/README-ind-dip-pol.md` describes the TTM4-Ind model and 
how MB-SPEC uses it to calculate induction contributions to the dipole moment and polarizability of water.

Some LAMMPS and .xyz input file templates can be found at `/src/unittests/inputs`.

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

## Scientific Limitation

Use only the `1b` and `1b+2b` polarizability surfaces for periodic bulk simulations. 
The `nb`, `1b+nb`, and `1b+2b+nb` modes currently lack Ewald summation, and have not been validated 
for periodic bulk simulations. We aim to address this limitation in a future release.

## References

The following publications describe the theoretical foundations of MB-SPEC and selected applications in vibrational spectroscopy.

1. Medders and Paesani, J. Chem. Theory Comput. 9, 4844 (2013)
2. Medders and Paesani, J. Chem. Theory Comput. 11, 1145 (2015)
3. Medders and Paesani, J. Chem. Phys. 142, 212411 (2015)
4. Medders and Paesani, J. Am. Chem. Soc. 138, 3912 (2016)
5. Straight and Paesani, J. Phys. Chem. B 120, 8359 (2016)
6. Moberg, Straight, Knight, and Paesani, J. Phys. Chem. Lett. 8, 2579 (2017) 
7. Reddy, Moberg, Straight, and Paesani, J. Chem. Phys. 147, 244504 (2017)
8. Moberg, Straight, and Paesani, J. Phys. Chem. B 122, 4356 (2018)
9. Moberg, Sharp, and Paesani, J. Phys. Chem. B 122, 10572 (2018)
10. Hunter, Shakib, and Paesani, J. Phys. Chem. B 122, 10754 (2018)
11. Reddy, Thiraux, Rudd, Lin, Adel, Joutsuka, Geiger, Allen, Morita, and Paesani, Chem 4, 1629 (2018)
12. Sengupta, Moberg, Paesani, and Tyrode, J. Phys. Chem. Lett. 9, 6744 (2018)
13. Ohno, Wang, Paesani, Skinner, and Geiger, J. Phys. Chem. A 122, 4457 (2018)
14. Sun, Tang, Imoto, Moberg, Ohto, Paesani, Bonn, Backus, and Nagata, Phys. Rev. Lett. 121, 246101 (2018)
15. Hudait, Moberg, Qiu, Odendahla, Paesani, and Molinero, Proc. Natl. Acad. Sci. U.S.A. 115, 8266 (2018)
16. Moberg, Li, Reddy, and Paesani, J. Chem. Phys. 150, 034701 (2019)
17. Moberg, Becker, Dierking, Zurheide, Bandow, Buck, Hudait, Molinero, Paesania, and Zeuch, Proc. Natl. Acad. Sci. U.S.A. 116, 24413 (2019)
18. Rashmi and Paesani, J. Am. Chem. Soc. 147, 1903 (2025)
19. Rashmi, Balogun, Azom, Agnew, Kumar, and Paesani, ACS Nano 19, 4876 (2025)
20. Agnew, Savoj, Rashmi, Savala, and Paesani, J. Phys. Chem. B 129, 6584 (2025)
21. Rashmi, Savala, Agnew, Savoj, and Paesani, J. Phys. Chem. B 129, 7301 (2025)
22. Dick, Rashmi, Agnew, Saha, Zhu, Kirchner, and Paesani, ACS Nano 20, 1304 (2026)
23. Balogun, Rashmi, Azom, Kumar, and Paesani, J. Phys. Chem. Lett. 17, 2464 (2026)
24. Rashmi, Chiang, Nagata, Bonn, and Paesani, J. Phys. Chem. Lett. 17 3013 (2026)
25. Solomon, Rashmi, Zhou, and Paesani, J. Phys. Chem. A 130, 6933 (2026)
