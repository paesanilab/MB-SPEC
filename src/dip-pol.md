# TTM4 dipole moment and polarizability

TTM4 represents water electrostatics using geometry-dependent permanent charges
and mutually induced dipoles. The dipole moment describes the system's charge
separation; the polarizability describes how that dipole changes under an
applied electric field.

In MB-SPEC, the relevant implementation is the TTM4-based induction model. The
MB-mu/MB-alpha paper calls this **TTM4-Ind**, distinguishing it from the complete
TTM4-F water potential. See [Medders and Paesani, Section 2](https://doi.org/10.1021/ct501131j).

## Dipole calculation

Each water molecule has four electrostatic sites: O, H, H, and a virtual M site
near the oxygen. The code obtains geometry-dependent charges from the
Partridge–Schwenke dipole surface and redistributes them onto the H and M sites.

The permanent charges produce an electric field at the polarizable sites.
Each site develops an induced dipole:

$$
\boldsymbol{\mu}_{i}^{\mathrm{ind}}
=
a_i\left[
\mathbf{E}_{i}^{\mathrm{perm}}
+\sum_{j\ne i}\mathbf{T}_{ij}^{\mathrm{damped}}
\,\boldsymbol{\mu}_{j}^{\mathrm{ind}}
\right].
$$

Here, $a_i$ is the site polarizability, $\mathbf E_i^{\mathrm{perm}}$ is the
field from the permanent charges, and $\mathbf T_{ij}^{\mathrm{damped}}$
describes the damped dipole interaction.

The dipoles are coupled: each induced dipole contributes to the field inducing
the others. MB-SPEC solves these simultaneous linear equations using a
Cholesky decomposition.

For the polarizable O and H sites, the damped interaction tensor is

$$
\mathbf T_{ij}^{\mathrm{damped}}
=\frac{3\left[1-\left(1+\frac43x_{ij}\right)e^{-x_{ij}}\right]}
{r_{ij}^{5}}\,
\mathbf r_{ij}\mathbf r_{ij}^{\mathsf T}
-\frac{1-e^{-x_{ij}}}{r_{ij}^{3}}\mathbf I,
$$

with

$$
x_{ij}=d_{ij}
\left[\frac{r_{ij}}{(a_i a_j)^{1/6}}\right]^4.
$$

Here, $\mathbf r_{ij}=\mathbf r_i-\mathbf r_j$,
$r_{ij}=|\mathbf r_{ij}|$, $\mathbf I$ is the identity tensor, and $d_{ij}$
is the dimensionless damping parameter. The implementation uses
$d_{ij}=0.055$ for intermolecular pairs and intramolecular H–H pairs, and
$d_{ij}=0.626$ for intramolecular O–H pairs. The M site has zero polarizability
and carries no induced dipole.

At large separation, the exponential terms vanish and the tensor reduces to
the usual point-dipole interaction. At short separation, damping suppresses
the singular behavior. These expressions use the model's electrostatic units;
with distances in Å and site polarizabilities in Å³, the tensor has units of
Å⁻³. See `smear2` in [ttm4-smear.cpp](tools/ttm4-smear.cpp) and the tensor
construction in [ttm4-dip.cpp](tools/ttm4-dip.cpp).

The total dipole is

$$
\mathbf M^{\mathrm{TTM4}}
=\underbrace{\sum_i q_i\mathbf r_i}_{\text{permanent}}
+\underbrace{\sum_i\boldsymbol\mu_i^{\mathrm{ind}}}_{\text{induced}}.
$$

The distinction matters in the code: `ttm4_es::dip` returns permanent plus
induced dipoles, whereas `dipind` and the separate `ttm4_dip` routine return
only the induced contribution. Dipole outputs are in Debye. See
[ttm4-es.cpp](tools/ttm4-es.cpp) and [ttm4-dip.cpp](tools/ttm4-dip.cpp).

## Polarizability calculation

At fixed atomic positions, a small uniform external field $\mathbf F$ changes
the induced dipoles throughout the system. The system polarizability tensor
measures this response:

$$
A_{ab}^{\mathrm{TTM4}}
=\left.\frac{\partial M_a}{\partial F_b}\right|_{\mathbf F=0}.
$$

When polarizability is expressed in SI units, this derivative has units of
C·m²/V. The Å³ convention used by MB-SPEC expresses polarizability as a
volume, obtained by dividing the SI polarizability by $4\pi\epsilon_0$ and
converting m³ to Å³. The coupled-dipole equations here are understood in the
model's consistent electrostatic units.

Differentiating the induced-dipole equations gives

$$
\frac{\partial\boldsymbol\mu_i^{\mathrm{ind}}}{\partial F_b}
=a_i\left[
\mathbf e_b+
\sum_{j\ne i}\mathbf T_{ij}^{\mathrm{damped}}
\frac{\partial\boldsymbol\mu_j^{\mathrm{ind}}}{\partial F_b}
\right],
$$

where $\mathbf e_b$ is the unit vector along Cartesian direction $b$.

Polarizability therefore follows from the same coupled dipole response.
MB-SPEC evaluates it analytically from the polarization matrix; it does not
perform separate finite-field calculations. See `ttm4_es::ptensor` in
[ttm4-es.cpp](tools/ttm4-es.cpp).

The site parameters are

$$
a_{\mathrm O}=1.310\ \text{Å}^3,\qquad
a_{\mathrm H}=0.294\ \text{Å}^3,\qquad
a_{\mathrm M}=0.
$$

Although these site parameters are isotropic and constant, their interactions
produce a geometry-dependent molecular or system polarizability tensor,
reported in Å³.

## Use in MB-SPEC

The total TTM4 polarizability already contains isolated-molecule
contributions. Raman's nonadditive contribution therefore subtracts them:

$$
\mathbf A_{\mathrm{nonadd}}^{\mathrm{TTM4}}
=\mathbf A_{\mathrm{system}}^{\mathrm{TTM4}}
-\sum_m\mathbf A_m^{\mathrm{TTM4}}.
$$

The full implemented Raman combination is

$$
\mathbf A^{\mathrm{MB\text{-}\alpha}}
=\sum_m\mathbf A_m^{\mathrm{Avila}}
+\mathbf A_{\mathrm{2B,corr}}
+\mathbf A_{\mathrm{nonadd}}^{\mathrm{TTM4}}.
$$

This subtraction prevents counting both the TTM4 and Avila monomer
polarizabilities. The nonadditive term includes two-body induction as well as
higher-body effects. See [polarizability.cpp](raman/polarizability.cpp) and
[alpha-ttm4.cpp](tools/alpha-ttm4.cpp).

The explicit short-range polarizability correction is added on top of
induction, whereas the MB-mu dipole construction replaces the corresponding
short-range induced-dimer contribution. See
[Medders and Paesani, Eqs. (7), (10), and (11)](https://doi.org/10.1021/ct501131j).

The internal cell-dependent implementation lacks Ewald electrostatics and is
not validated for periodic bulk calculations.

## References

- G. R. Medders and F. Paesani, *Infrared and Raman Spectroscopy of Liquid
  Water through “First-Principles” Many-Body Molecular Dynamics*, J. Chem.
  Theory Comput. **11**, 1145–1154 (2015).
  [DOI: 10.1021/ct501131j](https://doi.org/10.1021/ct501131j).
- C. J. Burnham, D. J. Anick, P. K. Mankoo, and G. F. Reiter, *The vibrational
  proton potential in bulk liquid water and ice*, J. Chem. Phys. **128**,
  154519 (2008). Introduces TTM4-F.
  [DOI: 10.1063/1.2895750](https://doi.org/10.1063/1.2895750).
