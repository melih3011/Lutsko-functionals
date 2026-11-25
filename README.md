# Lutsko functionals for hard spheres: LK, LK-WB and LK-WBII

C implementations of Lutsko's two-parameter version of **Fundamental Measure Theory (FMT)** for the
hard-sphere (HS) fluid. The repository contains the original Rosenfeld-based form **LK(A, B)** and the
White-Bear-based extensions **LK-WB(A, B)** and **LK-WBII(A, B)**. There is a density-functional solver for each of
several geometries: the test particle, the planar hard wall, the spherical slit and the spherical cavity. Each solver
can be linked against any of the functionals.

The free parameters $`A`$ and $`B`$ were fixed using two **test-particle sum rules**, for the excess chemical potential
and the isothermal compressibility. The code here produced the results of the following papers:

1. **M. Gül, R. Roth, R. Evans**, *Using test particle sum rules to construct accurate functionals in classical
   density functional theory*, [Phys. Rev. E **110**, 064115 (2024)](https://doi.org/10.1103/PhysRevE.110.064115).
2. **M. Gül, R. Roth, R. Evans**, *Using test particle sum rules to improve approximations in classical density
   functional theory: White-Bear and White-Bear mark II versions of the Lutsko functional*,
   [Phys. Rev. E **113**, 034104 (2026)](https://doi.org/10.1103/3vg8-mbf4),
   [arXiv:2604.07050](https://arxiv.org/abs/2604.07050).

The functional itself was introduced in

- **J. F. Lutsko**, *Explicitly stable fundamental-measure-theory models for classical density functional theory*,
  [Phys. Rev. E **102**, 062137 (2020)](https://doi.org/10.1103/PhysRevE.102.062137).

---

## Theory in brief

### Classical DFT and FMT

For a fluid in an external potential $`V_\text{ext}(\mathbf r)`$, the equilibrium density $`\rho(\mathbf r)`$
minimises the grand-potential functional

```math
\Omega[\rho] = \mathcal F_\text{id}[\rho] + \mathcal F_\text{ex}[\rho]
             + \int \mathrm d\mathbf r  \rho(\mathbf r)\big(V_\text{ext}(\mathbf r) - \mu\big).
```

In FMT the HS excess free energy is a local function of weighted densities,
$`\beta\mathcal F_\text{ex} = \int \mathrm d\mathbf r  \Phi(\{n_\alpha(\mathbf r)\})`$. Each weighted density is a
convolution $`n_\alpha = \rho * w_\alpha`$ with a weight function that encodes the geometry of a single sphere of radius $`R=\sigma/2`$:
the scalar weights $`n_0,\dots,n_3`$, the vector weights $`\mathbf n_1, \mathbf n_2`$, and the tensor
$`\mathbf T = \rho * \big(\hat{\mathbf r}\otimes\hat{\mathbf r} \delta(R-r)\big)`$.
The free-energy density is split as $`\Phi = \Phi_1 + \Phi_2 + \Phi_3`$, with the Rosenfeld (RF) terms

```math
\Phi_1 = -n_0\ln(1-n_3), \qquad
\Phi_2 = \frac{n_1 n_2 - \mathbf n_1\cdot\mathbf n_2}{1-n_3}.
```

### The Lutsko functional LK(A, B)

Lutsko (2020) generalised the third term by adding two free parameters:

```math
\Phi_3^\text{LK} = f(n_3)\Big[(A+B) n_2^3 - 3A n_2 \mathbf n_2\cdot\mathbf n_2
 + 3A \mathbf n_2\mathbf T\mathbf n_2 - 3B n_2 \mathrm{Tr} \mathbf T^2 + (2B-A) \mathrm{Tr} \mathbf T^3\Big],
\qquad f^\text{RF}(n_3) = \frac{1}{24\pi(1-n_3)^2}.
```

The choice $`A=-B=3/2`$ recovers Tarazona's tensor version of Rosenfeld's functional. In bulk the equation of state is

```math
\frac{\beta P^\text{LK}}{\rho} = \frac{\beta P^\text{PY}}{\rho} + (8A+2B-9) \frac{\eta^2}{3(1-\eta)^3},
\qquad \eta = \frac{\pi}{6}\rho\sigma^3,
```

so every point on the **PY line** $`8A+2B=9`$ reproduces the Percus–Yevick thermodynamics.

### White-Bear extensions LK-WB(A, B) and LK-WBII(A, B)  (paper 2)

The White-Bear (WB) and White-Bear mark II (WBII) functionals are built around the more accurate Carnahan–Starling (CS)
equation of state. Following the "FMT toolbox" idea, the Lutsko numerator is combined with their $`n_3`$-dependence:

```math
\Phi^\text{LK-WB} = \Phi_1 + \Phi_2 + \Phi_3^\text{LK}\big|_{f\to f^\text{WB}},
\qquad f^\text{WB}(n_3) = \frac{n_3 + (1-n_3)^2\ln(1-n_3)}{36\pi n_3^2(1-n_3)^2},
```

```math
\Phi^\text{LK-WBII} = \Phi_1 + \varphi_1(n_3) \Phi_2 + \varphi_2(n_3) \Phi_3^\text{LK},
```

```math
\varphi_1 = 1 + \frac{2n_3 - n_3^2 + 2(1-n_3)\ln(1-n_3)}{3n_3},\qquad
\varphi_2 = 1 - \frac{2n_3 - 3n_3^2 + 2n_3^3 + 2(1-n_3)^2\ln(1-n_3)}{3n_3^2}.
```

Both reduce to the tensor WB and WBII functionals for $`A=-B=3/2`$. Their equations of state are

```math
\frac{\beta P^\text{LK-WB}}{\rho} = \frac{\beta P^\text{CS}}{\rho} + (8A+2B-9)\frac{\eta^2(3-\eta)}{9(1-\eta)^3},
\qquad
\frac{\beta P^\text{LK-WBII}}{\rho} = \frac{\beta P^\text{CS}}{\rho} + (8A+2B-9)\frac{\eta^2(3-2\eta+\eta^2)}{9(1-\eta)^3},
```

so the line $`8A+2B=9`$ now becomes the **CS line**.

### Fixing A and B with test-particle sum rules

Fix one particle at the origin. It then acts as an external potential equal to the pair potential, and the
equilibrium profile satisfies $`\rho(r;\phi)=\rho_b g(r)`$ (Percus' test-particle route). This gives two independent DFT routes to bulk
thermodynamics, which can be compared with the bulk free energy $`\Phi(\rho_b)`$ of the same functional:

* **Excess chemical potential:** $`\beta\mu_\text{ex}^\text{DFT} = \beta \Delta\Omega`$, the grand-potential cost of
  inserting the test particle.
* **Reduced isothermal compressibility:**
  $`\chi_T^\text{DFT} = S(0) = 1 + 4\pi\rho_b\int_0^\infty \mathrm dr  r^2 \big(g(r)-1\big)`$.

An exact functional satisfies both sum rules. For an approximate one, the relative deviations

```math
\delta_\mu = \frac{\mu_\text{ex}^\text{DFT}-\mu_\text{ex}^\text{Bulk}}{\mu_\text{ex}^\text{Bulk}},\qquad
\delta_\chi = \frac{\chi_T^\text{DFT}-\chi_T^\text{Bulk}}{\chi_T^\text{Bulk}}
```

measure its internal inconsistency. Minimising $`M=\tfrac12(|\delta_\mu|+|\delta_\chi|)`$ in the $`(A,B)`$ plane,
averaged over the high packing fractions $`\eta = 0.35, 0.40, 0.45`$, gives the optimal parameters:

| Functional | $`A`$  | $`B`$   | reference        |
|------------|------|-------|------------------|
| LK         | 1.30 | −1.00 | paper 1 (PY-based) |
| LK-WB      | 1.35 | −0.85 | paper 2 (CS-based) |
| LK-WBII    | 1.25 | −0.20 | paper 2 (CS-based) |

LK-WBII(1.25, −0.20) improves markedly on WBII in consistency with the sum rules.

The functionals are also checked against simulation and exact results:
* hard-wall and test-particle density profiles against Monte Carlo data;
* the contact theorem and Gibbs adsorption theorem;
* the pair direct correlation function $`c^{(2)}(r)`$;
* virial coefficients;
* the equation of state.

A stringent stability test is HS confined in a small hard spherical cavity. There the optimised Lutsko functionals
converge, while plain RF and WB fail. Their $`q_3`$ (antisymmetrised) variants are included for comparison.

---

## Repository layout

```
.
├── src/                         C solvers: one folder per geometry
│   ├── rules.mk                 shared build rules
│   ├── test_particle/           hard sphere fixed at origin: sum rules for mu_ex and chi_T   (papers 1 & 2)
│   ├── hard_wall/               planar hard wall: profiles, contact & adsorption theorems
│   ├── slit/                    spherical slit: fluid between a hard sphere and a concentric hard shell
│   ├── spherical_cavity/        HS confined in a hard spherical cavity: stability test           (paper 2)
│   └── direct_correlation/      analytic bulk c2(r) for PY and LK(A,B)
│       each geometry folder:    <solver>.c, functionals/*.c, Makefile, helper scripts
├── mathematica/                 symbolic derivations for RF, WB, WBII, LK, LK-WB, LK-WBII
├── analysis/                    post-processing (e.g. optimal A,B vs eta)
├── docs/notes/                  working notes: formulas, early project write-up, numerics
├── legacy/                      older versions of the solvers kept for reference
└── data/                        raw solver output and MC reference data (local only, see data/README.md)
```

Raw data (`data/`) are git-ignored; only `data/README.md` is versioned, so the expected layout is still
documented. Manuscripts and figure data of the papers are not part of this repository - see the published articles above.

### Functional interface

All functionals in a `functionals/` folder share one C interface: the solver calls them and never needs to
know which functional is linked.

| function | meaning |
|---|---|
| `pressure`, `chemical_potential`, `sk0_bulk`, `detadmu` | bulk thermodynamics for given $`\eta, A, B`$ |
| `myphi1`, `myphi2`, `myphi3` | $`\Phi_1, \Phi_2, \Phi_3`$ |
| `myd2`, `myd3`, `mydv2`, `mydt1`, `mydt2` | $`\partial\Phi/\partial n_\alpha`$ (scalar, vector, tensor) |
| `filename` | tag used in output file names (`rf`, `wb`, `wb_II`, `lk`, `lkwb`, `lkwb_II`, `rfq3`, `wbq3`) |

The plain RF, WB and WBII files take $`A, B`$ in their signatures but ignore them: pass `0 0`. Each geometry keeps its
own copy of the functionals, because the copies differ in small solver-specific details (for example, which helper
functions are compiled in).

---

## Building

Requirements: a C compiler and `make`. The scripts also need Python 3 with `numpy` and `matplotlib` (plus
`scienceplots` for the plotting styles). The notebooks need Mathematica.

```bash
make -C src            # every geometry with every functional
make -C src/hard_wall  # one geometry
make -C src/test_particle lkwb2   # one executable, e.g. test particle + LK-WBII
make -C src clean
```

Executables keep their historical names so that the existing scripts still work:

| suffix | functional | suffix | functional |
|---|---|---|---|
| `_rf`  | Rosenfeld            | `_lk`    | LK(A,B) |
| `_wb`  | White-Bear           | `_lkwb`  | LK-WB(A,B) |
| `_wb2` | White-Bear mark II   | `_lkwb2` | LK-WBII(A,B) |
| `_rfq3`, `_wbq3` | $`q_3`$ versions of RF / WB (cavity only) | | |

The prefix is `spone_` for the spherical geometries and `hwone_` for the hard wall.

## Running

**Units.** Inside the codes all lengths are in units of the HS radius $`R=\sigma/2`$; the grid spacing is $`R/100`$,
or $`R/250`$ for the cavity. The first column of every output file is the distance in units of $`\sigma`$.
The last argument `a` is the Picard mixing parameter (default 0.1); reduce it at high packing fractions.

### Test particle (sum rules), `src/test_particle`

```bash
./spone_lkwb2  Rb  eta  A  B  mode  [a]
```
`Rb` is the radius of the fixed particle in units of $`R`$, so `1.0` means an identical test particle.

| `mode` | action | output |
|---|---|---|
| 0 | single state point | `sp_<f>_<eta>_<Rb>`: header with $`\mu_\text{ex}`$ and $`S(0)`$ from both routes, then $`r/\sigma`$, $`\rho/\rho_b`$, $`n_3`$, $`\Phi`$, $`c^{(1)}`$. Prints contact and adsorption checks. |
| 1 | scan $`\eta = 0.01 \dots 0.48`$ | `<f>_run_eta`: $`\eta`$, $`\mu_\text{ex}^\text{DFT}`$, $`\mu_\text{ex}^\text{Bulk}`$, $`\delta_\mu`$, $`\chi_T^\text{DFT}`$, $`\chi_T^\text{Bulk}`$, $`\delta_\chi`$ |
| 2 | scan the $`(A,B)`$ grid | stdout: `A B delta_mu delta_chi` (heat maps of paper 1) |

Example: `./spone_lkwb2 1.0 0.40 1.25 -0.2 0 0.01`.
`run_py_plane.py` runs mode 2 over $`\eta`$ and calls `find_min.py` on each scan. `avg_dev.py` averages the deviations over packing fractions.

### Planar hard wall, `src/hard_wall`

```bash
./hwone_lkwb2  eta  A  B  [a]
```
Writes `hw_<f>_<eta>_<A>_<B>` with columns $`z/\sigma`$, $`\rho\sigma^3`$, $`n_3`$, $`c^{(1)}`$, residual. Prints
the contact theorem ($`\rho(0^+)=\beta P`$) and the Gibbs adsorption theorem from both routes. MC reference profiles
are in `data/simulation/hard_wall_MC`.

### Spherical cavity, `src/spherical_cavity`

```bash
./spone_lkwb2  Rb  eta  A  B  print_mode  [a]
```
`Rb` is the cavity radius in units of $`R`$; paper 2 uses `Rb = 4`, i.e. $`2\sigma`$. `print_mode 1` writes the profile
`sp_<f>_<Rb>_<eta>_<A>_<B>` ($`r/\sigma`$, $`\rho\sigma^3`$, residual). It also prints adsorption, contact density, $`\langle N\rangle`$,
surface tension and excess volume. `print_mode 0` prints one line, `Rb gamma v_ex`, which is used by
`cavity_radius_scan.py`. Paper 2, Fig. 6:

```bash
./spone_lk     4 0.32 1.30 -1.00 1
./spone_lkwb   4 0.32 1.35 -0.85 1
./spone_lkwb2  4 0.32 1.25 -0.20 1
```

### Spherical slit, `src/slit`

```bash
./spone_lk  Nr  Rc  eta  A  B  [a]
```
Hard spheres between a hard sphere of radius `Rc` and a concentric hard spherical wall at `Nr`, both in units of
$`R`$. Writes `sp_<f>_<Nr>_<Rc>_<eta>_<A>_<B>`.

### Direct correlation function, `src/direct_correlation`

```bash
./corrf eta     # r/sigma, c2_PY, c2_LK(1.3,-1.0), c2_LK(1.0,0.0), c2_LK(0.0,4.0)
```

## Reproducibility notes and known issues

* The planar hard wall and the cavity solver reproduce stored profiles exactly. The test-particle solver reproduces
  them up to the density normalisation: older files store $`\rho`$, current ones $`\rho/\rho_b`$. Cavity profiles written
  before Feb 2026 store $`\rho/\rho_b`$; the current `sp_cavity.c` writes $`\rho\sigma^3`$.
* **Spherical slit:** the current build returns `nan` on macOS/clang. The stored profiles in `data/slit` came from a
  run that converged. An uninitialised-memory bug was fixed (arrays now use `calloc`), but the divergence remains
  and is being investigated.
* The hard-wall solver has no plain Rosenfeld functional.

## How to cite

If you use this code, please cite paper 1 and, for LK-WB/LK-WBII, paper 2, together with Lutsko's original work:

```bibtex
@article{Gul2024,
  title   = {Using test particle sum rules to construct accurate functionals in classical density functional theory},
  author  = {G\"ul, Melih and Roth, Roland and Evans, Robert},
  journal = {Phys. Rev. E}, volume = {110}, number = {6}, pages = {064115}, year = {2024},
  doi     = {10.1103/PhysRevE.110.064115}
}
@article{Gul2026,
  title   = {Using test particle sum rules to improve approximations in classical density functional theory:
             White-Bear and White-Bear mark II versions of the Lutsko functional},
  author  = {G\"ul, Melih and Roth, Roland and Evans, Robert},
  journal = {Phys. Rev. E}, volume = {113}, number = {3}, pages = {034104}, year = {2026},
  doi     = {10.1103/3vg8-mbf4}
}
@article{Lutsko2020,
  title   = {Explicitly stable fundamental-measure-theory models for classical density functional theory},
  author  = {Lutsko, James F.},
  journal = {Phys. Rev. E}, volume = {102}, number = {6}, pages = {062137}, year = {2020},
  doi     = {10.1103/PhysRevE.102.062137}
}
```

See also [`CITATION.cff`](CITATION.cff).

### Further background

* Y. Rosenfeld, *Phys. Rev. Lett.* **63**, 980 (1989): original FMT.
* P. Tarazona, *Phys. Rev. Lett.* **84**, 694 (2000): tensor weighted densities.
* R. Roth, R. Evans, A. Lang, G. Kahl, *J. Phys.: Condens. Matter* **14**, 12063 (2002); Y.-X. Yu, J. Wu,
  *J. Chem. Phys.* **117**, 10156 (2002): White-Bear.
* H. Hansen-Goos, R. Roth, *J. Phys.: Condens. Matter* **18**, 8413 (2006): White-Bear mark II.
* R. Roth, *J. Phys.: Condens. Matter* **22**, 063102 (2010): FMT review.
* J. K. Percus, *Phys. Rev. Lett.* **8**, 462 (1962): test-particle route.
* J. R. Henderson, in *Fundamentals of Inhomogeneous Fluids* (1992): statistical-mechanical sum rules.

## Authors

Melih Gül, Roland Roth (Institute for Theoretical Physics, University of Tübingen) and Robert Evans (H. H. Wills Physics
Laboratory, University of Bristol).

## License

[MIT](LICENSE)
