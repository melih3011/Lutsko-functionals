# data/ (local only, git-ignored)

Raw solver output and reference data. Only this README is versioned; regenerate the data with the solvers in `src/`.

| folder | content | produced by |
|---|---|---|
| `test_particle/` | test-particle profiles `sp_<f>_<eta>_<Rb>`, eta scans `<f>_run_eta`, (A,B) scans (`runs_*`, `lkwb*`, `mark2`, `deviations`), `averages` | `src/test_particle` |
| `hard_wall/` | hard-wall profiles `hw_<f>_<eta>_<A>_<B>`; `checks/` contact/adsorption tests | `src/hard_wall` |
| `slit/` | spherical-slit profiles `sp_<f>_<Nr>_<Rc>_<eta>_<A>_<B>` | `src/slit` |
| `spherical_cavity/` | `profiles/` cavity profiles, `scans/` gamma and v_ex vs cavity radius, `wb_sphere/`, `plots/`, `legacy_rr_code/` | `src/spherical_cavity` |
| `direct_correlation/` | bulk c2(r) from DFT and MC (`corrs_*`), analytic PY/LK curves | `src/direct_correlation` |
| `bulk_pressure/` | pressure of LK, LK-WB, LK-WBII vs eta (`pLK*.txt`) and plots | Mathematica notebooks |
| `simulation/hard_wall_MC/` | Monte Carlo HS density profiles at a planar hard wall (`hswall<eta*1000>.dat`) | MC simulation (external) |

Normalisation of the profile column: files written by older code versions store rho or rho/rho_b; current code writes
rho/rho_b (test particle) or rho*sigma^3 (hard wall, spherical cavity). See the main README.
