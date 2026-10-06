# Thermodynamics of Scale

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXXX)

Code, data, figures and manuscript for

> **R. Chen**, *Thermodynamics of Scale: Entropic c-Functions, Data Processing, Spectral Thermodynamics,
> Fluctuation Relations, and the Arrow of Time* (2026).
> DOI: [10.5281/zenodo.XXXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXXX)

## Summary

The paper sorts thermodynamic statements about the renormalization group (RG) into theorems, exact identities,
standard physics applied to scale, and conjecture, and checks every exact statement numerically.

- **RG irreversibility is information-theoretic.** The c, F and a theorems state c_UV > c_IR, i.e. a *decrease*
  towards the infrared. The entropic proofs use strong subadditivity and Lorentz invariance. For a massive lattice
  Dirac fermion the Casini–Huerta function c(ℓ) = 3ℓ dS/dℓ decreases monotonically to zero; at the critical point it
  reproduces c = 1 times the finite-size chord factor (πℓ/N)cot(πℓ/N).
- **H-theorem for RG channels.** Monotonicity of relative entropy under channels: if an RG step is a channel with an
  invariant state σ*, then D(Φⁿ(ρ)‖σ*) is non-increasing. This is checked on block restrictions of Gaussian
  fermion states.
- **Spectral thermodynamics of scale.** With heat-kernel time τ as inverse temperature,
  S(τ) = ln P + τ⟨λ⟩ = −D(π_τ‖ρ̄) is non-increasing, and the heat capacity C = τ²Var(λ) obeys
  C = d_s/2 − (1/2) dd_s/dln τ. Dimension is therefore twice a heat capacity on scaling plateaus. On ℤ^d,
  C = d x²(1 − r/x − r²) with x = 2τ and r = I₁/I₀; C has maximum 0.6799 d and tends to d/2.
- **Fluctuation relations for a physical scale protocol.** For a lattice field whose correlation length is changed
  in time (ξ: 10 → 1), the Jarzynski and Crooks relations hold. The sudden-quench Jarzynski estimator has finite
  variance if and only if ω_f²/ω_i² > 1/2 for every mode, so softening protocols fail. The Jarzynski and Bennett
  estimates agree with the exact ΔF = 7.1254 to within 0.015.
- **Landauer and the arrow of time.** Landauer's bound applies to physical erasure, not to the theorist's
  coarse-graining. The identification of the RG arrow with the arrow of time (dS/CFT, holographic c-theorems) is
  stated as a conjecture, which does not explain the low-entropy initial state.

This is the fifth paper of a series:
[1](https://doi.org/10.5281/zenodo.23189575) holographic entanglement of fractal regions
([code](https://github.com/Ruqing1963/holographic-fractal-entanglement)),
[2](https://doi.org/10.5281/zenodo.23193512) multifractal large deviations in quantum entanglement
([code](https://github.com/Ruqing1963/multifractal-quantum-entanglement)),
[3](https://doi.org/10.5281/zenodo.23194419) uncomputable fractal dimensions and the Chaitin barrier
([code](https://github.com/Ruqing1963/chaitin-fractal-complexity)),
[4](https://doi.org/10.5281/zenodo.23195881) crossover from discrete to continuous quantum geometry
([code](https://github.com/Ruqing1963/quantum-geometry-crossover)).

## Repository structure

| Path | Contents |
|---|---|
| `paper/` | LaTeX source and compiled PDF of the manuscript |
| `code/scale_thermodynamics.py` | Single script producing every number and figure in the paper |
| `figures/` | Figure as vector PDF (used by the paper) and PNG |
| `data/` | Numerical data as CSV (metadata in `#` header lines) |
| `results/` | Console output of the script (the numbers quoted in the paper) |

## Reproducing the results

Requirements: Python ≥ 3.10 and the packages in `requirements.txt`
(tested with Python 3.12.4, NumPy 1.26.4, SciPy 1.13.1, Matplotlib 3.8.4). Runtime is about 100 seconds; the
Langevin simulation uses 10⁵ trajectories per direction and a fixed random seed.

```bash
pip install -r requirements.txt
python code/scale_thermodynamics.py      # add --show to display the figure
```

Run from the repository root; the figure goes to `figures/` and data to `data/` (override with `ST_FIG_DIR`,
`ST_DATA_DIR`).

| Data file | Content | Paper |
|---|---|---|
| `entropic_c_function.csv` | c(ℓ) and S(ℓ) of the lattice Dirac fermion, N = 1200, m = 0, 0.05, 0.1, 0.2 | §2.1, Fig. 1 |
| `relative_entropy_vs_block.csv` | D(ρ_ℓ(m) ‖ σ_ℓ(0)) against block length | §2.2, Fig. 1 |
| `spectral_thermodynamics.csv` | ln P, ⟨λ⟩, Var(λ), S(τ), C(τ) of the 2D square lattice | §3, Fig. 1 |
| `jarzynski_crooks.csv` | mean and dissipated work, Jarzynski estimate ± s.e., Crooks slope, Bennett ΔF | §5, Table 1 |

To rebuild the paper (pdfLaTeX, two passes):

```bash
cd paper
pdflatex Chen_2026_Scale_Thermodynamics.tex
pdflatex Chen_2026_Scale_Thermodynamics.tex
```

## Citation

```bibtex
@misc{Chen2026ScaleThermodynamics,
  author = {Chen, Ruqing},
  title  = {Thermodynamics of Scale: Entropic c-Functions, Data Processing, Spectral Thermodynamics,
            Fluctuation Relations, and the Arrow of Time},
  year   = {2026},
  doi    = {10.5281/zenodo.XXXXXXXX},
  url    = {https://doi.org/10.5281/zenodo.XXXXXXXX}
}
```

## License

- **Code** (`code/`): [MIT License](LICENSE)
- **Manuscript, figures, data and results** (`paper/`, `figures/`, `data/`, `results/`):
  [CC BY 4.0](LICENSE-CC-BY-4.0.md)

## Contact

Ruqing Chen — GUT Geoservice Inc., Montreal — ruqing@hotmail.com
