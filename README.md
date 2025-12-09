# Althouse_Vector_Ecological_Specialization

## Ecological Specialization in Vectors  
### Reproducing the analyses and figures for  
**"Ecological specialization in vectors drives backward bifurcations and creates enzootic traps in multi-host, multi-vector systems."**

Benjamin M. Althouse  
*Submitted to PLoS Computational Biology* (2025)

---

## Repository structure

```
Althouse_Vector_Ecological_Specialization/
│
├─ Bifurcation/
│   ├─ Bifurcation_publication.R       # Generates Figure 2 (bifurcation diagram)
│   ├─ bifurcation_data_cache.rds      # Cached numerical bifurcation data
│   └─ output_2025-12-09/
│       └─ bifurcation_R0_two_panel_overlay.pdf  # Figure 2 output
│
├─ R0Compare/
│   ├─ CompareR0.R                              # Main text Figure 3
│   ├─ CompareR0_supplementary.R                # Supplement Fig S1 (normalized, unweighted biting)
│   ├─ CompareR0_supplementary_weighted_biting.R # Supplement Fig S2 (normalized, weighted biting)
│   ├─ r0numerical2specIntro.r                  # R₀ calculation: weighted FOI (Eq. 3)
│   └─ r0numerical2specIntroOldDenom.r          # R₀ calculation: unweighted FOI (Eq. 4)
│
└─ README.md
```

## Overview

This repository contains R code to reproduce the analytical results and figures from the manuscript on ecological specialization in vector-borne disease systems. The analysis examines how different force-of-infection (FOI) formulations affect the basic reproduction number (R₀) and system dynamics in multi-host, multi-vector disease models.

### Key analyses

1. **Bifurcation analysis** (`Bifurcation/`): Demonstrates backward bifurcations and multiple endemic equilibria when vectors exhibit ecological specialization (differential host preferences)

2. **R₀ surface comparison** (`R0Compare/`): Compares R₀ landscapes under weighted vs. unweighted FOI formulations across parameter space

## Quick-start

### Prerequisites
* **R ≥ 4.0** with packages: `RColorBrewer`, `fields`, `rgl`, `rstudioapi`

### Clone and run

```bash
git clone https://github.com/althouse/Althouse_Vector_Ecological_Specialization.git
cd Althouse_Vector_Ecological_Specialization
```

### Generate Figure 2 (Bifurcation diagram)

```bash
cd Bifurcation
Rscript Bifurcation_publication.R
```

Output: `output_YYYY-MM-DD/bifurcation_R0_two_panel_overlay.pdf`

### Generate Figure 3 (R₀ comparison)

```bash
cd R0Compare
Rscript CompareR0.R
```

Output: `output_YYYY-MM-DD/CompareR0-fig3-6x8in.pdf`

### Generate Supplement Figures

```bash
# Supplement Figure S1 (normalized, unweighted biting)
Rscript CompareR0_supplementary.R

# Supplement Figure S2 (normalized, weighted biting)  
Rscript CompareR0_supplementary_weighted_biting.R
```

Outputs: `output_YYYY-MM-DD_supp/` and `output_YYYY-MM-DD_supp_weighted/`

---

## Code details

### Bifurcation analysis

**`Bifurcation_publication.R`** generates a two-panel bifurcation diagram (Figure 2) showing endemic equilibria as a function of R₀, contrasting weighted (Eq. 3) and unweighted (Eq. 4) force-of-infection formulations.

The script:
1. Creates a date-stamped output folder (`output_YYYY-MM-DD/`)
2. Loads cached bifurcation data from previous ODE simulations
3. Calculates analytical R₀ for each transmission rate (β) value
4. Generates a 6.8 × 5-inch two-panel PDF showing prevalence vs. R₀
5. Demonstrates backward bifurcations and hysteresis effects

**Parameters** (based on Kédougou, Senegal sylvatic dengue system):
- Two host species: large primates (*Erythrocebus patas*), small primates (*Galago senegalensis*)
- Two vector species: *Aedes furcifer*, *Aedes taylori*
- Vectors exhibit strong ecological specialization (differential host preferences)
- Host populations: N_h = N_p = 1,000
- Vector populations: N_m1 = N_m2 = 25,000

### R₀ surface comparison

**`CompareR0.R`** produces Figure 3, a heat map comparing R₀ surfaces under the two FOI formulations across a grid of host population sizes.

Key features:
- Explores parameter space: large host (1,000–50,000) × small host (1,000–50,000)
- Calculates R₀ using both weighted and unweighted FOI denominators
- Highlights regions where formulations produce different epidemic predictions
- Includes field data ellipse showing empirical uncertainty

**`CompareR0_supplementary.R`** and **`CompareR0_supplementary_weighted_biting.R`** generate supplementary figures showing that when total biting is normalized (equal-bite constraint), the two formulations produce nearly identical R₀ surfaces, confirming their algebraic equivalence under conserved total biting.

### R₀ calculation functions

**`r0numerical2specIntro.r`** implements the weighted FOI formulation (Eq. 3):
- FOI denominator weighted by host abundance
- Function: `r0new()`
- Closed-form expression for spectral radius of next-generation matrix
- Based on Althouse et al. (2012) PLOS NTD analytical derivation

**`r0numerical2specIntroOldDenom.r`** implements the unweighted FOI formulation (Eq. 4):
- FOI denominator uses total bites (unweighted by host abundance)
- Function: `r0old()`
- Alternative formulation showing different biological constraints

Both functions take the same epidemiological parameters:
- Transmission probabilities (β)
- Vector feeding preferences (r)
- Host and vector abundances (N)
- Mortality and recovery rates (μ, γ)

---

## Model background

The code implements a multi-host, multi-vector SIR (Susceptible-Infected-Recovered) model for vector-borne disease transmission. The key biological insight is that **vector ecological specialization** (differential host-feeding preferences) can:

1. Create **backward bifurcations** where R₀ > 1 is necessary but not sufficient for disease persistence
2. Generate **multiple endemic equilibria** and hysteresis effects
3. Trap systems in **enzootic cycles** even when control efforts reduce R₀ below traditional epidemic thresholds

The choice of FOI denominator (weighted vs. unweighted by host density) affects quantitative predictions about disease dynamics, particularly when vectors exhibit strong host preferences.

---

## How to cite

If you use any part of this code or data, please cite the associated article:

> Althouse BM (2025) *Ecological specialization in vectors drives backward bifurcations and creates enzootic traps in multi-host, multi-vector systems.* PLoS Computational Biology (submitted).

And the original model derivation:

> Althouse BM, Lessler J, Sall AA, Diallo M, Hanley KA, et al. (2012) *Synchrony of Sylvatic Dengue Isolations: A Multi-Host, Multi-Vector SIR Model of Dengue Virus Transmission in Senegal.* PLOS Neglected Tropical Diseases 6(11): e1928. [doi:10.1371/journal.pntd.0001928](https://doi.org/10.1371/journal.pntd.0001928)

---

## Contact

Questions or pull requests are welcome!  
**Benjamin M. Althouse** — bma85@uw.edu

---

## License

This repository is released under the **MIT License**.

```
MIT License

Copyright (c) 2025 Benjamin M. Althouse

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## Acknowledgements

Field data collection by the Dakar Pasteur Institute vector team. Bifurcation analysis inspired by discussions with Derek Cummings and David Smith. Funding from NIH AI069145, NSF GRFP DGE-0707427, Santa Fe Institute, and the Global Good Fund.
