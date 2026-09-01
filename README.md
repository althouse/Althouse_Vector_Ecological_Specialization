# Ecological specialization in vectors

Reproducible code for the accepted PLOS Computational Biology article:

> Benjamin M. Althouse. **Ecological specialization in vectors alters transmission thresholds and endemic dynamics in multi-host, multi-vector systems.** PLOS Computational Biology (accepted, 2026; DOI pending).

## Overview

This repository compares two force-of-infection (FOI) formulations for a two-host, two-vector SIR-SI model:

- **Weighted FOI:** vectors disproportionately contact preferred hosts.
- **Unweighted FOI:** vectors encounter the full host community opportunistically.

With the same nominal parameters, the formulations can produce different epidemic thresholds, long-run prevalence, and oscillation amplitudes. Both autonomous formulations undergo a **forward (supercritical) transcritical bifurcation at `R0 = 1`**. The model does not support a backward bifurcation or hysteresis under the assumptions analyzed here.

An equal-bite rescaling makes the complete ODE systems, and therefore their `R0` surfaces, identical. This confirms that the unnormalized comparison reflects different total contact rates implied by the two ecological hypotheses.

## Repository contents

```text
.
├── Bifurcation/
│   ├── Bifurcation_publication.R
│   ├── bifurcation_data_cache.rds
│   └── compute_bifurcation_coefficients.R
├── R0Compare/
│   ├── CompareR0.R
│   ├── CompareR0_supplementary.R
│   ├── CompareR0_supplementary_weighted_biting.R
│   ├── figure_helpers.R
│   ├── r0numerical2specIntro.r
│   ├── r0numerical2specIntroOldDenom.r
│   └── test_normalization.R
├── outputs/
│   ├── Figure2_prevalence_vs_R0.pdf
│   ├── Figure3_R0_comparison.pdf
│   ├── FigureS1_equal_bite_R0_comparison.pdf
│   └── FigureS2_asymmetric_transmission.pdf
├── CITATION.cff
├── LICENSE
└── run_all.R
```

### Analysis-to-figure map

| Script | Manuscript result | Output |
|---|---|---|
| `Bifurcation/Bifurcation_publication.R` | Main Figure 2: periodically forced prevalence samples versus `R0` | `outputs/Figure2_prevalence_vs_R0.pdf` |
| `Bifurcation/compute_bifurcation_coefficients.R` | Center-manifold coefficients; forward transcritical bifurcation for both FOIs | Console report and numerical assertions |
| `R0Compare/CompareR0.R` | Main Figure 3: unnormalized weighted versus unweighted `R0` landscapes | `outputs/Figure3_R0_comparison.pdf` |
| `R0Compare/CompareR0_supplementary_weighted_biting.R` | Supplementary Figure S1: equal-bite identity | `outputs/FigureS1_equal_bite_R0_comparison.pdf` |
| `R0Compare/CompareR0_supplementary.R` | Supplementary Figure S2: asymmetric transmission (`beta_hv = beta_vh / 4`) | `outputs/FigureS2_asymmetric_transmission.pdf` |
| `R0Compare/test_normalization.R` | Regression test of the equal-bite identity | Pass/fail console report |

The vertical point clouds in Figure 2 are annual samples from a single periodically forced trajectory at each parameter value. They represent temporal variation, not coexisting equilibria.

## Requirements

- R 4.1 or newer
- R package `RColorBrewer`

Install the only non-base dependency with:

```r
install.packages("RColorBrewer")
```

The scripts use relative paths, run without RStudio, and do not require a graphical display.

## Reproduce everything

From the repository root:

```bash
Rscript run_all.R
```

This runs both numerical checks and regenerates all four PDFs in `outputs/`.

To run components individually:

```bash
Rscript Bifurcation/compute_bifurcation_coefficients.R
Rscript R0Compare/test_normalization.R
Rscript Bifurcation/Bifurcation_publication.R
Rscript R0Compare/CompareR0.R
Rscript R0Compare/CompareR0_supplementary_weighted_biting.R
Rscript R0Compare/CompareR0_supplementary.R
```

## Principal parameters

- Host populations for the baseline comparison: 1,000 per species
- Vector populations: 25,000 per species
- Primary-host biting rate: `0.5 day^-1`
- Cross-host biting rate: `0.001 day^-1`
- Host recovery rate: `1/4 day^-1`
- Large-host mortality: `1/(60 * 365) day^-1`
- Small-host mortality: `1/(15 * 365) day^-1`
- Vector mortality: `1/7 day^-1`

The Figure 3 grids cover large- and small-host abundances from 1 to 5,000, vector abundance from 1 to 50,000, and transmission probabilities from 0 to 0.5. Red ellipses show the field-based Kédougou parameter estimates and uncertainty ranges used in the article.

## Reproducibility notes

- `Bifurcation/bifurcation_data_cache.rds` contains the annual prevalence samples needed for Figure 2, so the original long ODE simulations are not required to reproduce that plot.
- `r0numerical2specIntro.r` and `r0numerical2specIntroOldDenom.r` implement the weighted and unweighted closed-form `R0` expressions, respectively.
- The equal-bite transformation uses `c_j = Nbar / D_j`, where `D_j` is calculated from normalized preference weights. The regression test fails if the two formulations differ by more than `1e-10` in the tested scenarios.
- Generated outputs use stable filenames so reruns are directly comparable.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Until the article DOI is assigned, please cite:

> Althouse BM (2026). *Ecological specialization in vectors alters transmission thresholds and endemic dynamics in multi-host, multi-vector systems.* PLOS Computational Biology, accepted.

The underlying model was introduced in:

> Althouse BM, Lessler J, Sall AA, Diallo M, Hanley KA, et al. (2012). Synchrony of Sylvatic Dengue Isolations: A Multi-Host, Multi-Vector SIR Model of Dengue Virus Transmission in Senegal. *PLOS Neglected Tropical Diseases* 6(11): e1928. https://doi.org/10.1371/journal.pntd.0001928

## License

The code is released under the [MIT License](LICENSE).

## Contact

Benjamin M. Althouse — bma85@uw.edu
