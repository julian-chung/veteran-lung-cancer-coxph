# Veteran Lung Cancer Survival Analysis

A survival analysis of the U.S. Veterans Administration lung cancer trial
(Kalbfleisch & Prentice, 1980), comparing standard vs experimental chemotherapy
in 137 male patients with advanced, inoperable lung cancer.

Built as a portfolio piece demonstrating applied survival analysis in R —
Kaplan–Meier estimation, Cox proportional hazards modelling, functional form
checking, spline transformations, interaction testing, diagnostics, and the
Aalen additive model.

---

## Deliverables

| File | Purpose |
|---|---|
| `notebooks/veteran-lung-cancer-overview.qmd` | 1-page visual summary: 4 plots telling the full story at a glance |
| `notebooks/veteran-lung-cancer-coxph-report.qmd` | Full written report: step-by-step analysis, model building, and interpretation |
| `notebooks/veteran-lung-cancer-coxph.qmd` | Teaching & reference notebook: same analysis with statistical theory and code explanations in collapsible callouts |
| `R/cox-ph.R` | All analysis code extracted from the notebooks — prose-free, executable, organized by section |

## Final Model

```r
coxph(Surv(time, status) ~ trt + prior + ns(age, df = 3) + karno + celltype,
      data = veteran)
```

Concordance: **0.735**. Key findings: no significant treatment effect
(HR ≈ 1.4, p = 0.097); Karnofsky score (HR = 0.966 per point, p < 0.001) and
cell type (adenocarcinoma HR ≈ 3.3, small cell HR ≈ 2.2 vs squamous, both
p < 0.01) are the dominant prognostic factors.

---

## Stack

- R with `survival`, `survminer`, `ggplot2`, `splines`, `forestmodel`, `ggfortify`
- Quarto for reproducible notebook rendering
- `renv` for dependency management

## Reproducing the Analysis

```r
# Restore the package environment
renv::restore()

# Then render the site
# quarto render   (from terminal)
# or use the Build pane in RStudio
```

---

## Project Structure

```
├── notebooks/
│   ├── veteran-lung-cancer-overview.qmd
│   ├── veteran-lung-cancer-coxph-report.qmd
│   └── veteran-lung-cancer-coxph.qmd
├── R/
│   └── cox-ph.R
├── references/
├── docs/                  # rendered site (GitHub Pages)
├── renv.lock
└── veteran-lung-cancer-coxph.Rproj
```

---

## License

MIT
