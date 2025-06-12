# Veteran Lung Cancer Survival Analysis

🩺 This project explores survival outcomes of veterans with advanced lung cancer using the **Cox Proportional Hazards model**.

This dataset is a classic survival analysis example from a randomized trial comparing two treatments. The project is built with:

- 📦 R + `{survival}`, `{survminer}`, `{ggplot2}`, `{broom}`, `{splines}`
- 📁 Project-local dependency management via `{renv}`
- 📓 Reproducible analysis with Quarto notebooks (`.qmd`)

---

## ✅ Project Overview

- **Teaching Notebook** (`notebooks/veteran-lung-cancer-coxph.qmd`):  
  An annotated Quarto notebook with callouts to explain each step and help interpret graphs and model outputs.
- **Final Report Script** (`R/cox-ph.R`):  
  The finalized R script containing the full analysis pipeline and model selection, suitable for reporting or publication.
- **Setup Script** (`scripts/setup.R`):  
  Ensures all dependencies are installed and the project environment is reproducible across platforms.

---

## 📂 Structure

```
├── notebooks/
│   └── veteran-lung-cancer-coxph.qmd
├── R/
│   └── cox-ph.R
├── scripts/
│   └── setup.R
├── renv/
├── renv.lock
├── .gitignore
├── README.md
└── veteran-lung-cancer-coxph.Rproj
```

---

## 📜 License

MIT