# =============================================================================
# Veteran Lung Cancer Survival Analysis — Reference Script
# =============================================================================
# All R code from the analysis notebooks in executable form.
# Use this for quick reference without scrolling through prose.
#
# Source notebooks:
#   notebooks/veteran-lung-cancer-coxph.qmd       (full analysis)
#   notebooks/veteran-lung-cancer-overview.qmd    (summary visualisations)
#
# Final model: trt + prior + ns(age, df=3) + karno + celltype
# =============================================================================


# ---- 1. Libraries ------------------------------------------------------------

library(survival)     # Surv(), coxph(), survfit(), aareg()
library(survminer)    # ggsurvplot(), ggforest(), ggcoxfunctional(), ggadjustedcurves()
library(dplyr)        # data wrangling
library(forcats)      # factor helpers
library(broom)        # tidy model output
library(purrr)        # map_dfr for univariable loop
library(ggplot2)      # plotting
library(splines)      # ns() for natural cubic splines
library(glue)         # glue() for inline stats output
library(forestmodel)  # forest_model() for Cox forest plots
library(ggfortify)    # autoplot() for Aalen model


# ---- 2. Data -----------------------------------------------------------------

# Used in full analysis (teaching notebook)
# trt: Control / Intervention
veteran <- survival::veteran %>%
  mutate(trt = factor(trt,
                      levels = c(1, 2),
                      labels = c("Control", "Intervention")))

# Used in overview (all variables fully labelled)
# trt: Standard / Experimental
veteran2 <- survival::veteran %>%
  mutate(
    trt      = factor(trt, levels = c(1, 2), labels = c("Standard", "Experimental")),
    celltype = factor(celltype, levels = c("squamous", "smallcell", "adeno", "large"),
                      labels = c("Squamous", "Small cell", "Adenocarcinoma", "Large cell")),
    prior    = factor(prior, levels = c(0, 10), labels = c("No", "Yes"))
  )


# ---- 3. Exploratory: Kaplan-Meier & Log-rank ---------------------------------

surv_obj <- Surv(time = veteran$time, event = veteran$status)

km_fit <- survfit(surv_obj ~ trt, data = veteran)
summary(km_fit)

ggsurvplot(
  km_fit,
  data             = veteran,
  conf.int         = TRUE,
  pval             = TRUE,
  pval.size        = 4,
  risk.table       = TRUE,
  risk.table.title = "Number at Risk",
  risk.table.col   = "strata",
  legend.title     = "Treatment Group",
  legend.labs      = c("Control", "Intervention"),
  xlab             = "Time (days)",
  ylab             = "Survival Probability",
  palette          = c("#A569BD", "#45B39D"),
  title            = "Kaplan-Meier Survival Curve: VA Lung Cancer Trial",
  ggtheme          = theme_minimal(base_size = 13)
)

# Log-rank test
log_rank_test <- survdiff(surv_obj ~ trt, data = veteran)
log_rank_test
1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1)


# ---- 4. Univariable Cox Models -----------------------------------------------

# Treatment only
cox_model <- coxph(surv_obj ~ trt, data = veteran)
summary(cox_model)

ggadjustedcurves(
  cox_model,
  data         = veteran,
  variable     = "trt",
  method       = "average",
  conf.int     = TRUE,
  legend.labs  = c("Control", "Intervention"),
  legend.title = "Treatment Group",
  xlab         = "Time (days)",
  ylab         = "Adjusted Survival Probability",
  palette      = c("#A569BD", "#45B39D"),
  ggtheme      = theme_minimal(base_size = 13)
)

# PH assumption
cox.zph(cox_model)
plot(cox.zph(cox_model))

# Martingale residuals
veteran$martingale  <- residuals(cox_model, type = "martingale")
veteran$linear_pred <- predict(cox_model, type = "lp")

ggplot(veteran, aes(x = linear_pred, y = martingale, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  labs(title = "Martingale Residuals vs Linear Predictor",
       x = "Linear Predictor (Log Hazard)", y = "Martingale Residual",
       color = "Treatment Group") +
  theme_minimal(base_size = 13)


# ---- 5. Model Building: Age & Prior ------------------------------------------

# trt + age
cox_age <- coxph(Surv(time, status) ~ trt + age, data = veteran)
summary(cox_age)

cox_zph_age <- cox.zph(cox_age)
print(cox_zph_age)
plot(cox_zph_age)

veteran$martingale_age  <- residuals(cox_age, type = "martingale")
veteran$linear_pred_age <- predict(cox_age, type = "lp")

ggplot(veteran, aes(x = age, y = martingale_age, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  labs(title = "Martingale Residuals vs Age",
       x = "Age (years)", y = "Martingale Residual", color = "Treatment Group") +
  theme_minimal(base_size = 13)

# trt + prior
cox_prior <- coxph(Surv(time, status) ~ trt + prior, data = veteran)
summary(cox_prior)

cox_zph_prior <- cox.zph(cox_prior)
print(cox_zph_prior)
plot(cox_zph_prior)

veteran$martingale_prior  <- residuals(cox_prior, type = "martingale")
veteran$linear_pred_prior <- predict(cox_prior, type = "lp")

ggplot(veteran, aes(x = prior, y = martingale_prior, color = trt)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  labs(title = "Martingale Residuals vs Prior Treatment",
       x = "Prior Treatment", y = "Martingale Residual", color = "Treatment Group") +
  theme_minimal(base_size = 13)

# trt + age + prior
cox_multivariable <- coxph(Surv(time, status) ~ trt + age + prior, data = veteran)
summary(cox_multivariable)

cox_zph_multivariable <- cox.zph(cox_multivariable)
print(cox_zph_multivariable)
plot(cox_zph_multivariable)

veteran$martingale_multivariable  <- residuals(cox_multivariable, type = "martingale")
veteran$linear_pred_multivariable <- predict(cox_multivariable, type = "lp")

ggplot(veteran, aes(x = linear_pred_multivariable, y = martingale_multivariable, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  labs(title = "Martingale Residuals vs Linear Predictor (Multivariable Model)",
       x = "Linear Predictor (Log Hazard)", y = "Martingale Residual",
       color = "Treatment Group") +
  theme_minimal(base_size = 13)


# ---- 6. Functional Form Check: Age ------------------------------------------

cox_age_only <- coxph(Surv(time, status) ~ age, data = veteran)
ggcoxfunctional(cox_age_only, data = veteran)

cox_age_for_form_check <- coxph(Surv(time, status) ~ age, data = veteran)
ggcoxfunctional(
  fit       = cox_age_for_form_check,
  data      = veteran,
  point.col = "#E74C3C",
  ggtheme   = theme_minimal(base_size = 13),
  title     = "Functional Form Check: Age"
)
# U-shaped curve -> non-linear -> use splines


# ---- 7. Spline Model & Comparison --------------------------------------------

cox_spline <- coxph(Surv(time, status) ~ ns(age, df = 3) + trt + prior, data = veteran)
summary(cox_spline)

cox.zph(cox_spline)
plot(cox.zph(cox_spline))

veteran$martingale_spline  <- residuals(cox_spline, type = "martingale")
veteran$linear_pred_spline <- predict(cox_spline, type = "lp")

ggplot(veteran, aes(x = linear_pred_spline, y = martingale_spline, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  labs(title = "Martingale Residuals vs Linear Predictor (Spline Model)",
       x = "Linear Predictor (Log Hazard)", y = "Martingale Residual",
       color = "Treatment Group") +
  theme_minimal(base_size = 13)

# LRT: linear age vs spline age (chi2 = 7.91, df = 2, p = 0.019 -> spline wins)
anova(cox_multivariable, cox_spline, test = "LRT")

termplot(cox_spline, se = TRUE, col.term = "blue")


# ---- 8. Extended Model: + Karnofsky + Cell Type ------------------------------

cox_karno_and_cell <- coxph(
  Surv(time, status) ~ ns(age, df = 3) + trt + prior + karno + celltype,
  data = veteran
)
summary(cox_karno_and_cell)
# Concordance: 0.735 (vs 0.559 for spline-only model)

forest_model(cox_karno_and_cell)

# KM by Karnofsky group
veteran$karno_group <- cut(veteran$karno,
                           breaks = c(0, 50, 70, 100),
                           labels = c("Low (≤50)", "Medium (51–70)", "High (>70)"),
                           right  = TRUE)

km_karno <- survfit(Surv(time, status) ~ karno_group, data = veteran)

ggsurvplot(km_karno,
           data         = veteran,
           risk.table   = TRUE,
           pval         = TRUE,
           conf.int     = FALSE,
           legend.title = "Karnofsky Group",
           legend.labs  = c("Low (≤50)", "Medium (51–70)", "High (>70)"),
           palette      = c("firebrick", "steelblue", "darkgreen"),
           xlab         = "Time (days)",
           ylab         = "Survival Probability",
           title        = "Kaplan-Meier Survival by Karnofsky Score Group")

# KM by cell type
km_cell <- survfit(Surv(time, status) ~ celltype, data = veteran)

ggsurvplot(km_cell,
           data         = veteran,
           risk.table   = TRUE,
           pval         = TRUE,
           conf.int     = FALSE,
           legend.title = "Cell Type",
           legend.labs  = levels(factor(veteran$celltype)),
           palette      = "Dark2",
           xlab         = "Time (days)",
           ylab         = "Survival Probability",
           title        = "Kaplan-Meier Survival by Cell Type")

# Combined side-by-side
plot_cell  <- ggsurvplot(km_cell,  data = veteran, risk.table = TRUE, pval = TRUE,
                         conf.int = FALSE, legend.title = "Cell Type",
                         palette = "Dark2", title = "Survival by Cell Type")
plot_karno <- ggsurvplot(km_karno, data = veteran, risk.table = TRUE, pval = TRUE,
                         conf.int = FALSE, legend.title = "Karnofsky Score",
                         palette = c("firebrick", "steelblue", "darkgreen"),
                         title = "Survival by Karnofsky Group")
arrange_ggsurvplots(list(plot_cell, plot_karno), ncol = 2, nrow = 1)


# ---- 9. Final Model Diagnostics ----------------------------------------------

cox_zph_final <- cox.zph(cox_karno_and_cell)
print(cox_zph_final)
plot(cox_zph_final)

veteran$martingale_final  <- residuals(cox_karno_and_cell, type = "martingale")
veteran$linear_pred_final <- predict(cox_karno_and_cell, type = "lp")

ggplot(veteran, aes(x = linear_pred_final, y = martingale_final, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  labs(title = "Martingale Residuals vs Linear Predictor (Final Model)",
       x = "Linear Predictor (Log Hazard)", y = "Martingale Residual",
       color = "Treatment Group") +
  theme_minimal(base_size = 13)

# Functional form check: Karnofsky (linear is fine — no spline needed)
cox_karno_for_form_check <- coxph(Surv(time, status) ~ karno, data = veteran)
ggcoxfunctional(
  fit       = cox_karno_for_form_check,
  data      = veteran,
  point.col = "#E74C3C",
  ggtheme   = theme_minimal(base_size = 13),
  title     = "Functional Form Check: Karnofsky Score"
)


# ---- 10. Interaction Models --------------------------------------------------

# trt x age (spline)
cox_interaction <- coxph(Surv(time, status) ~ trt * ns(age, df = 3) + prior,
                         data = veteran)
summary(cox_interaction)
cox.zph(cox_interaction)
plot(cox.zph(cox_interaction))

veteran$martingale_interaction  <- residuals(cox_interaction, type = "martingale")
veteran$linear_pred_interaction <- predict(cox_interaction, type = "lp")
veteran$deviance_interaction    <- residuals(cox_interaction, type = "deviance")

ggplot(veteran, aes(x = linear_pred_interaction, y = martingale_interaction, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  labs(title = "Martingale Residuals (Interaction: trt x age)",
       x = "Linear Predictor", y = "Martingale Residual", color = "Treatment Group") +
  theme_minimal(base_size = 13)

ggplot(veteran, aes(x = linear_pred_interaction, y = deviance_interaction, color = trt)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  labs(title = "Deviance Residuals (Interaction: trt x age)",
       x = "Linear Predictor", y = "Deviance Residual", color = "Treatment Group") +
  theme_minimal(base_size = 13)

# trt x prior
cox_interaction_prior <- coxph(Surv(time, status) ~ trt * prior + ns(age, df = 3),
                               data = veteran)
summary(cox_interaction_prior)
cox.zph(cox_interaction_prior)
plot(cox.zph(cox_interaction_prior))

veteran$martingale_interaction_prior  <- residuals(cox_interaction_prior, type = "martingale")
veteran$linear_pred_interaction_prior <- predict(cox_interaction_prior, type = "lp")
veteran$deviance_interaction_prior    <- residuals(cox_interaction_prior, type = "deviance")

ggplot(veteran, aes(x = linear_pred_interaction_prior, y = martingale_interaction_prior, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  labs(title = "Martingale Residuals (Interaction: trt x prior)",
       x = "Linear Predictor", y = "Martingale Residual", color = "Treatment Group") +
  theme_minimal(base_size = 13)

ggplot(veteran, aes(x = linear_pred_interaction_prior, y = deviance_interaction_prior, color = trt)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  labs(title = "Deviance Residuals (Interaction: trt x prior)",
       x = "Linear Predictor", y = "Deviance Residual", color = "Treatment Group") +
  theme_minimal(base_size = 13)


# ---- 11. Overview: Summary Visualisations (veteran2) -------------------------
# Uses veteran2 — fully labelled dataset (Standard/Experimental, all factors set)

# KM medians
km_all <- survfit(Surv(time, status) ~ 1,   data = veteran2)
km_by  <- survfit(Surv(time, status) ~ trt, data = veteran2)
surv_median(km_all)
surv_median(km_by)

lr   <- survdiff(Surv(time, status) ~ trt, data = veteran2)
lr_p <- signif(1 - pchisq(lr$chisq, length(lr$n) - 1), 3)

# KM by treatment
ggsurvplot(
  km_by,
  data         = veteran2,
  risk.table   = TRUE,
  pval         = TRUE,
  conf.int     = TRUE,
  palette      = c("#2E86AB", "#E84855"),
  legend.title = "Treatment",
  xlab         = "Time (days)",
  ylab         = "Survival probability",
  title        = "Kaplan-Meier Survival Curves by Treatment Arm",
  ggtheme      = theme_minimal(base_size = 13)
)

# KM by Karnofsky group
veteran2$karno_group <- cut(veteran2$karno,
                            breaks = c(0, 50, 70, 100),
                            labels = c("Low (≤50)", "Medium (51–70)", "High (>70)"),
                            right  = TRUE)

km_karno2 <- survfit(Surv(time, status) ~ karno_group, data = veteran2)

ggsurvplot(
  km_karno2,
  data         = veteran2,
  risk.table   = TRUE,
  pval         = TRUE,
  conf.int     = FALSE,
  palette      = c("#E84855", "#F4A259", "#2E86AB"),
  legend.title = "Karnofsky group",
  legend.labs  = c("Low (≤50)", "Medium (51–70)", "High (>70)"),
  xlab         = "Time (days)",
  ylab         = "Survival probability",
  title        = "Kaplan-Meier Survival Curves by Karnofsky Score Group",
  ggtheme      = theme_minimal(base_size = 13)
)

# Final Cox model forest plot
cox_final <- coxph(
  Surv(time, status) ~ trt + prior + ns(age, df = 3) + karno + celltype,
  data = veteran2
)

ggforest(
  cox_final,
  data      = veteran2,
  main      = "Adjusted Hazard Ratios — Final Cox PH Model",
  fontsize  = 0.9,
  noDigits  = 2,
  refLabel  = "Reference"
)

# Aalen additive model
aa_fit <- aareg(
  Surv(time, status) ~ trt + celltype + karno + age + prior,
  data = veteran2
)
aa_fit

autoplot(aa_fit) +
  theme_minimal(base_size = 12) +
  labs(title = "Aalen Additive Model: Cumulative Regression Functions")
