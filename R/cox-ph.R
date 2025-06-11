# Load required libraries
library(survival) # Core survival analysis tools
library(survminer)  # Visualization and diagnostics
library(dplyr)  # Data wrangling
library(ggplot2)  # Plotting
library(broom)  # Tidying model output
library(splines) # Splines for flexible modeling

# Load the dataset
veteran <- survival::veteran %>%
  mutate(trt = factor(trt,
                      levels = c(1, 2),
                      labels = c("Control", "Intervention")))

# Peek at the structure and summary
# glimpse(veteran)
# summary(veteran)

head(veteran)

# Create the survival object
surv_obj <- Surv(time = veteran$time, event = veteran$status)

# Fit Kaplan-Meier curve stratified by treatment
km_fit <- survfit(surv_obj ~ trt, data = veteran)

# View basic survival summary
# summary(km_fit)

# View basic survival summary
summary(km_fit)

# Plot the KM curve
ggsurvplot(
  km_fit,
  data = veteran,
  conf.int = TRUE,
  pval = TRUE,
  pval.size = 4,
  risk.table = TRUE,
  risk.table.title = "Number at Risk",
  risk.table.col = "strata",
  risk.table.y.text.col = TRUE,
  risk.table.fontsize = 4,
  legend.title = "Treatment Group",
  legend.labs = c("Control", "Intervention"),
  xlab = "Time (days)",
  ylab = "Survival Probability",
  palette = c("#A569BD", "#45B39D"),
  title = "Kaplan-Meier Survival Curve: VA Lung Cancer Trial",
  ggtheme = theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2)
    ),
  tables.theme = theme(  # Trying to tweak the risk table
    axis.text.y = element_text(size = 10, margin = margin(r = 12), lineheight = 2)
  )
)

# Perform the log-rank test
log_rank_test <- survdiff(surv_obj ~ trt, data = veteran)
log_rank_test

# Extract the p-value
1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1)

# Fit the Cox Proportional Hazards model
cox_model <- coxph(surv_obj ~ trt, data = veteran) # Where trt is a binary treatment variable with levels Control and Intervention

# Summarize the model
summary(cox_model)

# Use ggadjustedcurves for plotting adjusted survival curves directly from the Cox model
ggadjustedcurves(
  cox_model,                     # The fitted Cox model
  data = veteran,                # The original data used to fit the model
  variable = "trt",              # The variable to plot adjusted curves for
  method = "average",            # Method for adjusting (average effect of other covariates)
  conf.int = TRUE,
  legend.labs = c("Control", "Intervention"),
  legend.title = "Treatment Group",
  risk.table = TRUE,             # Include the risk table
  risk.table.title = "Number at Risk",
  xlab = "Time (days)",
  ylab = "Adjusted Survival Probability",
  title = "Cox Model Adjusted Survival Curves",
  palette = c("#A569BD", "#45B39D"),
  ggtheme = theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.2)
    ),
  tables.theme = theme(
    axis.text.y = element_text(size = 10, margin = margin(r = 12), lineheight = 2)
  )
)

# Test the PH assumption
cox.zph(cox_model)

plot(cox.zph(cox_model))

# Generate Martingale residuals
martingale_resid <- residuals(cox_model, type = "martingale")

# Add residuals to your data
veteran$martingale <- martingale_resid

# Plot residuals against a continuous covariate (e.g., time or another variable)
# Here we plot against the linear predictor for simplicity
veteran$linear_pred <- predict(cox_model, type = "lp")

ggplot(veteran, aes(x = linear_pred, y = martingale, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Martingale Residuals vs Linear Predictor",
    x = "Linear Predictor (Log Hazard)",
    y = "Martingale Residual",
    color = "Treatment Group"
  ) +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

cox_age <- coxph(Surv(time, status) ~ trt + age, data = veteran)
summary(cox_age)

# Test proportional hazards assumption for the model including age
cox_zph_age <- cox.zph(cox_age)
print(cox_zph_age)

# Visual check
plot(cox_zph_age)

# Generate Martingale residuals
veteran$martingale_age <- residuals(cox_age, type = "martingale")

# Generate linear predictor
veteran$linear_pred_age <- predict(cox_age, type = "lp")

# Plot residuals against age
ggplot(veteran, aes(x = age, y = martingale_age, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Martingale Residuals vs Age",
    x = "Age (years)",
    y = "Martingale Residual",
    color = "Treatment Group"
  ) +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Fit a univariable Cox model to assess functional form of age
cox_age_only <- coxph(Surv(time, status) ~ age, data = veteran)

# Check the functional form visually using Martingale residuals
ggcoxfunctional(cox_age_only, data = veteran)

# Fit the Cox model including prior treatment
cox_prior <- coxph(Surv(time, status) ~ trt + prior, data = veteran)
summary(cox_prior)

# Test proportional hazards assumption for the model including prior treatment
cox_zph_prior <- cox.zph(cox_prior)
print(cox_zph_prior)

# Visual check
plot(cox_zph_prior)

# Generate Martingale residuals
veteran$martingale_prior <- residuals(cox_prior, type = "martingale")
# Generate linear predictor
veteran$linear_pred_prior <- predict(cox_prior, type = "lp")
# Plot residuals against prior treatment
ggplot(veteran, aes(x = prior, y = martingale_prior, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Martingale Residuals vs Prior Treatment",
    x = "Prior Treatment (0 = No, 1 = Yes)",
    y = "Martingale Residual",
    color = "Treatment Group"
  ) +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.6) # Adding jitter to avoid overlap

# Fit the multivariable Cox model including age and prior treatment
cox_multivariable <- coxph(Surv(time, status) ~ trt + age + prior, data = veteran)
summary(cox_multivariable)

# Test proportional hazards assumption for the multivariable model
cox_zph_multivariable <- cox.zph(cox_multivariable)
print(cox_zph_multivariable)

# Visual check
plot(cox_zph_multivariable)

# Generate Martingale residuals
veteran$martingale_multivariable <- residuals(cox_multivariable, type = "martingale")
# Generate linear predictor
veteran$linear_pred_multivariable <- predict(cox_multivariable, type = "lp")
# Plot residuals against linear predictor
ggplot(veteran, aes(x = linear_pred_multivariable, y = martingale_multivariable, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Martingale Residuals vs Linear Predictor\n(Multivariable Model)",
    x = "Linear Predictor (Log Hazard)",
    y = "Martingale Residual",
    color = "Treatment Group"
  ) +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Refit a univariable model with age (again, for functional form check)
cox_age_for_form_check <- coxph(Surv(time, status) ~ age, data = veteran)

# Plot functional form using ggcoxfunctional
ggcoxfunctional(
  fit = cox_age_for_form_check,
  data = veteran,
  point.col = "#E74C3C",
  ggtheme = theme_minimal(base_size = 13),
  title = "Functional Form Check: Age"
)

# Fit a Cox model with restricted cubic splines for age
cox_spline <- coxph(Surv(time, status) ~ ns(age, df = 3) + trt + prior, data = veteran)
summary(cox_spline)

# Test proportional hazards assumption for the spline model
cox.zph(cox_spline)

# Visual check of Schoenfeld residuals for spline model
plot(cox.zph(cox_spline))

# Residuals vs linear predictor
veteran$martingale_spline <- residuals(cox_spline, type = "martingale")
veteran$linear_pred_spline <- predict(cox_spline, type = "lp")

ggplot(veteran, aes(x = linear_pred_spline, y = martingale_spline, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Martingale Residuals vs Linear Predictor\n(Spline Model)",
    x = "Linear Predictor (Log Hazard)",
    y = "Martingale Residual",
    color = "Treatment Group"
  ) +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Compare linear age model vs spline-transformed age model
anova(cox_multivariable, cox_spline, test = "LRT")

termplot(cox_spline, se = TRUE, col.term = "blue")

# Fit a cox model with interaction between treatment and age
cox_interaction <- coxph(Surv(time, status) ~ trt * ns(age, df = 3) + prior, data = veteran)
summary(cox_interaction)

# Test proportional hazards assumption for the interaction model
cox.zph(cox_interaction)

# Visual check of Schoenfeld residuals for interaction model
plot(cox.zph(cox_interaction))

# Residuals vs linear predictor
veteran$martingale_interaction <- residuals(cox_interaction, type = "martingale")
veteran$linear_pred_interaction <- predict(cox_interaction, type = "lp")
ggplot(veteran, aes(x = linear_pred_interaction, y = martingale_interaction, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Martingale Residuals vs Linear Predictor\n(Interaction Model)",
    x = "Linear Predictor (Log Hazard)",
    y = "Martingale Residual",
    color = "Treatment Group"
  ) +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Deviance residuals
veteran$deviance_interaction <- residuals(cox_interaction, type = "deviance")
# Plot deviance residuals
ggplot(veteran, aes(x = linear_pred_interaction, y = deviance_interaction, color = trt)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.6) + # Adding jitter to avoid overlap
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Deviance Residuals vs Linear Predictor\n(Interaction Model)",
    x = "Linear Predictor (Log Hazard)",
    y = "Deviance Residual",
    color = "Treatment Group"
  ) +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Fit a cox model with interaction between treatment and prior therapy
cox_interaction_prior <- coxph(Surv(time, status) ~ trt * prior + ns(age, df = 3), data = veteran)
summary(cox_interaction_prior)

# Test proportional hazards assumption for the interaction with prior therapy
cox.zph(cox_interaction_prior)

# Visual check of Schoenfeld residuals for interaction with prior therapy
plot(cox.zph(cox_interaction_prior))

# Residuals vs linear predictor
veteran$martingale_interaction_prior <- residuals(cox_interaction_prior, type = "martingale")
veteran$linear_pred_interaction_prior <- predict(cox_interaction_prior, type = "lp")

ggplot(veteran, aes(x = linear_pred_interaction_prior, y = martingale_interaction_prior, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Martingale Residuals vs Linear Predictor\n(Interaction with Prior Therapy)",
    x = "Linear Predictor (Log Hazard)",
    y = "Martingale Residual",
    color = "Treatment Group"
  ) +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Deviance residuals
veteran$deviance_interaction_prior <- residuals(cox_interaction_prior, type = "deviance")

# Plot deviance residuals
ggplot(veteran, aes(x = linear_pred_interaction_prior, y = deviance_interaction_prior, color = trt)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.6) + # Adding jitter to avoid overlap
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Deviance Residuals vs Linear Predictor\n(Interaction with Prior Therapy)",
    x = "Linear Predictor (Log Hazard)",
    y = "Deviance Residual",
    color = "Treatment Group"
  ) +
  scale_color_manual(values = c("Control" = "#A569BD", "Intervention" = "#45B39D")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Potential next steps:

# 1. AIC and BIC for model comparison
AIC(cox_model, cox_multivariable, cox_spline, cox_interaction, cox_interaction_prior)
BIC(cox_model, cox_multivariable, cox_spline, cox_interaction, cox_interaction_prior)

# 2. Influence diagnostics (e.g., dfbeta, dffit)
# Example: DFBETA plots for the main multivariable model
ggcoxdiagnostics(cox_multivariable, type = "dfbeta", linear.predictions = FALSE, ggtheme = theme_minimal())

# 3. Cox-Snell residuals for overall model fit
# Calculate Cox-Snell residuals
veteran$coxsnell <- -log(survfit(cox_multivariable)$surv)
# Plot Nelson-Aalen cumulative hazard vs Cox-Snell residuals
library(survival)
na_fit <- survfit(Surv(veteran$coxsnell, veteran$status) ~ 1)
plot(na_fit$time, -log(na_fit$surv), type = "s",
     xlab = "Cox-Snell Residuals", ylab = "Cumulative Hazard",
     main = "Cox-Snell Residuals Plot")
abline(0, 1, col = "red", lty = 2)

# 4. Baseline hazard estimation and plotting
basehaz_df <- basehaz(cox_multivariable, centered = FALSE)
plot(basehaz_df$time, basehaz_df$hazard, type = "l",
     xlab = "Time", ylab = "Baseline Cumulative Hazard",
     main = "Baseline Cumulative Hazard Function")

# 5. Use broom package to tidy model output for reporting
# Example: Tidy summary of the multivariable model
library(broom)
tidy_cox <- tidy(cox_multivariable, exponentiate = TRUE, conf.int = TRUE)
print(tidy_cox)

# 7. Session info for reproducibility
sessionInfo()

# Additional Covariates to Include:
# Karnoff score and cell type

# Stepwise inclusion of Karnofsky score and cell type, with diagnostics

# 1. Karnofsky score only
cox_karno <- coxph(Surv(time, status) ~ karno, data = veteran)
summary(cox_karno)

# Schoenfeld residuals (PH assumption)
cox_zph_karno <- cox.zph(cox_karno)
print(cox_zph_karno)
plot(cox_zph_karno)

# Martingale residuals
veteran$martingale_karno <- residuals(cox_karno, type = "martingale")
ggplot(veteran, aes(x = karno, y = martingale_karno)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal() +
  labs(title = "Martingale Residuals vs Karnofsky Score", x = "Karnofsky Score", y = "Martingale Residual")

# Functional form check
ggcoxfunctional(cox_karno, data = veteran)

# 2. Cell type only
cox_celltype <- coxph(Surv(time, status) ~ celltype, data = veteran)
summary(cox_celltype)

# Schoenfeld residuals
cox_zph_celltype <- cox.zph(cox_celltype)
print(cox_zph_celltype)
plot(cox_zph_celltype)

# Martingale residuals
veteran$martingale_celltype <- residuals(cox_celltype, type = "martingale")
ggplot(veteran, aes(x = as.numeric(celltype), y = martingale_celltype, color = celltype)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal() +
  labs(title = "Martingale Residuals vs Cell Type", x = "Cell Type (numeric)", y = "Martingale Residual")

# 3. Karnofsky score + cell type
cox_karno_celltype <- coxph(Surv(time, status) ~ karno + celltype, data = veteran)
summary(cox_karno_celltype)

# Schoenfeld residuals
cox_zph_karno_celltype <- cox.zph(cox_karno_celltype)
print(cox_zph_karno_celltype)
plot(cox_zph_karno_celltype)

# Martingale residuals
veteran$martingale_karno_celltype <- residuals(cox_karno_celltype, type = "martingale")
veteran$linear_pred_karno_celltype <- predict(cox_karno_celltype, type = "lp")
ggplot(veteran, aes(x = linear_pred_karno_celltype, y = martingale_karno_celltype, color = celltype)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal() +
  labs(title = "Martingale Residuals vs Linear Predictor (Karno + Cell Type)",
       x = "Linear Predictor", y = "Martingale Residual")

# 4. Karnofsky score + cell type + age (splines) + prior + trt
cox_full <- coxph(Surv(time, status) ~ trt + prior + ns(age, df = 3) + karno + celltype, data = veteran)
summary(cox_full)

# Schoenfeld residuals
cox_zph_full <- cox.zph(cox_full)
print(cox_zph_full)
plot(cox_zph_full)

# Martingale residuals
veteran$martingale_full <- residuals(cox_full, type = "martingale")
veteran$linear_pred_full <- predict(cox_full, type = "lp")
ggplot(veteran, aes(x = linear_pred_full, y = martingale_full, color = celltype)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  theme_minimal() +
  labs(title = "Martingale Residuals vs Linear Predictor (Full Model)",
       x = "Linear Predictor", y = "Martingale Residual")

# Functional form check for Karnofsky in the full model
ggcoxfunctional(cox_full, data = veteran, covariate = "karno")

