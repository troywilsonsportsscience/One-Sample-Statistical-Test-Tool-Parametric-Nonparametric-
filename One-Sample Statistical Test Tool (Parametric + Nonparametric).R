# ===============================================================
# One-Sample Statistical Test Tool (Parametric + Nonparametric)
# ===============================================================
# 
# This R script performs a fully automated one-sample hypothesis test
# for continuous numeric data. It is designed for flexible reuse across
# datasets and projects in data science, human performance, and research.
#
# --------------------------
# What the Script Does:
# --------------------------
# 1. Loads or simulates data with a numeric test metric.
# 2. Checks for normality using the Shapiro-Wilk test and visual inspection.
# 3. If data are non-normal:
#    - Applies three transformations: log, square root, and inverse.
#    - Selects the transformation that most improves normality (highest p-value).
#    - If none succeed, defaults to a Wilcoxon signed-rank test.
# 4. Executes the appropriate one-sample test:
#    - t-test (raw or transformed data)
#    - Wilcoxon signed-rank test (nonparametric)
# 5. Reports results in APA-style, including test direction and significance.
#
# -------------------------------
# Key Parameters (User-Editable):
# -------------------------------
# - `mu_reference`: the population value to compare against.
#     → If not set manually, the script will auto-calculate it as the sample mean.
# - `tail_type`: "two.sided", "greater", or "less" to specify the hypothesis direction.
# - `alpha_level`: significance threshold (default is 0.05).
# - `metric_label`: a friendly name for the test variable used in reporting.
#
# ---------------------------
# Transformations Included:
# ---------------------------
# - Log transformation: log(x)
# - Square root: sqrt(x)
# - Inverse: 1 / x
#
# -------------------------
# Visualization Features:
# -------------------------
# - Histograms with overlaid density curves for raw and transformed data
# - Faceted plots for comparison of transformation effects
#
# --------------------------
# Output Includes:
# --------------------------
# - Console-based APA-style test summary
# - Tidy result tables for parametric and nonparametric results
# - Optional back-transformed summary for transformed tests
#
# --------------------------
# Dependencies:
# --------------------------
# - tidyverse, car, broom, rstatix, BSDA, ggplot2, tibble, dplyr, janitor

# -------------------------------
# Helper: Install and Load Packages
# -------------------------------
# These packages are required for data manipulation, visualization, assumption checks, and APA-style reporting.

required_packages <- c("tidyverse", "car", "broom", "rstatix", "janitor")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(required_packages, install_if_missing))

# ===============================================================

# ----# 0. Set-up----
# Define target significance level
alpha_level <- 0.05  # Common value for hypothesis testing.

# Define directionality of the test: "two.sided", "less", or "greater"
tail_type <- "two.sided"  # change to "greater" or "less" if testing a directional hypothesis
# Note: For a one-sample t-test, we typically use "two.sided" unless there's a strong rationale for a directional test.

# Define metric you will be using so you can easily update reporting later.
metric_label <- "jump height"  # or "vertical displacement", "sprint time", etc.

# ----# 1. Real Data Input----
# Load or assign your real dataset here
# Example structure: your_data <- read_csv("your_file.csv") OR assign directly
# Your data frame must include a numeric variable called `test_metric`

# Uncomment and edit:
# path <- file.choose()  # Select CSV file
# df <- read_csv(path, show_col_types = FALSE) %>% clean_names() # Convert to snake_case for ease of use
# test_data <- df %>%
#   select(subject_id, test_metric)  # or rename columns to match this format

# Define population/reference value (mu) to compare against
mu_reference <- 55  # Edit this value as needed for your analysis

# If mu_reference is not predefined, calculate it from the data
# Defaulting to the mean of test_metric
if (!exists("mu_reference")) {
  mu_reference <- mean(test_data$test_metric, na.rm = TRUE)
  cat(sprintf("ℹ️  Auto-calculated mu_reference from data: %.2f\n", mu_reference))
}

# If you paste in a data frame, make sure it looks something like:
# test_data <- tibble(
#   subject_id = c(1, 2, 3, ...),
#   test_metric = c(57.2, 58.1, 55.6, ...)
# )

# ---- 1B. Simulated Example Data (optional)----

# Simulated Normal example: test metric from 30 subjects
# set.seed(123)
# test_data <- tibble(
#   subject_id = 1:30,
#   test_metric = rnorm(30, mean = 57, sd = 5)
# )

# # Simulated Non-Normal example: positively skewed or with outliers
# set.seed(202)
# test_data <- tibble(
#   subject_id = 1:30,
#   test_metric = rexp(30, rate = 1/55)  # Exponential distribution (positively skewed)
# )

# Simulated Non-Normal example: heavy right-tail outliers
set.seed(404)
test_data <- tibble(
  subject_id = 1:30,
  test_metric = c(rnorm(25, mean = 55, sd = 4), rep(100, 5))  # heavy right-tail outliers
)


# ----2A. Check Normality of the Sample----

# Visual check: Histogram with normal curve
ggplot(test_data, aes(x = test_metric)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(test_data$test_metric), 
                                         sd = sd(test_data$test_metric)), 
                col = "red", linetype = "dashed") +
  labs(title = "Histogram of Test Metric", x = "Test Metric", y = "Density")

# Shapiro-Wilk test for normality (recommended when n < 50)
shapiro_test <- shapiro.test(test_data$test_metric)
print(shapiro_test)

# Branch based on normality result
if (shapiro_test$p.value < 0.05) {
  cat("❗ Normality violated (Shapiro-Wilk p < 0.05). Trying transformations in Section 2B...\n\n")
  
  # Flow continues in Section 2B
  
} else {
  cat("✅ Normality confirmed (Shapiro-Wilk p ≥ 0.05). Proceeding with raw data.\n")
  cat("➡️ Move to Section 3 for results reporting. Sections 2B (transformations) and 2C (nonparametric) will be skipped.\n\n")
  
  # Prep for analysis
  analysis_var <- test_data$test_metric
  mu_test <- mu_reference
  label <- "raw"
  t_result <- t.test(analysis_var, mu = mu_test, alternative = tail_type)
  mean_val <- mean(analysis_var)
  sd_val   <- sd(analysis_var)
  test_label <- "One-sample t-test on raw data"
  analysis_method <- "t_standard"
}
# -----2B. If normality violated try transformations for parametric testing. Otherwise move to Wilcoxon Signed-Rank Test (Nonparametric)-----------

# Apply three transformations
test_data <- test_data %>%
  mutate(
    log_metric  = log(test_metric),
    sqrt_metric = sqrt(test_metric),
    inv_metric  = 1 / test_metric
  )

# Check normality of each
log_p  <- shapiro.test(test_data$log_metric)$p.value
sqrt_p <- shapiro.test(test_data$sqrt_metric)$p.value
inv_p  <- shapiro.test(test_data$inv_metric)$p.value

cat("🔎 Shapiro-Wilk p-values for transformations:\n",
    sprintf("- Log: %.4f\n- Sqrt: %.4f\n- Inverse: %.4f\n", log_p, sqrt_p, inv_p))

# --- Visualization of Original vs Transformed Distributions ---
plot_data <- test_data %>%
  select(test_metric, log_metric, sqrt_metric, inv_metric) %>%
  pivot_longer(cols = everything(), names_to = "transform_type", values_to = "value") %>%
  mutate(
    transform_type = dplyr::recode(transform_type,
                                   "test_metric" = "Raw",
                                   "log_metric" = "Log",
                                   "sqrt_metric" = "Square Root",
                                   "inv_metric" = "Inverse"
    )
  )

ggplot(plot_data, aes(x = value)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "skyblue", color = "black") +
  geom_density(col = "red", linetype = "dashed", size = 1) +
  facet_wrap(~transform_type, scales = "free", ncol = 2) +
  labs(title = "Histogram of Raw and Transformed Metrics",
       x = "Value",
       y = "Density") +
  theme_minimal()

# ---- Transformation Selection ----

# Create named list of p-values
transform_results <- list(
  log = list(p = log_p, var = test_data$log_metric, mu = log(mu_reference), label = "log-transformed"),
  sqrt = list(p = sqrt_p, var = test_data$sqrt_metric, mu = sqrt(mu_reference), label = "square root-transformed"),
  inv = list(p = inv_p, var = test_data$inv_metric, mu = 1 / mu_reference, label = "inverse-transformed")
)

# Filter for transformations that pass normality, then pick the best (highest p-value)
valid_transforms <- Filter(function(x) x$p > alpha_level, transform_results)

if (length(valid_transforms) > 0) {
  # Select the transformation with the highest p-value
  best_transform <- valid_transforms[[which.max(sapply(valid_transforms, function(x) x$p))]]
  
  analysis_var <- best_transform$var
  mu_test <- best_transform$mu
  label <- best_transform$label
  
  cat(sprintf("✔️ Normality restored using %s data. Proceeding with one-sample t-test.\n", label))
  
  t_result <- t.test(analysis_var, mu = mu_test, alternative = tail_type)
  mean_val <- mean(analysis_var)
  sd_val   <- sd(analysis_var)
  test_label <- paste("One-sample t-test on", label, "data")
  analysis_method <- "t_transformed"
  
} else {
  cat("❌ No transformation restored normality. Will use Wilcoxon test in Section 2C.\n")
  analysis_method <- "nonparametric"
}


# -----2C. If Normality Assumption is Violated Use Wilcoxon Signed-Rank Test (Nonparametric)--------

if (analysis_method == "nonparametric") {
  wilcox_result <- wilcox.test(test_data$test_metric, mu = mu_reference, alternative = tail_type)
  print(wilcox_result)
  
  # Reporting Guidance:
  # "A Wilcoxon signed-rank test indicated that the median test metric was [significantly/not significantly] 
  # different from mu_reference, V = XX, p = XX."
}
# If Warning - Warning message: In wilcox.test.default(test_data$test_metric, mu = mu_reference: cannot compute exact p-value with ties
# Indicates that there are tied ranks in the data, so the model swaps to another method of p-value calculation.

# -----# 3. Run One-Sample t-Test--------

# Only runs if parametric method was selected (normal or transformed)
if (analysis_method %in% c("t_standard", "t_transformed")) {
  
  # H0: Mean of test metric = mu_reference
  # H1: Mean of test metric ≠ mu_reference (or one-tailed, depending on tail_type)
  
  # Re-run t-test if not already executed
  if (!exists("t_result")) {
    t_result <- t.test(analysis_var, mu = mu_test, alternative = tail_type)
    mean_val <- mean(analysis_var)
    sd_val   <- sd(analysis_var)
  }
  
  print(t_result)
  
  # Tidy output for APA-style reporting
  t_result_tidy <- tidy(t_result)
  
  # Create readable summary table
  t_test_summary <- tibble::tibble(
    Method = t_result_tidy$method,
    `Test Type` = ifelse(analysis_method == "t_transformed", "Transformed", "Standard"),
    `Test Direction` = tail_type,
    `Sample Mean` = round(mean_val, 2),
    `SD` = round(sd_val, 2),
    `t(df)` = paste0("t(", round(t_result_tidy$parameter, 0), ") = ", round(t_result_tidy$statistic, 2)),
    `p-value` = round(t_result_tidy$p.value, 4),
    `95% CI` = paste0("[", round(t_result_tidy$conf.low, 2), ", ", round(t_result_tidy$conf.high, 2), "]")
  )
  
  print(t_test_summary)
  
}

# -----4A. APA-style Reporting (Parametric Results)-------------

if (analysis_method %in% c("t_standard", "t_transformed")) {
  
  # Optional: define a readable name for the tested metric
  # If not already defined in Section 1:
  # metric_label <- "jump height"
  
  # Determine significance wording
  sig_label <- ifelse(t_result$p.value < alpha_level, "significantly", "not significantly")
  
  # Directional comparison label
  direction_label <- switch(
    tail_type,
    "two.sided" = "different from",
    "greater"   = "greater than",
    "less"      = "less than"
  )
  
  # APA-style output
  cat(sprintf(
    "A one-sample t-test was conducted to compare %s to the reference value of %.2f (α = %.2f).\n",
    metric_label, mu_reference, alpha_level
  ))
  
  if (analysis_method == "t_transformed") {
    cat(sprintf(
      "The data were %s to meet normality assumptions prior to analysis.\n",
      label
    ))
  }
  
  cat(sprintf(
    "The mean %s (M = %.2f, SD = %.2f) was %s %.2f, t(%d) = %.2f, p = %.3f, 95%% CI [%.2f, %.2f].\n",
    metric_label,
    mean_val,
    sd_val,
    sig_label,
    mu_reference,
    round(t_result$parameter, 0),
    round(t_result$statistic, 2),
    round(t_result$p.value, 3),
    round(t_result$conf.int[1], 2),
    round(t_result$conf.int[2], 2)
  ))
}

# -----4B. APA-style Reporting (Wilcoxon Nonparametric Results)-------------

if (analysis_method == "nonparametric") {
  
  # Determine significance wording
  sig_label <- ifelse(wilcox_result$p.value < alpha_level, "significantly", "not significantly")
  
  # Directional label for APA-style language
  direction_label <- switch(
    tail_type,
    "two.sided" = "different from",
    "greater"   = "greater than",
    "less"      = "less than"
  )
  
  # Manually create tidy-style summary
  wilcox_summary <- tibble::tibble(
    Method = "Wilcoxon Signed-Rank Test",
    `Test Direction` = tail_type,
    `Reference Value` = mu_reference,
    `Statistic (V)` = wilcox_result$statistic,
    `p-value` = round(wilcox_result$p.value, 4),
    `Alpha Level` = alpha_level,
    `Significant?` = ifelse(wilcox_result$p.value < alpha_level, "Yes", "No")
  )
  
  print(wilcox_summary)
  
  # APA-style reporting
  cat(sprintf(
    "A Wilcoxon signed-rank test was conducted to assess whether %s was %s the reference value of %.2f (α = %.2f).\n",
    metric_label,
    direction_label,
    mu_reference,
    alpha_level
  ))
  
  cat("This nonparametric test was used because the normality assumption could not be satisfied, even after transformation.\n")
  
  cat(sprintf(
    "The result showed that the median %s was %s %.2f, V = %.2f, p = %.3f.\n",
    metric_label,
    sig_label,
    mu_reference,
    wilcox_result$statistic,
    wilcox_result$p.value
  ))
}
# End of Script
