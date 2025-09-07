# ===============================================================
# One-Sample Statistical Test Tool (Parametric + Nonparametric)
# ===============================================================
# Author: Troy Wilson
#
# Purpose:
#   Automated one-sample hypothesis testing for continuous data.
#   Performs assumption checks (Shapiro-Wilk), optional transforms,
#   and runs either a t-test or Wilcoxon signed-rank test.
#   Outputs APA-style results with effect sizes (Cohen's d, Hedges' g).
#
# Documentation:
#   Full usage instructions, examples, and notes are in the README.md
#   file included in this repository.
# ===============================================================

# -------------------------------
# Helper: Install and Load Packages
# -------------------------------
# These packages are required for data manipulation, visualization, assumption checks, and APA-style reporting.

required_packages <- c("tidyverse", "car", "broom", "rstatix", "janitor", "tibble", "ggplot2", "dplyr", "stringr")

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

# ---- Transform Settings (user-configurable) ----
try_transform         <- TRUE                 # attempt transforms for normality, FALSE default to Wilcoxon signed-rank test if normailty fails
allowed_transforms    <- c("log", "sqrt", "inv")  # any subset of: "log", "sqrt", "inv"
show_transform_plot   <- TRUE                 # show the 4-panel plot (raw + transforms)

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
set.seed(202)
test_data <- tibble(
  subject_id = 1:30,
  test_metric = rexp(30, rate = 1/55)  # Exponential distribution (positively skewed)
)

# Simulated Non-Normal example: heavy right-tail outliers
# set.seed(404)
# test_data <- tibble(
#   subject_id = 1:30,
#   test_metric = c(rnorm(25, mean = 55, sd = 4), rep(100, 5))  # heavy right-tail outliers
# )


# ---- 2A. Normality check ----
analysis_method <- NA_character_  # "t_standard", "t_transformed", "nonparametric"

# Visual check
ggplot(test_data, aes(x = test_metric)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm,
                args = list(mean = mean(test_data$test_metric),
                            sd   = sd(test_data$test_metric)),
                col = "red", linetype = "dashed") +
  labs(title = "Histogram of Test Metric", x = "Test Metric", y = "Density")

# Shapiro-Wilk
shapiro_test <- shapiro.test(test_data$test_metric)
print(shapiro_test)

normality_failed <- (shapiro_test$p.value < alpha_level)

if (!normality_failed) {
  cat("Normality confirmed (Shapiro-Wilk p ≥ alpha). Using raw data.\n\n")
  analysis_var <- test_data$test_metric
  mu_test <- mu_reference
  label <- "raw"
  t_result <- t.test(analysis_var, mu = mu_test, alternative = tail_type)
  mean_val <- mean(analysis_var); sd_val <- sd(analysis_var)
  test_label <- "One-sample t-test on raw data"
  analysis_method <- "t_standard"
} else {
  cat("Normality violated (Shapiro-Wilk p < alpha). Assessing next steps...\n\n")
  if (!isTRUE(try_transform)) {
    cat("ℹ️ Transformations disabled (try_transform = FALSE). Proceeding directly to Wilcoxon.\n\n")
    analysis_method <- "nonparametric"
  }
}

# ---- 2B. Transformations (only if normality failed AND toggled on) ----
if (normality_failed && isTRUE(try_transform)) {
  cat("↪ Attempting transformations...\n")
  
  available_transforms <- list(
    log  = function(x) log(x),
    sqrt = function(x) sqrt(x),
    inv  = function(x) 1 / x
  )
  selected_transforms <- available_transforms[names(available_transforms) %in% allowed_transforms]
  
  # Filter invalid transforms
  is_valid <- function(name, x, mu) {
    if (name == "log")  return(all(x > 0,  na.rm = TRUE) && is.finite(log(mu)))
    if (name == "sqrt") return(all(x >= 0, na.rm = TRUE) && is.finite(sqrt(mu)))
    if (name == "inv")  return(!any(abs(x) < 1e-12, na.rm = TRUE) && is.finite(1/mu))
    TRUE
  }
  valid_names <- names(selected_transforms)[purrr::map_lgl(names(selected_transforms),
                                                           ~ is_valid(.x, test_data$test_metric, mu_reference))]
  selected_transforms <- selected_transforms[valid_names]
  
  if (length(selected_transforms) == 0) {
    cat("⚠️ No valid transformations available. Using Wilcoxon.\n\n")
    analysis_method <- "nonparametric"
  } else {
    transformed_cols <- purrr::map_dfc(selected_transforms, ~ .x(test_data$test_metric))
    transformed_cols <- setNames(transformed_cols, paste0(names(selected_transforms), "_metric"))
    
    transformed <- bind_cols(test_data, transformed_cols) |>
      pivot_longer(
        cols = ends_with("_metric"),
        names_to = "transform_type",
        values_to = "transformed_value"
      ) |>
      mutate(transform_type = str_remove(transform_type, "_metric"))
    
    # Shapiro on each transformed series
    transform_results <- transformed |>
      group_by(transform_type) |>
      summarise(p_value = shapiro.test(transformed_value)$p.value, .groups = "drop")
    print(transform_results)
    
    if (isTRUE(show_transform_plot)) {
      cat("↪ Showing histograms for transformed values...\n")
      facet_data <- transformed |>
        left_join(transform_results, by = "transform_type") |>
        mutate(facet_label = paste0(transform_type, " (SW p = ", sprintf("%.3f", p_value), ")"))
      ggplot(facet_data, aes(x = transformed_value)) +
        geom_histogram(aes(y = ..density..), bins = 12, fill = "skyblue", color = "black") +
        geom_density(col = "red", linetype = "dashed", size = 1) +
        facet_wrap(~facet_label, scales = "free", ncol = 2) +
        labs(title = "Distribution of Transformed Metric", x = metric_label, y = "Density") +
        theme_minimal() |>
        print()
    }
    
    viable <- transform_results |>
      filter(p_value > alpha_level) |>
      arrange(desc(p_value))
    
    if (nrow(viable) > 0) {
      best_transform <- viable$transform_type[1]
      
      analysis_var <- transformed |>
        filter(transform_type == best_transform) |>
        pull(transformed_value)
      
      mu_test <- switch(best_transform,
                        log  = log(mu_reference),
                        sqrt = sqrt(mu_reference),
                        inv  = 1 / mu_reference)
      
      label <- case_when(
        best_transform == "log"  ~ "log-transformed",
        best_transform == "sqrt" ~ "square root-transformed",
        best_transform == "inv"  ~ "inverse-transformed",
        TRUE ~ best_transform
      )
      
      cat("✔️ Transformation successful using:", best_transform, "→ t-test on transformed data.\n\n")
      t_result <- t.test(analysis_var, mu = mu_test, alternative = tail_type)
      mean_val <- mean(analysis_var); sd_val <- sd(analysis_var)
      test_label <- paste("One-sample t-test on", label, "data")
      analysis_method <- "t_transformed"
    } else {
      cat("❌ No transformation restored normality. Using Wilcoxon.\n\n")
      analysis_method <- "nonparametric"
    }
  }
}

# ---- 2C. Wilcoxon (nonparametric) if selected path ----
if (analysis_method == "nonparametric") {
  wilcox_result <- wilcox.test(test_data$test_metric, mu = mu_reference, alternative = tail_type)
  print(wilcox_result)
}

  
  # Reporting Guidance:
  # "A Wilcoxon signed-rank test indicated that the median test metric was [significantly/not significantly] 
  # different from mu_reference, V = XX, p = XX."

# If Warning - Warning message: In wilcox.test.default(test_data$test_metric, mu = mu_reference: cannot compute exact p-value with ties
# Indicates that there are tied ranks in the data, so the model swaps to another method of p-value calculation.

# -----# 3A. Run One-Sample t-Test--------

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

# ---3B. Effect size for one-sample t (parametric)---
n <- length(analysis_var)
cohen_d <- (mean_val - mean(mu_test)) / sd_val                 # one-sample d = (M - mu) / SD
J <- 1 - (3 / (4*n - 1))                                       # small-sample correction
hedges_g <- J * cohen_d

# (Optional) quick descriptor
interpret_d <- function(d) {
  d_abs <- abs(d)
  if (d_abs < 0.20) "trivial"
  else if (d_abs < 0.50) "small"
  else if (d_abs < 0.80) "medium"
  else "large"
}
d_label <- interpret_d(cohen_d)

# -----4A. APA-style Reporting (Parametric Results)-------------
if (analysis_method %in% c("t_standard", "t_transformed")) {
  
  sig_label <- ifelse(t_result$p.value < alpha_level, "significantly different", "not significantly different")
  direction_label <- switch(
    tail_type,
    "two.sided" = "different from",
    "greater"   = "greater than",
    "less"      = "less than"
  )
  
  cat(sprintf(
    "A one-sample t-test was conducted to compare %s to the reference value of %.2f (α = %.2f).\n",
    metric_label, mu_reference, alpha_level
  ))
  
  if (analysis_method == "t_transformed") {
    cat(sprintf("The data were %s to meet normality assumptions prior to analysis.\n", label))
  }
  
  cat(sprintf(
    "The mean %s (M = %.2f, SD = %.2f) was %s %.2f, t(%d) = %.2f, p = %.3f, 95%% CI [%.2f, %.2f].\n",
    metric_label,
    mean_val, sd_val,
    sig_label, mu_reference,
    round(t_result$parameter, 0),
    round(t_result$statistic, 2),
    round(t_result$p.value, 3),
    round(t_result$conf.int[1], 2), round(t_result$conf.int[2], 2)
  ))
  
  # Effect size line
  if (analysis_method == "t_transformed") {
    cat(sprintf("Effect size (on %s scale): Cohen's d = %.2f (Hedges' g = %.2f), %s.\n",
                label, cohen_d, hedges_g, d_label))
  } else {
    cat(sprintf("Effect size: Cohen's d = %.2f (Hedges' g = %.2f), %s.\n",
                cohen_d, hedges_g, d_label))
  }
}


# -----4B. APA-style Reporting (Wilcoxon Nonparametric Results)-------------
if (analysis_method == "nonparametric") {
  
  sig_label <- ifelse(wilcox_result$p.value < alpha_level, "significantly", "not significantly")
  direction_label <- switch(
    tail_type,
    "two.sided" = "different from",
    "greater"   = "greater than",
    "less"      = "less than"
  )
  
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
  
  cat(sprintf(
    "A Wilcoxon signed-rank test was conducted to assess whether %s was %s the reference value of %.2f (α = %.2f).\n",
    metric_label, direction_label, mu_reference, alpha_level
  ))
  cat("This nonparametric test was used because the normality assumption could not be satisfied.\n")
  cat(sprintf(
    "The result showed that the median %s was %s %.2f, V = %.2f, p = %.3f.\n",
    metric_label, sig_label, mu_reference, wilcox_result$statistic, wilcox_result$p.value
  ))
}
# End of Script
