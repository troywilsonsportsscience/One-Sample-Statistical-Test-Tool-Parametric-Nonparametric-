# ===============================================================
# One-Sample Statistical Test Tool (Parametric + Nonparametric)
# ===============================================================
# Author: Troy Wilson
# Purpose: Automated one-sample hypothesis testing (t or Wilcoxon),
# with normality checks, safe transforms (incl. Yeo–Johnson), APA output, and effect sizes.
# ===============================================================

# -------------------------------
# Helper: Install and Load Packages
# -------------------------------
required_packages <- c("tidyverse", "broom", "janitor", "ggplot2")
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
invisible(lapply(required_packages, install_if_missing))

# Optional (for AD normality; effect-size CIs; Yeo–Johnson)
optional_pkgs <- c("nortest", "MBESS", "bestNormalize")
invisible(lapply(optional_pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) {
  try(install.packages(p), silent = TRUE)
}))

# ===============================================================
# 0) Settings
alpha_level <- 0.05
tail_type   <- "two.sided"        # "two.sided", "greater", "less"
metric_label <- "jump height"

# Controls
try_transform          <- FALSE
allowed_transforms     <- c("log","sqrt","inv","yeojohnson")
method_normality       <- "shapiro"  # "shapiro", "ad", "none"
show_plots             <- TRUE
verbose                <- TRUE
report_cis_effect_size <- FALSE      # TRUE if using MBESS and you want CI for g
report_transform_normality_table <- TRUE   # print p-values for raw + transforms
show_transform_panel <- TRUE               # facet histograms for all transforms

# Reference value (mu)
mu_reference <- 55

# ===============================================================
# ---- 1A) Normality helpers ----
safe_sw <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  shapiro.test(x)$p.value
}
safe_ad <- function(x) {
  if (!requireNamespace("nortest", quietly = TRUE)) return(NA_real_)
  x <- x[is.finite(x)]
  if (length(x) < 8) return(NA_real_)
  nortest::ad.test(x)$p.value
}
norm_fun <- switch(method_normality,
                   "shapiro" = safe_sw,
                   "ad"      = safe_ad,
                   "none"    = function(x) NA_real_,
                   safe_sw)

# ---- 1B) Safe transforms ----
safe_log  <- function(x) {
  s <- -min(x, na.rm = TRUE) + 1e-6
  list(v = log(x + s), s = s, e = 0, info = list(type = "log"))
}
safe_sqrt <- function(x) {
  s <- -min(x, na.rm = TRUE); if (s < 0) s <- 0
  list(v = sqrt(x + s), s = s, e = 0, info = list(type = "sqrt"))
}
safe_inv  <- function(x) {
  s <- -min(x, na.rm = TRUE); if (s < 0) s <- 0
  e <- 1e-6
  list(v = 1/(x + s + e), s = s, e = e, info = list(type = "inv"))
}
safe_yj <- function(x) {
  if (!requireNamespace("bestNormalize", quietly = TRUE)) return(NULL)
  obj <- bestNormalize::yeojohnson(x)
  list(v = obj$x.t, s = 0, e = 0, info = list(type = "yeojohnson", object = obj))
}

# ===============================================================
# 2) Data
library(readr); library(dplyr); library(janitor)

path <- file.choose()
df_raw <- read_csv(path, show_col_types = FALSE) %>% clean_names()

if (ncol(df_raw) < 2) stop("CSV must have at least 2 columns: id, metric")

df <- df_raw %>%
  rename(
    subject_id  = 1,
    test_metric = 2
  ) %>%
  select(subject_id, test_metric) %>%
  drop_na()

test_data <- df
n <- nrow(test_data)
if (n < 3) stop("Need at least 3 observations.")

# For testing: simulate data
# set.seed(202)
# test_data <- tibble(subject_id = 1:30, test_metric = rexp(30, rate = 1/55))

stopifnot(all(c("subject_id","test_metric") %in% names(test_data)))
test_data <- test_data %>% select(subject_id, test_metric) %>% drop_na()

n <- nrow(test_data)
if (n < 3) stop("Need at least 3 observations.")

# ===============================================================
# 2A) Transform normality audit
centered <- test_data$test_metric - mu_reference
sw_p <- norm_fun(centered)
normality_failed <- if (method_normality == "none") FALSE else (!is.na(sw_p) && sw_p < alpha_level)
if (verbose && method_normality != "none") {
  cat(sprintf("Normality (on x - μ): p = %s → %s\n",
              ifelse(is.na(sw_p), "NA", sprintf("%.4f", sw_p)),
              ifelse(normality_failed, "❌ fail", "✅ pass")))
}

# ---- 2B) Raw histogram + Q–Q (optional) ----
if (show_plots) {
  p_hist_raw <-
    ggplot(test_data, aes(x = test_metric)) +
    geom_histogram(aes(y = ..density..), bins = 12, fill = "grey80", color = "black") +
    stat_function(fun = dnorm,
                  args = list(mean = mean(test_data$test_metric),
                              sd   = sd(test_data$test_metric)),
                  linetype = "dashed") +
    labs(title = "Histogram of Test Metric (raw)",
         x = metric_label, y = "Density") +
    theme_minimal()
  print(p_hist_raw)
  
  p_qq_raw <-
    ggplot(tibble(val = test_data$test_metric), aes(sample = val)) +
    stat_qq() + stat_qq_line(linetype = "dashed") +
    labs(title = "Q–Q Plot (raw scale)",
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  print(p_qq_raw)
}

# ---- 2C) Transform comparison panels (histograms + Q–Q) ----
if (show_plots && show_transform_panel) {
  # Build candidate list based on allowed_transforms and available helpers
  cand <- list()
  if ("log" %in% allowed_transforms)        cand$log  <- safe_log
  if ("sqrt" %in% allowed_transforms)       cand$sqrt <- safe_sqrt
  if ("inv" %in% allowed_transforms)        cand$inv  <- safe_inv
  if ("yeojohnson" %in% allowed_transforms) cand$yeojohnson <- safe_yj
  
  # Start with raw; p-values test normality of (value - mu_reference)
  comp_list <- list(
    list(name = "raw",
         v = test_data$test_metric,
         p = suppressWarnings(norm_fun(test_data$test_metric - mu_reference)))
  )
  
  # Add each transform (skip if helper returns NULL, e.g., YJ when bestNormalize isn't installed)
  for (nm in names(cand)) {
    out <- cand[[nm]](test_data$test_metric)
    if (is.null(out)) next
    v <- out$v
    v <- v[is.finite(v)]
    p <- if (length(v) < 3) NA_real_ else suppressWarnings(norm_fun(v - mu_reference))
    comp_list[[length(comp_list) + 1]] <- list(name = nm, v = v, p = p)
  }
  
  # Tidy long data for faceting; label facets with transform name + p-value
  comp_df <- dplyr::bind_rows(lapply(comp_list, function(z) {
    tibble::tibble(
      scale = dplyr::case_when(
        z$name == "raw" ~ "Raw",
        TRUE ~ paste0(stringr::str_to_title(z$name), " transformed")
      ),
      val = z$v,
      p = z$p
    )
  })) %>%
    dplyr::mutate(facet_label = paste0(
      scale, " (p = ", ifelse(is.finite(p), sprintf("%.3f", p), "NA"), ")"
    ))
  
  # Panel: histograms + density overlays
  p_transform_panel <-
    ggplot2::ggplot(comp_df, ggplot2::aes(x = val)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..), bins = 12, fill = "grey90", color = "black") +
    ggplot2::geom_density(linetype = "dotted") +
    ggplot2::facet_wrap(~ facet_label, scales = "free") +
    ggplot2::labs(
      title = "Raw vs Candidate Transforms: Normality Check (centered to μ for test)",
      x = "Value (scale-specific)", y = "Density"
    ) +
    ggplot2::theme_minimal()
  print(p_transform_panel)
  
  # Panel: Q–Q plots
  p_transform_qq <-
    ggplot2::ggplot(comp_df, ggplot2::aes(sample = val)) +
    ggplot2::stat_qq() +
    ggplot2::stat_qq_line(linetype = "dashed") +
    ggplot2::facet_wrap(~ facet_label, scales = "free") +
    ggplot2::labs(
      title = "Q–Q Plots: Raw and Candidate Transforms",
      x = "Theoretical Quantiles", y = "Sample Quantiles"
    ) +
    ggplot2::theme_minimal()
  print(p_transform_qq)
}

# ===============================================================
# 3) Choose analysis scale
analysis_method <- NA_character_  # "t_raw","t_transformed","wilcoxon"
analysis_var <- NULL
mu_test <- mu_reference
transform_meta <- list(name = "none", shift = 0, eps = 0, info = NULL)

if (!normality_failed) {
  analysis_var <- test_data$test_metric
  mu_test <- mu_reference
  analysis_method <- "t_raw"
} else if (try_transform) {
  if (verbose) cat("↪ Trying transforms...\n")
  cand <- list()
  if ("log" %in% allowed_transforms)  cand$log  <- safe_log
  if ("sqrt" %in% allowed_transforms) cand$sqrt <- safe_sqrt
  if ("inv" %in% allowed_transforms)  cand$inv  <- safe_inv
  if ("yeojohnson" %in% allowed_transforms) cand$yeojohnson <- safe_yj
  
  tf_info2 <- list()
  for (nm in names(cand)) {
    out <- cand[[nm]](test_data$test_metric)
    if (is.null(out)) next
    v <- out$v
    v_clean <- v[is.finite(v)]
    p <- if (length(v_clean) < 3) NA_real_ else norm_fun(v - mu_reference)
    tf_info2[[nm]] <- list(name = nm, p = p, v = v, s = out$s, e = out$e, info = out$info)
  }
  viable <- Filter(function(z) !is.na(z$p) && z$p > alpha_level, tf_info2)
  if (length(viable) > 0) {
    best <- viable[[which.max(sapply(viable, function(z) z$p))]]
    analysis_var <- best$v
    mu_test <- switch(best$name,
                      log  = log(mu_reference + best$s),
                      sqrt = sqrt(mu_reference + best$s),
                      inv  = 1/(mu_reference + best$s + best$e),
                      yeojohnson = {
                        if (!is.null(best$info$object)) {
                          as.numeric(stats::predict(best$info$object, newdata = mu_reference))
                        } else stop("Yeo–Johnson object missing; cannot transform μ.")
                      })
    transform_meta <- list(name = best$name, shift = best$s, eps = best$e,
                           info = best$info, norm_p = best$p)   # NEW: record p
    analysis_method <- "t_transformed"
    if (verbose) cat(sprintf("✔️ Using %s transform (normality p = %.4f).\n", best$name, best$p))
  } else {
    if (verbose) cat("❌ No transform restored normality → Wilcoxon.\n")
    analysis_method <- "wilcoxon"
    max_tf_p <- suppressWarnings(max(sapply(tf_info2, function(z) z$p), na.rm = TRUE))
    if (!is.finite(max_tf_p)) max_tf_p <- NA_real_
  }
} else {
  if (verbose) cat("ℹ️ try_transform = FALSE → Wilcoxon.\n")
  analysis_method <- "wilcoxon"
}

if (is.na(analysis_method)) stop("Internal error: analysis_method not set.")
if (analysis_method == "wilcoxon") {
  analysis_var <- test_data$test_metric
}

# ---- 3B) Final histogram + Q–Q on analyzed scale ----
if (show_plots) {
  final_label <- if (analysis_method == "t_transformed")
    paste0(transform_meta$name, "-transformed") else "raw"
  
  p_hist_final <-
    ggplot(tibble(val = analysis_var), aes(x = val)) +
    geom_histogram(aes(y = ..density..), bins = 12, fill = "grey80", color = "black") +
    stat_function(
      fun  = dnorm,
      args = list(mean = mean(analysis_var, na.rm = TRUE),
                  sd   = sd(analysis_var,   na.rm = TRUE)),
      linetype = "dashed"
    ) +
    labs(
      title = paste("Final Histogram (", final_label, " scale)", sep = ""),
      x = if (analysis_method == "t_transformed")
        paste(metric_label, "(", final_label, ")", sep = "")
      else metric_label,
      y = "Density"
    ) +
    theme_minimal()
  print(p_hist_final)
  
  p_qq_final <-
    ggplot(tibble(val = analysis_var), aes(sample = val)) +
    stat_qq() + stat_qq_line(linetype = "dashed") +
    labs(title = paste0("Q–Q Plot (", final_label, " scale)"),
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
  print(p_qq_final)
}

# ===============================================================
# 4) Tests
t_out <- NULL; w_out <- NULL; z <- NA_real_
if (analysis_method %in% c("t_raw","t_transformed")) {
  t_out <- t.test(analysis_var, mu = mu_test, alternative = tail_type)
} else {
  diffs <- test_data$test_metric - mu_reference
  any_ties <- any(diffs == 0)
  exact_ok <- (n <= 50) && !any_ties
  used_correct <- !exact_ok
  w_out <- wilcox.test(test_data$test_metric, mu = mu_reference,
                       alternative = tail_type, exact = exact_ok, correct = used_correct)
  V  <- as.numeric(w_out$statistic)
  EV <- n*(n+1)/4
  VV <- n*(n+1)*(2*n+1)/24
  cc <- if (used_correct) 0.5*sign(V - EV) else 0
  z  <- (V - EV - cc) / sqrt(VV)
}

# ===============================================================
# 5) Descriptives & Effect sizes
effects <- list()
mean_raw <- mean(test_data$test_metric, na.rm = TRUE)
sd_raw   <- sd(test_data$test_metric,   na.rm = TRUE)
med_raw  <- median(test_data$test_metric, na.rm = TRUE)

# --- helper to label effect sizes ---
label_d <- function(d) {
  if (is.na(d)) return("NA")
  if (abs(d) < 0.2) return("negligible")
  if (abs(d) < 0.5) return("small")
  if (abs(d) < 0.8) return("medium")
  return("large")
}
label_r <- function(r) {
  if (is.na(r)) return("NA")
  if (abs(r) < 0.1) return("negligible")
  if (abs(r) < 0.3) return("small")
  if (abs(r) < 0.5) return("medium")
  return("large")
}

if (analysis_method %in% c("t_raw","t_transformed")) {
  mean_val <- mean(analysis_var); sd_val <- sd(analysis_var)
  cohen_d  <- (mean_val - mu_test) / sd_val
  J        <- 1 - 3/(4*length(analysis_var) - 1)
  hedges_g <- J * cohen_d
  effects$cohen_d  <- round(cohen_d, 3)
  effects$hedges_g <- round(hedges_g, 3)
  effects$d_label  <- label_d(cohen_d)
  # optional CI for d
  if (report_cis_effect_size && requireNamespace("MBESS", quietly = TRUE)) {
    ci <- try(MBESS::ci.smd(smd = cohen_d, n.1 = length(analysis_var), n.2 = 1e9,
                            conf.level = 0.95, Unbiased = TRUE), silent = TRUE)
    if (!inherits(ci, "try-error")) {
      effects$hedges_g_95CI <- sprintf("[%.3f, %.3f]", ci$Lower.Conf.Limit.smd, ci$Upper.Conf.Limit.smd)
    }
  }
} else {
  effects$wilcoxon_V <- as.numeric(w_out$statistic)
  effects$effect_r   <- round(as.numeric(z / sqrt(n)), 3)
  effects$r_label    <- label_r(effects$effect_r)
}

# ===============================================================
# 6) APA-style output with WHY + interpretation
if (verbose) {
  mu_symbol <- "\u03bc"
  norm_name <- toupper(method_normality)
  sw_piece  <- if (is.finite(sw_p)) sprintf("%s on (x - %s): p = %.3f",
                                            norm_name, mu_reference, sw_p) else paste(norm_name,"NA")
  
  if (analysis_method == "t_raw") {
    cat("\n--- APA (t, raw) ---\n")
    cat(sprintf("One-sample t-test (raw data): M = %.2f, SD = %.2f vs %s = %.2f, t(%d) = %.2f, p = %.3f, 95%% CI [%.2f, %.2f].\n",
                mean_raw, sd_raw, mu_symbol, mu_test,
                round(t_out$parameter), t_out$statistic, t_out$p.value,
                t_out$conf.int[1], t_out$conf.int[2]))
    cat(sprintf("Effect size: Cohen's d = %.3f (%s), Hedges' g = %.3f%s.\n",
                effects$cohen_d, effects$d_label, effects$hedges_g,
                if (!is.null(effects$hedges_g_95CI)) paste0(" (95% CI ", effects$hedges_g_95CI, ")") else ""))
    cat(sprintf("Assumptions: Normality passed on raw data (%s).\n", sw_piece))
    if (t_out$p.value < alpha_level) {
      cat(sprintf("Interpretation: %s was significantly different from %s = %.2f at α = %.2f.\n",
                  metric_label, mu_symbol, mu_reference, alpha_level))
    } else {
      cat(sprintf("Interpretation: No significant evidence that %s differed from %s = %.2f at α = %.2f.\n",
                  metric_label, mu_symbol, mu_reference, alpha_level))
    }
    
  } else if (analysis_method == "t_transformed") {
    cat("\n--- APA (t, transformed) ---\n")
    cat(sprintf("One-sample t-test (%s transform): Descriptives on raw scale: M = %.2f, SD = %.2f. Test conducted on transformed data vs %s = %.2f, t(%d) = %.2f, p = %.3f.\n",
                transform_meta$name,
                mean_raw, sd_raw, mu_symbol, mu_reference,
                round(t_out$parameter), t_out$statistic, t_out$p.value))
    cat(sprintf("Effect size: Cohen's d = %.3f (%s), Hedges' g = %.3f.\n",
                effects$cohen_d, effects$d_label, effects$hedges_g))
    tfp <- if (!is.null(transform_meta$norm_p)) sprintf("p = %.3f", transform_meta$norm_p) else "NA"
    cat(sprintf("Assumptions: Normality failed on raw (%s). %s transform improved normality (%s). Parametric test run on transformed scale.\n",
                sw_piece, transform_meta$name, tfp))
    if (t_out$p.value < alpha_level) {
      cat(sprintf("Interpretation: %s was significantly different from %s = %.2f at α = %.2f.\n",
                  metric_label, mu_symbol, mu_reference, alpha_level))
    } else {
      cat(sprintf("Interpretation: No significant evidence that %s differed from %s = %.2f at α = %.2f.\n",
                  metric_label, mu_symbol, mu_reference, alpha_level))
    }
    
  } else { # Wilcoxon
    cat("\n--- APA (Wilcoxon) ---\n")
    cat(sprintf("Wilcoxon signed-rank test: Median = %.2f vs %s = %.2f, V = %.0f, p = %.3f.\n",
                med_raw, mu_symbol, mu_reference, effects$wilcoxon_V, w_out$p.value))
    cat(sprintf("Effect size: r = %.3f (%s). Note: r is based on the normal approximation; CIs are not provided by default.\n",
                effects$effect_r, effects$r_label))
    if (exists("max_tf_p") && is.finite(max_tf_p)) {
      cat(sprintf("Assumptions: Normality failed on raw (%s). Transforms did not restore normality (best p = %.3f). Nonparametric test chosen.\n",
                  sw_piece, max_tf_p))
    } else if (!try_transform) {
      cat(sprintf("Assumptions: Normality failed on raw (%s). try_transform = FALSE → nonparametric test chosen.\n", sw_piece))
    } else {
      cat(sprintf("Assumptions: Normality failed on raw (%s). Transforms ineffective → nonparametric test chosen.\n", sw_piece))
    }
    if (w_out$p.value < alpha_level) {
      cat(sprintf("Interpretation: %s was significantly different from %s = %.2f at α = %.2f.\n",
                  metric_label, mu_symbol, mu_reference, alpha_level))
    } else {
      cat(sprintf("Interpretation: No significant evidence that %s differed from %s = %.2f at α = %.2f.\n",
                  metric_label, mu_symbol, mu_reference, alpha_level))
    }
  }
}

# ===============================================================
# 7) Tidy table
tidy_out <-
  if (analysis_method %in% c("t_raw","t_transformed")) {
    tibble::tibble(
      Method       = broom::tidy(t_out)$method,
      `Test Type`  = ifelse(analysis_method == "t_transformed", "t (transformed)", "t (raw)"),
      `Direction`  = tail_type,
      n            = length(analysis_var),
      `Sample Mean`= round(mean(analysis_var), 2),
      `SD`         = round(sd(analysis_var), 2),
      `t(df)`      = paste0("t(", round(t_out$parameter, 0), ") = ", round(t_out$statistic, 2)),
      `p-value`    = round(t_out$p.value, 4),
      `95% CI`     = paste0("[", round(t_out$conf.int[1], 2), ", ", round(t_out$conf.int[2], 2), "]"),
      `Transform`  = transform_meta$name,
      `Shift`      = transform_meta$shift,
      `Eps`        = transform_meta$eps,
      `Cohen_d`    = effects$cohen_d,
      `Hedges_g`   = effects$hedges_g,
      `Hedges_g_95CI` = if (!is.null(effects$hedges_g_95CI)) effects$hedges_g_95CI else NA_character_
    )
  } else {
    tibble::tibble(
      Method        = "Wilcoxon signed-rank",
      `Direction`   = tail_type,
      n             = n,
      `Reference μ` = mu_reference,
      `V`           = effects$wilcoxon_V,
      `p-value`     = round(w_out$p.value, 4),
      `Effect r`    = effects$effect_r
    )
  }

if (verbose) print(tidy_out)
