# ===============================================================
# One-Sample Statistical Test Tool (Parametric + Nonparametric)
# ===============================================================
# Author: Troy Wilson
# Purpose: Automated one-sample hypothesis testing (t or Wilcoxon),
# with normality checks, safe transforms, APA output, and effect sizes.
# ===============================================================

# -------------------------------
# Helper: Install and Load Packages
# -------------------------------
required_packages <- c("tidyverse", "broom", "janitor")
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
invisible(lapply(required_packages, install_if_missing))

# Optional (for AD normality; CIs)
optional_pkgs <- c("nortest", "MBESS")
invisible(lapply(optional_pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) try(install.packages(p), silent = TRUE)))

# ===============================================================
# 0) Settings
alpha_level <- 0.05
tail_type   <- "two.sided"        # "two.sided", "greater", "less"
metric_label <- "jump height"

# Controls
try_transform          <- TRUE
allowed_transforms     <- c("log","sqrt","inv")
method_normality       <- "shapiro"  # "shapiro", "ad", "none"
show_plots             <- TRUE
verbose                <- TRUE
report_cis_effect_size <- FALSE      # set TRUE if using MBESS and you want CI for g
report_transform_normality_table <- TRUE   # print p-values for raw + transforms
show_transform_panel <- TRUE               # facet histograms for all transforms


# Reference value (mu)
mu_reference <- 55

# ===============================================================
# 1) Data
# path <- file.choose()  # Select CSV file
# df <- read_csv(path, show_col_types = FALSE) %>% clean_names() # Convert to snake_case for ease of use
# test_data <- df %>%
#   select(subject_id, test_metric)  # or rename columns to match this format


# Example (non-normal; will trigger transform or Wilcoxon)
set.seed(202)
test_data <- tibble(subject_id = 1:30, test_metric = rexp(30, rate = 1/55))

stopifnot(all(c("subject_id","test_metric") %in% names(test_data)))
test_data <- test_data %>% select(subject_id, test_metric) %>% drop_na()

n <- nrow(test_data)
if (n < 3) stop("Need at least 3 observations.")

# ===============================================================
# ---- 2A. Transform normality audit ----
if (report_transform_normality_table) {
  if (verbose) cat(sprintf(
    "\nTransform normality audit using %s (pass if p > %.2f):\n",
    toupper(method_normality), alpha_level
  ))
  
  # candidate transforms consistent with your allowed list
  avail <- list(log = safe_log, sqrt = safe_sqrt, inv = safe_inv)
  sel   <- avail[names(avail) %in% allowed_transforms]
  
  # evaluate raw
  raw_p <- norm_fun(test_data$test_metric)
  
  # evaluate each transform
  tf_info <- purrr::imap(sel, function(fun, nm) {
    out <- fun(test_data$test_metric)   # returns list(v, s, e)
    p   <- norm_fun(out$v)
    list(name = nm, p = p, v = out$v, shift = out$s, eps = out$e)
  })
  
  # build tidy table
  transform_normality_table <- tibble::tibble(
    Transform = c("raw", sapply(tf_info, `[[`, "name")),
    Shift     = c(0,     sapply(tf_info, `[[`, "shift")),
    Eps       = c(0,     sapply(tf_info, `[[`, "eps")),
    p_value   = c(raw_p, sapply(tf_info, `[[`, "p"))
  ) %>%
    mutate(Pass = if (method_normality == "none") NA else p_value > alpha_level) %>%
    arrange(desc(p_value))
  
  print(transform_normality_table)
  
  # Optional: facet histograms for *all* transforms side-by-side
  if (show_plots && show_transform_panel) {
    comp_df <- dplyr::bind_rows(
      tibble::tibble(scale = "Raw", val = test_data$test_metric, p = raw_p),
      dplyr::bind_rows(lapply(tf_info, function(z)
        tibble::tibble(scale = paste0(stringr::str_to_title(z$name), " transformed"),
                       val = z$v, p = z$p)))
    ) %>%
      mutate(facet_label = paste0(scale, " (p = ", sprintf("%.3f", p), ")"))
    
    p_transform_panel <-
      ggplot2::ggplot(comp_df, ggplot2::aes(x = val)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ..density..), bins = 12,
                              fill = "grey90", color = "black") +
      ggplot2::geom_density(linetype = "dotted") +
      ggplot2::facet_wrap(~ facet_label, scales = "free") +
      ggplot2::labs(title = "Raw vs Candidate Transforms: Normality Check",
                    x = "Value (scale-specific)", y = "Density") +
      ggplot2::theme_minimal()
    
    print(p_transform_panel)
  }
  
  # stash in result later
}

# 2) Normality checks
safe_sw <- function(x) { x <- x[is.finite(x)]; if (length(x) < 3) return(NA_real_); shapiro.test(x)$p.value }
safe_ad <- function(x) {
  if (!requireNamespace("nortest", quietly = TRUE)) return(NA_real_)
  x <- x[is.finite(x)]; if (length(x) < 8) return(NA_real_)
  nortest::ad.test(x)$p.value
}
norm_fun <- switch(method_normality,
                   "shapiro" = safe_sw,
                   "ad"      = safe_ad,
                   "none"    = function(x) NA_real_,
                   safe_sw)

sw_p <- norm_fun(test_data$test_metric)
normality_failed <- if (method_normality == "none") FALSE else (!is.na(sw_p) && sw_p < alpha_level)

if (verbose && method_normality != "none") {
  cat(sprintf("Normality: p = %s → %s\n",
              ifelse(is.na(sw_p), "NA", sprintf("%.4f", sw_p)),
              ifelse(normality_failed, "❌ fail", "✅ pass")))
}

if (show_plots) {
  p_hist <-
    ggplot(test_data, aes(x = test_metric)) +
    geom_histogram(aes(y = ..density..), bins = 12, fill = "grey80", color = "black") +
    stat_function(fun = dnorm,
                  args = list(mean = mean(test_data$test_metric), sd = sd(test_data$test_metric)),
                  linetype = "dashed") +
    labs(title = "Histogram of Test Metric", x = metric_label, y = "Density") +
    theme_minimal()
  print(p_hist)
}

# ---- 2A.1) Safe transforms (with consistent mu-transform) ----
analysis_method <- NA_character_  # "t_raw","t_transformed","wilcoxon"
analysis_var <- NULL
mu_test <- mu_reference
transform_meta <- list(name = "none", shift = 0, eps = 0)

safe_log  <- function(x) { s <- -min(x, na.rm=TRUE) + 1e-6; list(v = log(x + s), s = s, e = 0) }
safe_sqrt <- function(x) { s <- -min(x, na.rm=TRUE); if (s < 0) s <- 0; list(v = sqrt(x + s), s = s, e = 0) }
safe_inv  <- function(x) { s <- -min(x, na.rm=TRUE); if (s < 0) s <- 0; e <- 1e-6; list(v = 1/(x + s + e), s = s, e = e) }

if (!normality_failed) {
  analysis_var <- test_data$test_metric
  mu_test <- mu_reference
  analysis_method <- "t_raw"
} else if (try_transform) {
  if (verbose) cat("↪ Trying transforms...\n")
  avail <- list(log = safe_log, sqrt = safe_sqrt, inv = safe_inv)
  sel   <- avail[names(avail) %in% allowed_transforms]
  tf_info <- purrr::imap(sel, function(fun, nm) {
    out <- fun(test_data$test_metric)
    p   <- norm_fun(out$v)
    list(name = nm, p = p, v = out$v, s = out$s, e = out$e)
  })
  # choose best normality (max p above alpha)
  viable <- Filter(function(z) !is.na(z$p) && z$p > alpha_level, tf_info)
  if (length(viable) > 0) {
    best <- viable[[which.max(sapply(viable, function(z) z$p))]]
    analysis_var <- best$v
    # transform mu accordingly
    mu_test <- switch(best$name,
                      log  = log(mu_reference + best$s),
                      sqrt = sqrt(mu_reference + best$s),
                      inv  = 1/(mu_reference + best$s + best$e))
    transform_meta <- list(name = best$name, shift = best$s, eps = best$e)
    analysis_method <- "t_transformed"
    if (verbose) cat(sprintf("✔️ Using %s transform (SW p = %.4f). Shift=%.3g eps=%.3g\n",
                             best$name, best$p, best$s, best$e))
  } else {
    if (verbose) cat("❌ No transform restored normality → Wilcoxon.\n")
    analysis_method <- "wilcoxon"
  }
} else {
  if (verbose) cat("ℹ️ try_transform = FALSE → Wilcoxon.\n")
  analysis_method <- "wilcoxon"
}

# ---- 2A.2 Visualize chosen scale (raw vs transformed) ----
if (show_plots) {
  # Final histogram on the scale actually analyzed (raw or transformed)
  final_label <- if (analysis_method == "t_transformed")
    paste0(transform_meta$name, "-transformed")
  else
    "raw"
  
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
      else
        metric_label,
      y = "Density"
    ) +
    theme_minimal()
  
  print(p_hist_final)
  
  # Side-by-side comparison ONLY if we used a transform
  if (analysis_method == "t_transformed") {
    comp_df <- dplyr::bind_rows(
      tibble(scale = "Raw",                  val = test_data$test_metric),
      tibble(scale = paste0(stringr::str_to_title(transform_meta$name), " transformed"), val = analysis_var)
    )
    
    p_transform_compare <-
      ggplot(comp_df, aes(x = val)) +
      geom_histogram(aes(y = ..density..), bins = 12, fill = "grey90", color = "black") +
      geom_density(linetype = "dotted") +
      facet_wrap(~ scale, scales = "free") +
      labs(
        title = "Raw vs Transformed Distributions",
        x = "Value (scale-specific)", y = "Density"
      ) +
      theme_minimal()
    
    print(p_transform_compare)
  }
}
plots <- list(
  p_hist_raw = if (show_plots) p_hist else NULL,
  p_hist_final = if (show_plots && exists("p_hist_final")) p_hist_final else NULL,
  p_transform_compare = if (show_plots && exists("p_transform_compare")) p_transform_compare else NULL
)

# ===============================================================
# 3) Tests
t_out <- NULL; w_out <- NULL
if (analysis_method %in% c("t_raw","t_transformed")) {
  t_out <- t.test(analysis_var, mu = mu_test, alternative = tail_type)
} else { # Wilcoxon signed-rank vs mu_reference
  diffs <- test_data$test_metric - mu_reference
  any_ties <- any(diffs == 0)
  # exact only makes sense small n and no ties
  exact_ok <- (n <= 50) && !any_ties
  w_out <- wilcox.test(test_data$test_metric, mu = mu_reference,
                       alternative = tail_type, exact = exact_ok, correct = !exact_ok)
}

# ===============================================================
# 4) Descriptives & Effect sizes
tidy_row <- NULL
effects  <- list()

if (analysis_method %in% c("t_raw","t_transformed")) {
  mean_val <- mean(analysis_var); sd_val <- sd(analysis_var)
  tidy_row <- broom::tidy(t_out)
  # one-sample d = (M - mu)/SD ; unbiased g
  cohen_d  <- (mean_val - mu_test) / sd_val
  J        <- 1 - 3/(4*length(analysis_var) - 1)
  hedges_g <- J * cohen_d
  effects$cohen_d  <- round(cohen_d, 3)
  effects$hedges_g <- round(hedges_g, 3)
  if (report_cis_effect_size && requireNamespace("MBESS", quietly = TRUE)) {
    # MBESS::ci.smd for one-sample: use n1=n, n2=Inf approximates one-sample; pragmatic shortcut:
    ci <- try(MBESS::ci.smd(smd = cohen_d, n.1 = length(analysis_var), n.2 = 1e9,
                            conf.level = 0.95, Unbiased = TRUE), silent = TRUE)
    if (!inherits(ci, "try-error")) {
      effects$hedges_g_95CI <- sprintf("[%.3f, %.3f]", ci$Lower.Conf.Limit.smd, ci$Upper.Conf.Limit.smd)
    }
  }
} else {
  # Nonparametric effect size: r = Z / sqrt(n)
  # For one-sample signed-rank, V = sum of positive ranks.
  V <- as.numeric(w_out$statistic)
  # expected value and variance (no ties, no zeros). Use standard formulas:
  EV <- n*(n+1)/4
  VV <- n*(n+1)*(2*n+1)/24
  # continuity correction as used by wilcox.test when not exact:
  use_cc <- isTRUE(attr(w_out, "parameter")$exact == FALSE) || isTRUE(is.nan(w_out$p.value))
  z <- (V - EV - if (use_cc) 0.5*sign(V - EV) else 0) / sqrt(VV)
  r <- z / sqrt(n)
  effects$wilcoxon_V <- V
  effects$effect_r   <- round(as.numeric(r), 3)
}

if (verbose) {
  if (analysis_method %in% c("t_raw","t_transformed")) {
    cat("\n--- APA (t) ---\n")
    cat(sprintf(
      "One-sample t-test (%s): M = %.2f, SD = %.2f vs μ = %.2f; t(%d) = %.2f, p = %.3f, 95%% CI [%.2f, %.2f].\n",
      ifelse(analysis_method=="t_transformed", paste0(transform_meta$name,"-transformed"), "raw"),
      mean(analysis_var), sd(analysis_var), mu_test,
      round(t_out$parameter), t_out$statistic, t_out$p.value,
      t_out$conf.int[1], t_out$conf.int[2]
    ))
    cat(sprintf("Effect size: Cohen's d = %.3f; Hedges' g = %.3f%s\n",
                effects$cohen_d, effects$hedges_g,
                if (!is.null(effects$hedges_g_95CI)) paste0(" (95% CI ", effects$hedges_g_95CI, ")") else ""))
  } else {
    cat("\n--- APA (Wilcoxon) ---\n")
    cat(sprintf("Wilcoxon signed-rank vs μ = %.2f: V = %.0f, p = %.3f.\n", mu_reference, effects$wilcoxon_V, w_out$p.value))
    cat(sprintf("Effect size: r = %.3f (Z/√n).%s\n", effects$effect_r,
                if (any(test_data$test_metric == mu_reference)) " Note: ties present; exact p-value not computed." else ""))
  }
}

# ===============================================================
# 5) Tidy table
tidy_out <-
  if (analysis_method %in% c("t_raw","t_transformed")) {
    tibble::tibble(
      Method = tidy_row$method,
      `Test Type` = ifelse(analysis_method=="t_transformed","t (transformed)","t (raw)"),
      `Direction` = tail_type,
      n = length(analysis_var),
      `Sample Mean` = round(mean(analysis_var), 2),
      `SD` = round(sd(analysis_var), 2),
      `t(df)` = paste0("t(", round(t_out$parameter, 0), ") = ", round(t_out$statistic, 2)),
      `p-value` = round(t_out$p.value, 4),
      `95% CI` = paste0("[", round(t_out$conf.int[1], 2), ", ", round(t_out$conf.int[2], 2), "]"),
      `Transform` = transform_meta$name,
      `Shift` = transform_meta$shift,
      `Eps` = transform_meta$eps,
      `Cohen_d` = effects$cohen_d,
      `Hedges_g` = effects$hedges_g,
      `Hedges_g_95CI` = if (!is.null(effects$hedges_g_95CI)) effects$hedges_g_95CI else NA_character_
    )
  } else {
    tibble::tibble(
      Method = "Wilcoxon signed-rank",
      `Direction` = tail_type,
      n = n,
      `Reference μ` = mu_reference,
      `V` = effects$wilcoxon_V,
      `p-value` = round(w_out$p.value, 4),
      `Effect r` = effects$effect_r
    )
  }
if (verbose) print(tidy_out)

# ===============================================================
# 6) Structured return (invisible)
plots <- list(
  p_hist_raw = if (show_plots && exists("p_hist")) p_hist else NULL,
  p_hist_final = if (show_plots && exists("p_hist_final")) p_hist_final else NULL,
  p_transform_compare = if (show_plots && exists("p_transform_compare")) p_transform_compare else NULL,
  p_transform_panel = if (show_plots && exists("p_transform_panel")) p_transform_panel else NULL
)

diagnostics <- list(
  transform_normality = if (exists("transform_normality_table")) transform_normality_table else NULL
)

result <- list(
  settings = list(alpha = alpha_level, tail = tail_type, metric_label = metric_label,
                  method_normality = method_normality, try_transform = try_transform,
                  allowed_transforms = allowed_transforms),
  analysis_method = analysis_method,
  transform = transform_meta,
  test = if (analysis_method %in% c("t_raw","t_transformed")) t_out else w_out,
  effects = effects,
  tidy = tidy_out,
  plots = plots,
  diagnostics = diagnostics
)
invisible(result)
#====End of Script====
