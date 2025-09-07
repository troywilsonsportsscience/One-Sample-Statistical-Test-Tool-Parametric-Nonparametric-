# One-Sample Statistical Test Tool (Parametric + Nonparametric)
- Automated one-sample hypothesis testing for continuous data. The tool:
- Checks normality (Shapiro–Wilk by default; Anderson–Darling optional).
- Tries safe transforms (log/sqrt/inv with shifts/epsilons) to rescue normality.
- Chooses one-sample t-test (raw or transformed) or Wilcoxon signed-rank.
- Prints APA-style text and a tidy results table.
- Reports effect sizes: Cohen’s d, Hedges’ g (t path) or r = Z/√n (Wilcoxon).
- Returns a structured result object (settings, test, effects, tidy table, plots, diagnostics).

---
## Setup
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/one-sample-test-tool.git
2. Open the R script (one_sample_test.R) in RStudio or your preferred R environment.
3. Run the script. The first block will automatically install and load all required packages.

---

## User-Configuration Settings
```r
alpha_level  <- 0.05                 # significance threshold
tail_type    <- "two.sided"          # "two.sided", "greater", "less"
metric_label <- "jump height"        # label used in reporting

try_transform      <- TRUE           # attempt transforms if non-normal
allowed_transforms <- c("log","sqrt","inv")
method_normality   <- "shapiro"      # "shapiro", "ad", or "none"

show_plots                         <- TRUE
verbose                            <- TRUE
report_cis_effect_size             <- FALSE   # CI for g via MBESS if installed
report_transform_normality_table   <- TRUE    # print p-values for raw + transforms
show_transform_panel               <- TRUE    # facet histograms for all transforms
```
---

## Dependencies
```r
# The script handles all required packages automatically:
# Required
install.packages(c("tidyverse", "broom", "janitor"))

# Optional (enables Anderson–Darling; CI for Hedges' g)
install.packages(c("nortest", "MBESS"))
```

---
## What it produces
- Console output (if verbose=TRUE)
     - Normality p-value and interpretation
     - If transformed: chosen transform + shift/epsilon used
     - APA-style report for t or Wilcoxon
     - Effect sizes
- Plots (if show_plots=TRUE)
     - Raw histogram with density/normal overlay (p_hist_raw)
     - Final-scale histogram (raw or transformed) (p_hist_final)
     - Raw vs Transformed comparison (only if a transform was used) (p_transform_compare)
     - Transform audit panel: raw + each allowed transform (optional) (p_transform_panel)
- Tidy results table (res$tidy)
     - Method, direction, n, mean/SD (on analyzed scale)
     - t(df) & CI (t path) or V (Wilcoxon)
     - Effect sizes (Cohen’s d, Hedges’ g or r)
     - Transform metadata (name, shift, epsilon)

---
## Example Data Sets - For learning and Testing

### Noraml
```r
set.seed(1)
test_data <- tibble(
  subject_id  = 1:40,
  test_metric = rnorm(40, mean = 56, sd = 4)
)
mu_reference <- 55
```

### Transform-rescue
```r
set.seed(2)
test_data <- tibble(
  subject_id  = 1:35,
  test_metric = rgamma(35, shape = 2, scale = 5)  # skewed
)
mu_reference <- 55
```
### Heavy outliers - likely nonparametric - Wilcoxon
```r
set.seed(3)
test_data <- tibble(
  subject_id  = 1:30,
  test_metric = c(rnorm(25, mean = 55, sd = 4), rep(100, 5))
)
mu_reference <- 55
```
---

## Output Example
### Parametric (Classic, transformed):
A one-sample t-test compared jump height to μ = 55 (α = 0.05).
Data were log-transformed to meet normality.
M = 4.01, SD = 0.10; t(29) = 2.52, p = .018, 95% CI [3.97, 4.06].
Effect size: Cohen’s d = 0.46; Hedges’ g = 0.45.

### Nonparametric (Wilcoxon):
Wilcoxon signed-rank vs μ = 55 (α = 0.05): V = 229, p = .952.
Effect size: r = Z/√n.

---

## Effect Sizes
- Cohen’s d: mean difference in standard deviation units.
- Hedges’ g: small-sample corrected d.
- Both are reported for parametric tests; interpret as:
   - < 0.20 = trivial
   - 0.20–0.49 = small
   - 0.50–0.79 = medium
   - ≥ 0.80 = large

---

## License
MIT License. Free to use, modify, and share.

---

## Contributing
Contributions welcome! If you'd like to suggest improvements, fix bugs, feel free to open a pull request or start a discussion.

---

## Contact

Created and maintained by **Troy Wilson** with the assistance of OpenAI's ChatGPT to help structure code, explain statistical logic, and improve documentation clarity. For issues or feature requests, please use the GitHub Issues tab.
