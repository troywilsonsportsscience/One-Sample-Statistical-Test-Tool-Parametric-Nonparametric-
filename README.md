# One-Sample Statistical Testing Tool (R)
This project provides a reusable R script for conducting **one-sample hypothesis tests** with full automation, assumption checks, visual diagnostics, and APA-style reporting. It is designed for use in **data science, athlete performance analysis, and research applications**.

---

## Features
- **Normality assessment** using the Shapiro–Wilk test with histogram & density overlay  
- **Optional transformations**: log, square root, inverse  
- **Automatic best transform selection** (highest Shapiro p-value above α)  
- **Fallback to Wilcoxon signed-rank test** if normality cannot be achieved or transformations are disabled  
- **Customizable test direction**: two-sided, greater, or less   
- **APA-style reporting** with test statistic, p-value, CI, and significance wording  
- **Effect sizes for t-tests**: Cohen’s *d* and Hedges’ *g* (with qualitative interpretation)  
- **Faceted ggplot diagnostics** for raw and transformed data  
- **Automatic package installation and loading**

---

## Setup
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/one-sample-test-tool.git
2. Open the R script (one_sample_test.R) in RStudio or your preferred R environment.
3. Run the script. The first block will automatically install and load all required packages.

---

## User-Configuration Settings
- alpha_level      <- 0.05 (significance threshold)
- tail_type        <- "two.sided" (options: "two.sided", "greater", "less")
- metric_label     <- "jump height" (friendly name for reporting)
- try_transform      <- TRUE (set FALSE to skip transformations and forces the script to skip transformations and go directly to Wilcoxon if normality fails)
- allowed_transforms <- c("log", "sqrt", "inv") (choose which transforms to allow)
- show_transform_plot <- TRUE (show diagnostic histograms/densities)
- mu_reference <- 55 (population/reference value)

---

## Dependencies
The script handles all required packages automatically:
- tidyverse (ggplot2, dplyr, tidyr, purrr, stringr)
- broom
- car
- rstatix
- tibble
- janitor

---

## WorkFlow
1. Example test data
set.seed(123)
test_data <- tibble::tibble(
  subject_id = 1:30,
  test_metric = rnorm(30, mean = 57, sd = 5)
)
2. Assumption check (Section 2A)  
   - Histogram with density curve
   - Shapiro–Wilk test
3. Transformation (optional) (Section 2B)
   - Applies allowed transforms (log, sqrt, inv)
   - Keeps only valid ones for your data
   - Picks the best transform if it restores normality
4. Statistical test
   - One-sample t-test (raw or transformed)
   - Wilcoxon signed-rank if normality fails
5. Reporting (Section 4)
      - APA-style test summary
      - Tidy results table
      - Effect sizes (d and g) if parametric

---

## Output Example
### Parametric
A one-sample t-test was conducted to compare jump height to the reference value of 55 (α = 0.05).
The data were log-transformed to meet normality assumptions before analysis.
The mean jump height (M = 4.01, SD = 0.10) was significantly different from 55,
t(29) = 2.52, p = 0.018, 95% CI [3.97, 4.06].
Effect size (on log-transformed scale): Cohen's d = 0.46 (Hedges' g = 0.45), small.

### Nonparametric - Wilcoxon signed-rank test
A Wilcoxon signed-rank test was conducted to assess whether jump height was different from 55 (α = 0.05).
This nonparametric test was used because the normality assumption could not be satisfied.
The result showed that the median jump height was not significantly different from 55,
V = 229, p = 0.952.

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

Created and maintained by **Troy Wilson**. For issues or feature requests, please use the GitHub Issues tab.
This project was developed with the assistance of OpenAI's ChatGPT to help structure code, explain statistical logic, and improve documentation clarity.
