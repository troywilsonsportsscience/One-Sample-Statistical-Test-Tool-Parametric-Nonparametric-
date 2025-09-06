# One-Sample Statistical Testing Tool (R)
This project provides a reusable R script for conducting **one-sample hypothesis tests** with full automation, assumption checks, visual diagnostics, and APA-style reporting. It is designed for use in **data science, athlete performance analysis, and research applications**.

---

## Features
- **Normality assessment** using the Shapiro-Wilk test and histogram
- **Automatic data transformation** (log, square root, inverse)
- **Best-fit transformation selection** based on Shapiro-Wilk p-values
- **Fallback to Wilcoxon signed-rank test** if transformations fail
- **Customizable test direction**: two-sided, greater, or less
- **Auto-calculation of mu_reference** (or user-defined)
- **APA-style reporting** (including tail direction and alpha level)
- **Faceted ggplot diagnostics** for all transformed variants
- **Auto-installs required packages**

---

## Setup
1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/one-sample-test-tool.git
2. Open the R script (one_sample_test.R) in RStudio or your preferred R environment.
3. Run the script. The first block will automatically install and load all required packages.

---

## Dependencies
The script handles all required packages automatically:
- tidyverse
- broom
- car
- rstatix
- janitor

---

## Example Use Case
#Example test data
set.seed(123)
test_data <- tibble::tibble(
  subject_id = 1:30,
  test_metric = rnorm(30, mean = 57, sd = 5)
)

### The script will:
- Check for normality (Section 2A)
- Apply and evaluate transformations if needed (Section 2B)
- Run a one-sample t-test or Wilcoxon signed-rank test (Section 3)
- Report APA-style results (Section 4)
- Display diagnostic plots
  
---

## Output Example
A one-sample t-test was conducted to compare jump height to the reference value of 55 (α = 0.05).
The log-transformed mean jump height (M = 4.01, SD = 0.10) was significantly greater than 55,
t(29) = 2.52, p = 0.018, 95% CI [3.97, 4.06].

---

## Customization
- You can edit these user-configurable parameters in Section 1 of the script:
- mu_reference: Population value to test against (optional; auto-calculated if omitted)
- tail_type: Direction of test ("two.sided", "greater", "less")
- alpha_level: Significance threshold (default = 0.05)
- metric_label: Friendly variable name used in output
You may also replace the simulated data section (1B) with your real dataset via CSV or direct assignment.

---

## License
MIT License. Free to use, modify, and share.

---

## Contributing
Contributions welcome! If you'd like to suggest improvements, fix bugs, or add support for other test types (e.g., two-sample or paired tests), feel free to open a pull request or start a discussion.

---

## Contact

Created and maintained by **Troy Wilson**. For issues or feature requests, please use the GitHub Issues tab.
This project was developed with the assistance of OpenAI's ChatGPT to help structure code, explain statistical logic, and improve documentation clarity.
