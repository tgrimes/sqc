# sqc

Companion R package for a Statistical Quality Control course.

- **Constants:** `d2()`, `d3()`, `c4()`, and `qc_constants()` (a table of `d2`, `d3`, `c4`, `A2`, `A3`, `B3`, `B4`, `D3`, `D4`), computed exactly for any subgroup size.
- **Charts:** `qc_chart()` computes an $\bar X$ or $R$ chart (center line, limits, signals); `autoplot()` draws it with ggplot2.
- **Laboratory QC:** `lj_chart()` for Levey-Jennings charts with the Westgard rules, and `apply_rules()` for the rules themselves.
- **Run lengths:** `sim_run_length()`, `arl()`, and `arl_table()` simulate the average run length of control rules, with normal or skewed data.
- **Acceptance sampling:** `sampling_plan()`, `oc_curve()`, `plan_risks()`, `find_plan()`, `aoq()`, `ati()`, `aoql()`, and `sim_acceptance()`.
- **Data:** `filler` and `filler_rotating`, the four-head filler examples.

Dependencies are limited to ggplot2, rlang, and tibble.

## Install

Install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("tgrimes/sqc")
```

To also build the vignette, use `remotes::install_github("tgrimes/sqc", build_vignettes = TRUE)`.

## Example

```r
library(sqc)

qc_chart(filler, type = "xbar", value = fill_ml, sample = cycle) |>
  autoplot(subtitle = "Four-head filler, one bottle from each head")

qc_constants(2:6)
```

See the vignette for a walkthrough: `vignette("sqc")`.

The example data are simulated; `data-raw/make-data.R` regenerates them.
