# sqc

Companion R package for a Statistical Quality Control course.

- **Constants:** `d2()`, `d3()`, `c4()`, and `qc_constants()` (a table of `d2`, `d3`, `c4`, `A2`, `A3`, `B3`, `B4`, `D3`, `D4`), computed exactly for any subgroup size.
- **Charts:** `qc_chart()` computes an $\bar X$, $R$, individuals (I), or moving range (MR) chart (center line, limits, signals, and Phase I `exclude` of runs with an assignable cause); `autoplot()` draws it with ggplot2.
- **Laboratory QC:** `lj_chart()` for Levey-Jennings charts with the Westgard rules, and `apply_rules()` for the rules themselves.
- **Run lengths:** `sim_run_length()`, `arl()`, and `arl_table()` simulate the average run length of control rules, with normal or skewed data.
- **Acceptance sampling:** `sampling_plan()`, `oc_curve()`, `plan_risks()`, `find_plan()`, `aoq()`, `ati()`, `aoql()`, and `sim_acceptance()`.
- **Process capability:** `capability()` (Cp, Cpk, Cpm, Pp, Ppk with intervals), `sim_capability()`, `cpk_oc()`, `capability_plan()`, `dpmo()`, and `sigma_level()`.
- **Variance components:** `vc_oneway()` (ANOVA estimators and the intraclass correlation of a one-way random effects model) and `sim_vc()`; `gage_rr()` (crossed Gage R&R: repeatability, reproducibility, %GRR, distinct categories) and `sim_gage()`.
- **Small shifts:** `cusum_chart()` and `ewma_chart()` with `autoplot()`; `cusum_arl()`, `ewma_arl()`, `shewhart_arl()`, `find_cusum_h()`, `find_ewma_L()` (average run lengths by integral equations) and `sim_arl()` (by simulation).
- **Counts and risk-adjusted outcomes:** `cusum_poisson()` (observed minus expected) and `cusum_risk()` (risk-adjusted Bernoulli), with `sim_cusum_poisson()`, `sim_cusum_risk()`, and `cusum_limit()` to choose the decision interval.
- **Count data:** `count_chart()` (p, np, c, u charts, normal or exact limits, Laney adjustment), `funnel_plot()` and `sim_funnel()` (funnel plots with overdispersion adjustments), and `pareto_chart()`.
- **Dependent data:** `acf_plot()`, `sim_ar1()`, `sim_dependent()`, `i_chart_alarm_rate()`, and `ar1_chart()` (AR(1) residual charts).
- **Multivariate data:** `hotelling_chart()` (Hotelling's $T^2$ with Phase I and Phase II limits, contributions, and the control ellipse), `t2_ellipse()`, `sim_box_ellipse()`, and `sim_hotelling()`.
- **Designed experiments:** `doe_plan()`, `sim_blocking()`, `sim_randomization()`, `sim_pseudorep()` (principles); `factorial_design()`, `alias_structure()`, `factorial_effects()` (with `predict()` and `autoplot()`), `interaction_plot()`, `doe_power()`, `doe_n()`, `sim_ofat()` (factorials); `fishbone()`; the datasets `helicopter` and `helicopter_confirm`.
- **Data** (simulated), by course topic: `filler` and `filler_rotating` (Topic 1, rational subgrouping); `tablet_press` and `tablet_supplier` (Topic 5, process capability; `tablet_press` also in Topic 6, measurement systems); `gp_deaths` and `surgeon_cases` (Topic 7, sequential monitoring); `hospital_infections`, `hand_hygiene`, and `clabsi_unit` (Topic 8, count data); `fridge_temp` (Topic 9, dependent data); `bp_home` (Topic 10, multivariate monitoring); `helicopter` and `helicopter_confirm` (Topic 12, factorial designs); and `assay_study`, an HbA1c precision study kept for reference and for the examples of `gage_rr()` and `vc_oneway()`.

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
