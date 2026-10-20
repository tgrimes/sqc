# sqc 0.4.0

## Run lengths with estimated limits

* `sim_run_length()` and `arl_table()` have a new argument `scale`, the standard deviation of the standardized results (default 1). Together with `shift` it describes a chart whose limits were estimated: with true mean and SD mu and sigma and estimated limits mu-hat and sigma-hat, use `shift = (mu - mu-hat) / sigma-hat` and `scale = sigma / sigma-hat`. The help page gives the formulas and an example.

# sqc 0.3.0

## Acceptance sampling

* `sampling_plan()` creates a single sampling plan `(n, c)`, optionally with a lot size `N` (hypergeometric; the binomial approximation otherwise). `oc_curve()` computes the OC curve of one plan or a list of plans, `autoplot()` draws it and can mark the AQL and LTPD, and `plan_risks()` returns the producer's and consumer's risks.
* `find_plan()` finds the smallest plan that meets a producer's and a consumer's risk.
* `aoq()`, `ati()`, and `aoql()` give the average outgoing quality, the average total inspection, and the AOQ limit under rectifying inspection.
* `sim_acceptance()` estimates the OC curve by simulating the inspection of lots.

# sqc 0.2.0

## Rules, Levey-Jennings charts, and run lengths

* `lj_chart()` computes a Levey-Jennings chart for the control results of one or more levels, given the established mean and SD of each control material, and applies the Westgard rules. `autoplot()` draws one panel per level with lines at the mean and at 1, 2, and 3 SD, and shades rejected runs. `print()` and `as.data.frame()` methods are provided.
* `apply_rules()` applies the rules `"1-2s"`, `"1-3s"`, `"2-2s"`, `"R-4s"`, `"4-1s"`, and `"10x"` (or `"multirule"`, the five rejection rules) to a matrix of z-scores with one row per run and one column per control level, and returns the violations. The 4-1s and 10x rules count consecutive runs on each level separately. The help page lists the null probability of each rule.
* `sim_run_length()` simulates the number of runs until a rejection for one or more rules, any shift, any number of levels, and normal or lognormal controls. `arl()` returns the average run length, and `arl_table()` tabulates the ARL of several rules at several shifts. The simulator reproduces the exact ARLs for 1-3s, 1-2s, 2-2s, 4-1s, and 10x.

# sqc 0.1.0

First release.

## Constants

* `d2()`, `d3()`, and `c4()` return the control chart constants for subgroups of size `n`, and `qc_constants()` returns the full table (`d2`, `d3`, `c4`, `A2`, `A3`, `B3`, `B4`, `D3`, `D4`). They are computed exactly by numerical integration, so any subgroup size of at least 2 is allowed, and each size is computed once per R session.
* Each constant has its own help page with its formula, its use in chart limits, and examples.

## Charts

* `qc_chart()` computes an X-bar or R chart (center line, limits, and signals) from a matrix with one row per subgroup, or from a data frame with `value` and `sample` columns. Supplying `center` and `sigma` gives fixed limits for Phase II monitoring.
* A point beyond the limits is flagged as a signal. Patterns within the limits are not checked by default; set `run_length` (a whole number of at least 2) to also flag runs of that many or more points on one side of the center line.
* `qc_chart()` checks its arguments: `k`, `sigma`, `center`, and `run_length` must be valid numbers, subgroups must have equal size (at least 2) with no missing values, and a warning is given when the estimated sigma is 0.
* `autoplot()` draws a `qc_chart` with ggplot2 and returns a ggplot object. `plot()`, `print()`, and `as.data.frame()` methods are also provided.

## Data

* `filler` (mixed-head sampling plan) and `filler_rotating` (rational subgroup plan) are the four-head bottle filler examples used in Week 1.

## Documentation

* The vignette "Getting started with sqc" walks through the constants, charts, Phase I and Phase II, and the example data.
