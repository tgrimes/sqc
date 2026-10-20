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
