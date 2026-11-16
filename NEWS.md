# sqc 0.13.1

## Course renumbered by topic

* The course is now organized in twelve **topics** instead of weeks (1 rational subgrouping, 2 operating characteristics of control charts, 3 operating characteristics of sampling plans, 4 Phase I estimation, 5 process capability, 6 measurement systems, 7 sequential monitoring, 8 count data, 9 dependent data, 10 multivariate monitoring, 11 DOE principles, 12 factorial designs). Every reference to a course week in the documentation of the datasets and functions, `data-raw/make-data.R`, the README, and the vignette now names the topic instead (old Weeks 5 and 6 are Topic 6, Week 7 is Topic 5, Weeks 8 and 9 are Topic 7, and Weeks 10 to 14 are Topics 8 to 12).
* `assay_study` is no longer used by any topic; its documentation says that it is kept for reference and for the examples of `gage_rr()` and `vc_oneway()`. No current topic calls `vc_oneway()` or `sim_vc()`; both are kept for reference.
* The vignette's overview lists the functions and datasets of each topic, and the README lists the datasets by topic.
* No data values, seeds, or function behavior changed.

## Documentation fixes

* `sim_gage()` example: the true %GRR of the example design is about 18, not 17.
* `capability()`: `expected_out` (and `dpmo`) use the SD used for the indices (`sigma` when it is supplied), not always the sample SD.
* `funnel_plot()`: the additive overdispersion limits for proportions are now written correctly as $\theta_0 \pm z_L\sqrt{\theta_0(1 - \theta_0)/s_i + \tau^2}$ (the formula used $s_i$, the size, as if it were a standard error).
* `sim_dependent()`: the residual chart uses the Yule-Walker (lag-1 autocorrelation) estimate of phi, not `ar1_chart()`'s `arima()` fit; the `rep` column and the `n - 1` denominator of `share` for the residual chart are documented.
* `cusum_poisson()` and `cusum_risk()`: the design value is returned as `design` (a named number), and the `observed` and `expected` columns of `cusum_risk()` hold the outcomes and the risks.
* `gage_rr()` documents the returned `alpha`; `sim_box_ellipse()` documents the returned `rho` and `rho_true`.
* Example lines longer than 100 characters (`count_chart()`, `cusum_risk()`, `fishbone()`, `funnel_plot()`, `pareto_chart()`, `sim_gage()`) are wrapped, which clears the `R CMD check` note on Rd line widths.

# sqc 0.13.0

## Designed experiments

* `doe_plan()` makes a randomized plan: completely randomized, or a randomized complete block design with the order of the treatments randomized within each block (reproducible with a `seed` that does not disturb the random number generator).
* `sim_blocking()` simulates the power of the **paired** and the **two-sample** t-test as the unit-to-unit variation grows (the shoe experiment); `sim_randomization()` shows the **bias** and the false alarm probability when a systematic effect (the left foot) is not randomized; `sim_pseudorep()` shows that duplicate measurements of the same units are **not replicates**.
* `factorial_design()` builds the run sheet of a two-level factorial (replicates, random run order, optional level labels, and fractions from `generators`); `alias_structure()` gives the defining relation, the resolution, and the groups of aliased effects.
* `factorial_effects()` estimates main effects and interactions by contrasts (standard errors, t-tests, and confidence intervals from the pure error of a replicated experiment); `predict()` gives predictions; `autoplot()` draws the effects with their intervals or the half-normal plot; `interaction_plot()` draws the interaction plot.
* `doe_power()` and `doe_n()` give the power of an effect test in a two-level factorial and the number of replicates needed.
* `sim_ofat()` compares one-factor-at-a-time experimentation with the factorial by simulation: the probability of choosing the best combination when there is an interaction.
* `fishbone()` draws a cause-and-effect diagram.
* New simulated datasets `helicopter` (a replicated 2^3 experiment with a strong interaction) and `helicopter_confirm` (confirmation flights), stand-ins for the data collected in class.

# sqc 0.12.0

## Multivariate monitoring

* `hotelling_chart()` charts several correlated variables with **Hotelling's T-squared** (the standardized distance from the center, which accounts for the standard deviations and the correlations). Known parameters give the chi-square limit; estimated parameters give the beta limit for Phase I points and the F limit for Phase II points (Tracy, Young, and Mason, 1992). The object has the T-squared statistics, the signals, the marginal z-scores (what separate charts see), and the **contributions** of each variable (the extra distance a variable adds given the others). `autoplot()` draws the T-squared chart or, for two variables, the points with the control **ellipse** and the **box** of the separate charts.
* `t2_ellipse()` gives the points of the ellipse of constant T-squared.
* `sim_box_ellipse()` compares the box of two separate charts with the ellipse by simulation: false alarms, points flagged by one and not the other, a change of the correlation, and a shift.
* `sim_hotelling()` simulates the false alarm probability of the chi-square, beta, and F limits when the mean and covariance are estimated from `m` points.
* New simulated dataset `bp_home`: one patient's morning home blood pressure for 90 days (systolic and diastolic, correlation about 0.78), with day 75 (138/70) inside the box but outside the ellipse and day 82 (144/91) the reverse.

# sqc 0.11.0

## Dependent data

* `acf_plot()` draws the autocorrelation function of a series with the limits that an independent series stays inside.
* `sim_ar1()` simulates AR(1) series with a given marginal standard deviation, optionally with a step shift.
* `i_chart_alarm_rate()` gives the false alarm probability of an individuals chart for AR(1) data: the moving range estimates sigma times sqrt(1 - rho), so the limits are too narrow by that factor (3.4% of the points for rho = 0.5 and 34% for rho = 0.9 are beyond the 3-sigma limits).
* `sim_dependent()` measures by simulation the share of points beyond the limits of the I chart (with the moving-range estimate and with the sample standard deviation) and of the AR(1) residual chart.
* `ar1_chart()` fits an AR(1) model (with `arima()`, or known parameters for Phase II) and draws the **residual chart** (special-cause chart) or the series with its one-step-ahead prediction band; `autoplot()`.
* New simulated dataset `fridge_temp`: a vaccine refrigerator logged every 15 minutes for 10 days (AR(1) with phi = 0.9, all readings within 2 to 8 degrees C, with a door-open event on day 8).

# sqc 0.10.0

## Count data

* `count_chart()` computes **p, np, c, and u charts** for counts and proportions, with the limits that vary with the sample size (or exposure) for p and u charts. `limits = "exact"` gives binomial or Poisson probability limits (each side with a false alarm probability of at most 0.00135 for `k = 3`), which matter for small expected counts, and `adjust = "laney"` gives the Laney p' and u' charts for overdispersed data. `autoplot()` draws the chart with stepped limits.
* `funnel_plot()` draws a **funnel plot** of institutions (observed over expected counts, or proportions) against their size, with 95% and 99.8% limits from the Poisson or binomial distribution (mid-p z-scores), and optional **overdispersion** adjustments (multiplicative or additive, Spiegelhalter 2005). `sim_funnel()` simulates how many institutions are flagged for a given between-institution variation, with or without the adjustment.
* `pareto_chart()` draws a Pareto chart (sorted counts and cumulative percentage).
* New simulated datasets `hospital_infections` (120 hospitals with predicted and observed infections, in the style of CMS Care Compare; **simulated, to be replaced with the published figures**), `hand_hygiene` (weekly audits, p chart), and `clabsi_unit` (monthly central-line infections and line-days, u chart).

# sqc 0.9.0

## CUSUM and EWMA charts

* `cusum_chart()` computes the two-sided **tabular CUSUM** of results in time order (reference value `k`, decision interval `h`, optional head start, optional reset after a signal) and the estimated start of the shift for each signal (the observation after the last time the sum was 0). `ewma_chart()` computes the **EWMA** with exact (time-varying) or asymptotic limits. `autoplot()` draws both.
* `cusum_arl()` and `ewma_arl()` compute the **average run length** for any shift by the integral-equation method (Goel and Wu for the CUSUM, Crowder for the EWMA, solved by Gauss-Legendre quadrature); `shewhart_arl()` gives the ARL of a Shewhart chart for comparison. The two-sided CUSUM ARL combines the two one-sided charts, and a head start gives the fast initial response. `find_cusum_h()` and `find_ewma_L()` find the decision interval or limit width that gives a chosen in-control ARL.
* `sim_arl()` simulates run lengths of a Shewhart, CUSUM, or EWMA chart, and `arl()` summarizes them (the same function as for `sim_run_length()`).

## Likelihood-ratio CUSUMs for counts and risk-adjusted outcomes

* `cusum_poisson()` is the observed-minus-expected CUSUM for counts with expected values (a relative risk to detect), and `cusum_risk()` is the risk-adjusted CUSUM for 0/1 outcomes with predicted risks (an odds ratio to detect; Steiner et al.). Each reports the estimated start of the change and the accumulated excess of observed over expected events since then. `autoplot()` draws them.
* `sim_cusum_poisson()` and `sim_cusum_risk()` simulate the statistic for the series to be monitored, and `cusum_limit()` chooses the decision interval that gives a chosen false alarm probability over the length of the series.
* New simulated datasets `gp_deaths` (a Shipman-style series of observed and expected deaths among a general practitioner's patients over 24 years; **simulated, to be replaced with published figures**) and `surgeon_cases` (three surgeons with predicted risks, one of whom's results worsen).

# sqc 0.8.0

## Process capability

* `capability()` computes Cp, Cpl, Cpu, Cpk (and Cpm with a target) of a sample against one or two specification limits, with the exact chi-square interval for Cp and the Bissell normal-approximation or a percentile bootstrap interval for Cpk, the one-sided lower bound for Cpk, the performance indices Pp and Ppk, the observed and normal-theory fraction outside the specification, and the defects per million. A within-subgroup standard deviation can be supplied with `sigma` to separate the capability indices from the performance indices.
* `sim_capability()` simulates the sampling distribution of Cp and Cpk (and the coverage of the interval) for samples of size `n` from a normal process.
* `dpmo()` and `sigma_level()` convert between Cpk (with an optional shift of the mean) and defects per million, and between defects per million and the sigma level (the 1.5 sigma convention: 3.4 per million is six sigma).
* `cpk_oc()` gives the OC curve of a supplier acceptance rule "accept if the sample Cpk is at least k" (normal approximation), and `capability_plan()` finds the smallest sample size and cut-off that meet a producer's and a consumer's risk, in the spirit of `find_plan()` for sampling plans.
* New datasets `tablet_press` (a press in control that makes 2% of its tablets out of specification) and `tablet_supplier` (a supplier's sample of 30 tablets with Cpk 1.40), used in Week 7.

# sqc 0.7.0

## Gage R&R

* `gage_rr()` analyzes a crossed, balanced measurement system study (parts x operators x replicates) with the two-factor random effects model: the ANOVA table (parts and operators tested against the interaction), the ANOVA estimators of the repeatability, operator, part-by-operator, and part-to-part components, and the Gage R&R summary (reproducibility, gage R&R, total variation, percent contribution, percent study variation or "%GRR", percent of tolerance when a specification width is given, and the number of distinct categories). The interaction is dropped when its p-value is above `alpha` (0.25, the AIAG rule) unless `interaction = "yes"` or `"no"` is given. Negative estimates are shown and treated as 0.
* `sim_gage()` simulates crossed studies with given component SDs and returns the estimated components, %GRR, and number of distinct categories of each, to see how precisely a study of a given size measures the gage and how often an operator component comes out negative.

# sqc 0.6.0

## Variance components

* `vc_oneway()` fits the one-way random effects model by the ANOVA (method of moments) estimators, for balanced or unbalanced data: the ANOVA table, the variance components with their SDs and percent of the total, the intraclass correlation with an F-based confidence interval, and a chi-square interval for the SD within groups. A negative estimate of the between-group variance is reported as it is, flagged, and treated as 0 for the SD, percent, and ICC.
* `sim_vc()` simulates the ANOVA estimators of the one-way model for `a` groups of `n` observations, to see their precision and how often the between-group estimate is negative.
* New dataset `assay_study`: a precision study of a point-of-care HbA1c analyzer (10 samples, 3 technicians, 2 replicates), used in Weeks 5 and 6.

# sqc 0.5.0

## Individuals and moving range charts

* `qc_chart()` has two new chart types, `"I"` and `"MR"`, for individual results in time order. The standard deviation is estimated from the moving ranges, `MRbar / d2` with `d2 = 1.128`, so the I chart limits are `xbar +/- 2.66 MRbar` and the MR chart limits are `0` and `3.267 MRbar` (for `k = 3`). `autoplot()` draws both, and `sigma` and `center` can be supplied for a Phase II chart as before.
* `qc_chart()` has a new argument `exclude` for Phase I: the subgroups (or individual results) with a known assignable cause are left out of the estimation of the center line and sigma, but are still plotted, drawn as crosses, and flagged if they are beyond the limits. For individual results a moving range is dropped when either of its two results is excluded. It works for all four chart types.
* The `qc_chart` object has a new element `index` (the observation or subgroup number of each plotted statistic) and `exclude`, and `as.data.frame()` has an `excluded` column. `signals` now hold observation numbers (for the MR chart, which starts at the second result, these are observation numbers, not positions in the moving-range vector).

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
