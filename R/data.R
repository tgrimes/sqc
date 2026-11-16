#' Four-head bottle filler: one bottle from each head per sample
#'
#' Simulated fills from a four-head filler, used throughout Topic 1 (rational subgrouping).
#' Target is 500 mL and the bottle-to-bottle SD within a head is 1 mL. Head 3 runs 3 mL high
#' (an assignable cause); the other heads are centered on 500 mL.
#'
#' The data follow the *mixed-head* sampling plan: every 5th filling cycle, one bottle from
#' each of the four heads is measured, so each sample (subgroup) has size 4. The lecture
#' uses the first 20 samples (`cycle <= 20`).
#'
#' @format A data frame with 160 rows (40 samples by 4 heads) and 3 columns:
#' \describe{
#'   \item{cycle}{Sample number, in time order (1 to 40). Sample `k` is the `k`th sampled filling cycle.}
#'   \item{head}{Filler head (1 to 4).}
#'   \item{fill_ml}{Measured fill volume in mL.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 6).
#' @seealso [filler_rotating]
"filler"

#' Four-head bottle filler: rational subgroups of 5 consecutive bottles from one head
#'
#' The same simulated process as [filler] (target 500 mL, within-head SD 1 mL, head 3 runs 3 mL
#' high), sampled with the *rational subgroup plan*: 5 consecutive cycles from head 1, then
#' 5 from head 2, head 3, head 4, then head 1 again, and so on. There are 100 sampled
#' bottles, forming 20 subgroups of 5.
#'
#' Subgroup `i` is rows `5 * (i - 1) + 1` to `5 * i`, or `(cycle - 1) \%/\% 5 + 1`.
#'
#' @format A data frame with 100 rows and 3 columns:
#' \describe{
#'   \item{cycle}{Sampled bottle number in time order (1 to 100).}
#'   \item{head}{Filler head (1 to 4) that filled the bottle.}
#'   \item{fill_ml}{Measured fill volume in mL.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 6).
#' @seealso [filler]
"filler_rotating"

#' Precision study of a point-of-care HbA1c analyzer
#'
#' Simulated results of a precision study. It was the example of the variance components and Gage R&R lessons
#' before the course was renumbered by topic; no current topic uses it (Topic 6, measurement systems, uses its own
#' gauge studies), and it is kept for reference and for the examples of [gage_rr()] and [vc_oneway()]. Ten patient
#' blood samples were each measured twice by each of three technicians (60 results). The true HbA1c
#' values of the samples span 5.5 to 9.7 percent (SD 1.36); the technicians differ slightly (A reads
#' about 0.15 high and C about 0.25 low relative to B) and the repeatability SD is 0.15. Results
#' are displayed to 0.1 percent. The clinically meaningful change in a patient's HbA1c is 0.5
#' percentage points.
#'
#' @format A data frame with 60 rows and 4 columns:
#' \describe{
#'   \item{sample}{Patient sample (1 to 10).}
#'   \item{technician}{Technician (`"A"`, `"B"`, `"C"`).}
#'   \item{replicate}{Replicate measurement (1 or 2).}
#'   \item{hba1c}{Measured HbA1c, in percent of total hemoglobin.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 371).
#' @seealso [vc_oneway()], [gage_rr()]
"assay_study"

#' Tablet press in statistical control
#'
#' Simulated tablet weights of a press that is in statistical control but makes 2% of its tablets outside the
#' specification of 475 to 525 mg. Used in Topic 5 (process capability), Part A, and again in Topic 6 (measurement
#' systems), which asks how much of the observed variation is due to the gauge. The 100 tablets are 20 subgroups of 5
#' consecutive tablets. The true mean is 500 mg and the standard deviation is 10.75 mg
#' (\eqn{C_p = C_{pk} = 0.78}{Cp = Cpk = 0.78}); the X-bar and R charts show no signals.
#'
#' @format A data frame with 100 rows and 3 columns:
#' \describe{
#'   \item{subgroup}{Subgroup (1 to 20), in time order.}
#'   \item{tablet}{Tablet within the subgroup (1 to 5).}
#'   \item{weight_mg}{Tablet weight in mg.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 28).
#' @seealso [capability()], [tablet_supplier]
"tablet_press"

#' A supplier's capability sample of 30 tablets
#'
#' Simulated tablet weights used in Topic 5 (process capability), Part B: a supplier reports a capability index of
#' 1.4 from 30 tablets, and the contract requires 1.33 against the specification of 475 to 525 mg. The sample Cpk is
#' 1.40. The true mean is 500.4 mg and the standard deviation is 5.86 mg.
#'
#' @format A data frame with 30 rows and 2 columns:
#' \describe{
#'   \item{tablet}{Tablet number (1 to 30).}
#'   \item{weight_mg}{Tablet weight in mg.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 138).
#' @seealso [capability()], [tablet_press]
"tablet_supplier"

#' Deaths in a GP practice: a simulated stand-in patterned on the Shipman case
#'
#' **Simulated data**, used in Topic 7 (sequential monitoring), with the structure of the case of Harold Shipman, a general
#' practitioner whose patients died at higher rates than expected for their age and sex over many years. They are **not** the real
#' data of the case, which should be substituted if the published figures are available.
#'
#' The practice has about 3,000 patients. For each of 24 years, `expected` is the number of deaths expected from the age and sex of
#' the patients (it rises slowly as the list ages) and `observed` is the number that occurred. The true relative risk was 1 in years 1 to 7,
#' 1.2 in years 8 to 16, and 1.35 in years 17 to 24, about 180 excess deaths in all. No annual count is more than 3 standard
#' deviations above its expectation (a Shewhart chart never signals), but a Poisson CUSUM signals in year 12.
#'
#' @format A data frame with 24 rows and 3 columns:
#' \describe{
#'   \item{year}{Year (1 to 24).}
#'   \item{expected}{Expected deaths.}
#'   \item{observed}{Observed deaths.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 133).
#' @seealso [cusum_poisson()], [sim_cusum_poisson()]
"gp_deaths"

#' Three surgeons with predicted risks: simulated
#'
#' **Simulated data**, used in Topic 7 (sequential monitoring: the risk-adjusted CUSUM). Three surgeons each have 400 consecutive cases with the predicted
#' probability of death of each patient (given, as from a risk model). Surgeon A operates on high-risk patients (mean predicted risk
#' 9%) and surgeon B on low-risk patients (3%), and both perform exactly as predicted. Surgeon C has average-risk patients (5%)
#' and performs as predicted for 250 cases, after which the odds of death double.
#'
#' @format A data frame with 1200 rows and 4 columns:
#' \describe{
#'   \item{surgeon}{Surgeon (`"A"`, `"B"`, `"C"`).}
#'   \item{case}{Case number in time order (1 to 400).}
#'   \item{risk}{Predicted probability of death.}
#'   \item{death}{Observed outcome: 1 for death, 0 for survival.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 6422).
#' @seealso [cusum_risk()], [sim_cusum_risk()]
"surgeon_cases"

#' Hospital-acquired infections, observed and predicted: simulated
#'
#' **Simulated data** in the style of the CMS Care Compare healthcare-associated infection measures, used in Topic 8 (count data: funnel plots). For
#' each of 120 hospitals, the number of infections **predicted** by a risk model for its patients and the number **observed**; the
#' standardized infection ratio is `observed / predicted`. Most hospitals have a true ratio of 1; four medium and large hospitals (H017,
#' H058, H091, H103) have a true ratio of about 2.1 and two (H034, H077) of about 0.4. The hospital with the highest ratio is a small one
#' (H028) that is inside the 99.8% funnel limits. The real Care Compare figures can replace these without changing any code.
#'
#' @format A data frame with 120 rows and 3 columns:
#' \describe{
#'   \item{hospital}{Hospital identifier (`"H001"` to `"H120"`).}
#'   \item{predicted}{Predicted (expected) number of infections.}
#'   \item{observed}{Observed number of infections.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 318).
#' @seealso [funnel_plot()], [sim_funnel()]
"hospital_infections"

#' Hand-hygiene audits by week: simulated
#'
#' **Simulated data**, used in Topic 8 (count data: the p chart with varying sample sizes). Weekly audits of hand hygiene on a hospital unit: the number of
#' observed opportunities (which varies from 38 to 120) and the number in which hand hygiene was performed. The compliance is 78% in
#' weeks 1 to 18 and 90% from week 19, after a campaign.
#'
#' @format A data frame with 26 rows and 3 columns:
#' \describe{
#'   \item{week}{Week (1 to 26).}
#'   \item{audits}{Number of observed opportunities.}
#'   \item{compliant}{Number of opportunities with hand hygiene.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 991).
#' @seealso [count_chart()]
"hand_hygiene"

#' Central-line infections in an intensive care unit: simulated
#'
#' **Simulated data**, used in Topic 8 (count data: c and u charts). Central-line associated infections in 24 months, with the line-days of each month
#' (the exposure). The rate is 1.4 per 1000 line-days except in month 17, which has a real cluster (9 infections in 1450 line-days), and
#' month 6, which has 4 infections in only 600 line-days (not unusual for so few line-days, although the normal-limit u chart flags it).
#'
#' @format A data frame with 24 rows and 3 columns:
#' \describe{
#'   \item{month}{Month (1 to 24).}
#'   \item{line_days}{Central-line days in the month.}
#'   \item{infections}{Number of central-line associated infections.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 1).
#' @seealso [count_chart()]
"clabsi_unit"

#' Vaccine refrigerator temperature log: simulated
#'
#' **Simulated data**, used in Topic 9 (dependent data). The temperature of a vaccine refrigerator logged every 15 minutes for 10 days (960
#' readings, starting on Monday 5 January 2026 at 00:00). The readings follow an AR(1) series with parameter 0.9 (consecutive readings are
#' highly correlated), a mean of 5.0 degrees C and a marginal standard deviation of 0.5, so they stay within the 2 to 8 degrees C required for
#' vaccines. On day 8 at 14:00 the door is left open and the temperature rises by 1.6 degrees (a sudden departure that then decays slowly).
#'
#' @format A data frame with 960 rows and 2 columns:
#' \describe{
#'   \item{time}{Time of the reading (POSIXct, UTC).}
#'   \item{temp_c}{Temperature in degrees C.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed given there).
#' @seealso [acf_plot()], [ar1_chart()], [sim_dependent()]
"fridge_temp"

#' Home blood pressure of one patient: simulated
#'
#' **Simulated data**, used in Topic 10 (multivariate monitoring). The morning home blood pressure (systolic and diastolic, in mmHg) of one patient on
#' 90 consecutive days from Monday 2 March 2026. The two readings are positively correlated (correlation about 0.78): when the systolic reading is
#' high the diastolic reading tends to be high too. Days 1 to 60 are the baseline (Phase I) and days 61 to 90 are monitored (Phase II). Day 75
#' (138/70) is inside the separate 3-sigma limits of both variables but far outside the control ellipse, and day 82 (144/91) is beyond the 3-sigma
#' limit of the systolic chart but inside the ellipse. The charts describe how unusual a reading is **for this patient**; they are not clinical
#' thresholds.
#'
#' @format A data frame with 90 rows and 4 columns:
#' \describe{
#'   \item{day}{Day (1 to 90).}
#'   \item{date}{Date of the reading.}
#'   \item{systolic}{Systolic blood pressure (mmHg).}
#'   \item{diastolic}{Diastolic blood pressure (mmHg).}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 48).
#' @seealso [hotelling_chart()], [sim_box_ellipse()]
"bp_home"

#' Paper helicopter factorial experiment: simulated stand-in
#'
#' **Simulated data**, a stand-in for the data that the class collects in Topic 12 (factorial designs; replace it with the class data). A replicated \eqn{2^3} factorial
#' experiment with the paper helicopter: the flight time (seconds) of 16 flights, in the random order in which they were made, at every combination of
#' three factors coded -1 (low) and +1 (high). The true model has a grand mean of 2.4 s, effects of +0.30 (wing), +0.04 (body), -0.10 (clip), and a
#' strong wing by clip interaction of +0.40: a long wing is worse than a short wing without the clip and much better with it, so that changing one factor
#' at a time from the short-wing, no-clip baseline never finds the best combination. The error standard deviation is 0.15 s.
#'
#' @format A data frame with 16 rows and 5 columns:
#' \describe{
#'   \item{run}{Order of the flight (1 to 16).}
#'   \item{wing}{Wing length: -1 short, +1 long.}
#'   \item{body}{Body length: -1 short, +1 long.}
#'   \item{clip}{Paper clip: -1 none, +1 one clip.}
#'   \item{flight_time}{Flight time in seconds.}
#' }
#' @source Simulated; see `data-raw/make-data.R` (seed 983).
#' @seealso [factorial_effects()], [factorial_design()], [interaction_plot()], [helicopter_confirm]
"helicopter"

#' Paper helicopter confirmation flights: simulated stand-in
#'
#' **Simulated data**: 12 confirmation flights, in time order, with the best settings of the [helicopter] experiment (long wing, long body, one clip).
#' The true mean flight time is 2.72 s and the standard deviation 0.15 s. They are charted with an individuals chart to confirm that the process is
#' stable and that the predicted improvement is real.
#'
#' @format A data frame with 12 rows and 2 columns: `flight` (1 to 12) and `flight_time` (seconds).
#' @source Simulated; see `data-raw/make-data.R` (seed 36).
#' @seealso [helicopter]
"helicopter_confirm"
