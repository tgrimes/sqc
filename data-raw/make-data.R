# Example datasets for the sqc package. Run from the package root:  Rscript data-raw/make-data.R
# Homework datasets are deliberately not included here.

# ---- filler: four-head filler, mixed-head sampling plan (Topic 1, rational subgrouping) ---
# 40 sampled cycles; one bottle per head each. Target 500 mL, within-head SD 1 mL,
# head 3 runs 3 mL (= 3 SD) high.
set.seed(6)
head_mean <- c(500, 500, 503, 500)
filler <- expand.grid(cycle = 1:40, head = 1:4)
filler$fill_ml <- round(rnorm(nrow(filler), head_mean[filler$head], 1), 1)
filler <- filler[order(filler$cycle, filler$head), ]
rownames(filler) <- NULL
attr(filler, "out.attrs") <- NULL   # left over from expand.grid()

# ---- filler_rotating: same process, rational subgroup plan (Topic 1) ---------
# 5 consecutive cycles from head 1, then 5 from head 2, head 3, head 4, then head 1 again
# (100 sampled bottles, 20 subgroups of 5).
set.seed(6)
filler_rotating <- data.frame(cycle = 1:100)
filler_rotating$head <- ((filler_rotating$cycle - 1) %/% 5) %% 4 + 1
filler_rotating$fill_ml <- round(rnorm(100, head_mean[filler_rotating$head], 1), 1)

save(filler, file = "data/filler.rda", compress = "xz", version = 2)
save(filler_rotating, file = "data/filler_rotating.rda", compress = "xz", version = 2)

# ---- assay_study: precision study of a point-of-care HbA1c analyzer ---------
# Not used by the current topics; kept for reference and for the examples of gage_rr() and vc_oneway().
# 10 patient blood samples (true HbA1c chosen to span the clinical range), 3 technicians each measure
# every sample twice. Technician offsets (A +0.15, B 0, C -0.25) and repeatability SD 0.15 (% HbA1c),
# values displayed to 0.1. The seed was chosen so that the one-way SD within samples is about 0.23,
# giving a smallest detectable change of about 0.65 (above the clinically meaningful 0.5).
true_hba1c <- c(5.5, 5.9, 6.3, 6.7, 7.1, 7.5, 7.9, 8.4, 9.0, 9.7)
tech_offset <- c(A = 0.15, B = 0, C = -0.25)
set.seed(371)
assay_study <- expand.grid(replicate = 1:2, technician = c("A", "B", "C"), sample = 1:10, stringsAsFactors = FALSE)
assay_study$hba1c <- round(true_hba1c[assay_study$sample] + tech_offset[assay_study$technician] +
                             rnorm(nrow(assay_study), 0, 0.15), 1)
assay_study <- assay_study[order(assay_study$sample, assay_study$technician, assay_study$replicate), c("sample", "technician", "replicate", "hba1c")]
rownames(assay_study) <- NULL
attr(assay_study, "out.attrs") <- NULL
save(assay_study, file = "data/assay_study.rda", compress = "xz", version = 2)

# ---- tablet_press and tablet_supplier: Topic 5 (process capability; tablet_press also in Topic 6) ------------------------------------------------
# Specification 475 to 525 mg (target 500). tablet_press: a press in statistical control (20 subgroups of 5
# consecutive tablets), true mean 500 and SD 10.75 mg (Cp about 0.78, 2% outside the specification); the seed
# was chosen so that the X-bar and R charts have no signals and exactly 2 of the 100 tablets are outside.
set.seed(28)
w <- matrix(round(rnorm(100, 500, 10.75), 1), 20, 5)
tablet_press <- data.frame(subgroup = rep(1:20, each = 5), tablet = rep(1:5, times = 20), weight_mg = as.vector(t(w)))
# tablet_supplier: 30 tablets from a supplier whose true mean is 500.4 and SD 5.86 mg (true Cpk about 1.4);
# the seed was chosen so that the sample Cpk is 1.40 and the weights look normal.
set.seed(138)
tablet_supplier <- data.frame(tablet = 1:30, weight_mg = round(rnorm(30, 500.4, 5.86), 1))
save(tablet_press, file = "data/tablet_press.rda", compress = "xz", version = 2)
save(tablet_supplier, file = "data/tablet_supplier.rda", compress = "xz", version = 2)

# ---- gp_deaths: Topic 7 (sequential monitoring). SIMULATED stand-in patterned on the Shipman case ---------
# A GP practice with about 3,000 patients. `expected` is the number of deaths expected each year from the age and sex of the
# patients (rising slowly as the list ages). The true relative risk is 1 in years 1 to 7, 1.2 in years 8 to 16, and 1.35 in
# years 17 to 24 (about 180 excess deaths in all). The seed was chosen so that no annual count is beyond 3 standard
# deviations of the Poisson approximation, while the Poisson CUSUM signals in year 12. These are NOT the real data of the case:
# replace them with the published figures if they are available.
yrs <- 1:24
gp_expected <- round(32 + 0.55 * yrs + c(0.6, -0.4, 0.2, 0.9, -0.7, 0.3, -0.2, 0.5, -0.6, 0.1, 0.8, -0.3, 0.4, -0.9, 0.2, 0.6, -0.5, 0.3, -0.1, 0.7, -0.8, 0.2, 0.4, -0.3), 1)
gp_rr <- c(rep(1, 7), rep(1.2, 9), rep(1.35, 8))
set.seed(133)
gp_deaths <- data.frame(year = yrs, expected = gp_expected, observed = rpois(24, gp_rr * gp_expected))

# ---- surgeon_cases: Topic 7 (sequential monitoring: risk-adjusted CUSUM). SIMULATED -----------------------------------------------
# Three surgeons, 400 consecutive cases each, with the predicted risk of death of each patient (given, as from a risk model).
# Surgeon A operates on high-risk patients (mean predicted risk 9%) and B on low-risk patients (3%); both perform exactly as
# predicted (equal skill). Surgeon C has average-risk patients (5%) and performs as predicted for 250 cases, then the
# odds of death double (case 251 onward). The seed was chosen (searching with risks rounded and floored as below) so that the
# risk-adjusted CUSUM (odds ratio 2, h about 4.2) signals only for C, at about case 344 with the estimated start at case 252,
# while the unadjusted CUSUM (p0 = 0.059, h about 4.46) falsely signals for A (case 136) and misses C, with clear margins on both.
set.seed(6422)
rA <- pmax(round(rbeta(400, 1.5, 15), 4), 1e-4); rB <- pmax(round(rbeta(400, 1.5, 48.5), 4), 1e-4); rC <- pmax(round(rbeta(400, 1.5, 26), 4), 1e-4)
yA <- rbinom(400, 1, rA); yB <- rbinom(400, 1, rB)
orC <- c(rep(1, 250), rep(2, 150)); yC <- rbinom(400, 1, orC * rC / (1 - rC + orC * rC))
surgeon_cases <- data.frame(surgeon = rep(c("A", "B", "C"), each = 400), case = rep(1:400, 3), risk = c(rA, rB, rC), death = c(yA, yB, yC))
save(gp_deaths, file = "data/gp_deaths.rda", compress = "xz", version = 2)
save(surgeon_cases, file = "data/surgeon_cases.rda", compress = "xz", version = 2)

# ---- hospital_infections: Topic 8 (count data: funnel plots). SIMULATED, in the style of the CMS Care Compare healthcare-associated infection data
# 120 hospitals with a predicted number of infections (from a risk model: the predicted counts are skewed, most below 3) and an observed number.
# All hospitals have a true standardized infection ratio (SIR) of 1 except four with a true SIR of about 2.1 (H017, H058, H091, H103) and
# two with about 0.4 (H034, H077), all of medium or large size. The seed was chosen so that the hospital with the highest SIR is a small one
# (H028: 3 observed, 0.54 predicted, SIR 5.6) that is inside the 99.8% funnel limits, so that the four truly worse hospitals are all beyond the
# 99.8% limits, and so that most of the top 10 and bottom 10 hospitals by SIR are small. These are NOT the real data: replace them with the
# published predicted and observed counts (columns "Predicted Cases" and "Observed Cases" of Care Compare) if they are available.
set.seed(318)
m <- 120
pred <- pmin(pmax(round(exp(rnorm(m, log(2.0), 0.9)), 2), 0.15), 9)
theta <- rep(1, m)
idx <- c(17, 58, 91, 103, 34, 77)
pred[idx] <- c(15.4, 12.8, 7.9, 6.3, 14.2, 11.5); theta[idx] <- c(2.1, 2.0, 2.2, 2.4, 0.4, 0.45)
obs <- rpois(m, theta * pred)
hospital_infections <- data.frame(hospital = sprintf("H%03d", 1:m), predicted = pred, observed = obs)

# ---- hand_hygiene: Topic 8 (count data: p chart with varying n). SIMULATED
# 26 weeks of hand-hygiene audits on a hospital unit. The number of observed opportunities varies from 38 to 120 a week and the compliance is
# 78% in weeks 1 to 18 and 90% from week 19 (after a campaign). Week 7 has only 38 audits and a low compliance that is inside the limits.
set.seed(991)
n_aud <- sample(40:120, 26, replace = TRUE); n_aud[7] <- 38
hand_hygiene <- data.frame(week = 1:26, audits = n_aud, compliant = rbinom(26, n_aud, c(rep(0.78, 18), rep(0.90, 8))))

# ---- clabsi_unit: Topic 8 (count data: c and u charts). SIMULATED
# 24 months of central-line infections in an intensive care unit, with the line-days (the exposure) of each month. The rate is 1.4 per 1000
# line-days, except that month 17 has 9 infections in 1450 line-days (a real cluster) and month 6 has 4 infections in 600 line-days, which
# the normal-limit u chart flags and the exact limits do not (a count of 4 is not unusual when only 0.8 are expected).
set.seed(1)
ld <- round(runif(24, 650, 1500), -1); ld[6] <- 600; ld[17] <- 1450
inf <- rpois(24, 1.4 * ld / 1000); inf[17] <- 9; inf[6] <- 4
clabsi_unit <- data.frame(month = 1:24, line_days = ld, infections = inf)

save(hospital_infections, file = "data/hospital_infections.rda", compress = "xz", version = 2)
save(hand_hygiene, file = "data/hand_hygiene.rda", compress = "xz", version = 2)
save(clabsi_unit, file = "data/clabsi_unit.rda", compress = "xz", version = 2)

# ---- fridge_temp: Topic 9 (dependent data). SIMULATED
# A vaccine refrigerator logged every 15 minutes for 10 days (960 readings, from Monday 2026-01-05 00:00). The readings are AR(1) with phi = 0.9,
# mean 5.0 and marginal SD 0.5 (innovation SD 0.218), so they stay within the required 2 to 8 degrees C (range 3.4 to 6.6). On day 8 at 14:00
# (reading 729) the door is left open: the innovation is 1.6 degrees (about 7 innovation SDs) and the temperature then decays back slowly.
# The seed was chosen so that the lag-1 autocorrelation is 0.90, the individuals chart (limits from the moving range) has about 35% of the readings
# beyond its limits, the residual chart with 4-sigma limits flags only reading 729 (with 3-sigma limits one other reading), and an individuals
# chart of every 20th reading (every 5 hours) has no signal.
set.seed(2692)
phi <- 0.9; sig <- 0.5; n_f <- 960
e <- rnorm(n_f, 0, sig * sqrt(1 - phi^2)); e[729] <- e[729] + 1.6
xf <- numeric(n_f); xf[1] <- rnorm(1, 0, sig); for (t in 2:n_f) xf[t] <- phi * xf[t - 1] + e[t]
fridge_temp <- data.frame(time = as.POSIXct("2026-01-05 00:00:00", tz = "UTC") + 15 * 60 * (seq_len(n_f) - 1), temp_c = round(5 + xf, 2))
save(fridge_temp, file = "data/fridge_temp.rda", compress = "xz", version = 2)

# ---- bp_home: Topic 10 (multivariate monitoring). SIMULATED
# One patient's morning home blood pressure (systolic and diastolic, mmHg) for 90 days from Monday 2026-03-02. The two readings are positively
# correlated (true correlation 0.78): systolic mean 125, SD 6; diastolic mean 80, SD 4. Days 1 to 60 are the stable baseline (Phase I) and days
# 61 to 90 are monitored (Phase II). Two readings are planted: day 75 is 138/70, which is inside the separate 3-sigma limits of both variables but
# far outside the control ellipse (the systolic is high while the diastolic is low), and day 82 is 144/91, which is beyond the 3-sigma limit of the
# systolic chart alone but inside the ellipse (both are high together, as they often are). The seed was chosen so that the Phase I estimates are
# close to the truth, T-squared has no other signal in the 90 days, and the separate individuals charts (limits from the moving range) have no
# other signal.
set.seed(48)
n_b <- 90; rho_b <- 0.78
z1 <- rnorm(n_b); z2 <- rnorm(n_b)
sys_b <- round(125 + 6 * z1); dia_b <- round(80 + 4 * (rho_b * z1 + sqrt(1 - rho_b^2) * z2))
sys_b[75] <- 138; dia_b[75] <- 70; sys_b[82] <- 144; dia_b[82] <- 91
bp_home <- data.frame(day = 1:n_b, date = as.Date("2026-03-02") + 0:(n_b - 1), systolic = sys_b, diastolic = dia_b)
save(bp_home, file = "data/bp_home.rda", compress = "xz", version = 2)

# ---- helicopter and helicopter_confirm: Topic 12 (factorial designs). SIMULATED STAND-INS for the data of the paper helicopter run in class
# A replicated 2^3 experiment (16 flights in random order) on three factors: wing length (-1 short, +1 long), body length (-1 short, +1 long), and a
# paper clip (-1 none, +1 one clip). The true model for the flight time (seconds) has a grand mean of 2.4, effects wing +0.30, body +0.04, clip -0.10, and a
# strong wing x clip interaction +0.40 (a long wing needs the clip: the effect of the wing is -0.10 without the clip and +0.70 with it), with an error
# standard deviation of 0.15. The seed was chosen by search so that the estimated effects are close to the truth and the other interactions are small.
# helicopter_confirm: 12 confirmation flights at wing +1, body +1, clip +1 (the true mean is 2.72), in time order (seed 36).
grid <- expand.grid(wing = c(-1, 1), body = c(-1, 1), clip = c(-1, 1))
set.seed(983)
hd <- grid[rep(1:8, 2), ]; hd <- hd[sample(16), ]
mu_h <- 2.4 + 0.15 * hd$wing + 0.02 * hd$body - 0.05 * hd$clip + 0.20 * hd$wing * hd$clip
helicopter <- data.frame(run = 1:16, hd, flight_time = round(mu_h + rnorm(16, 0, 0.15), 2), row.names = NULL)
set.seed(36)
helicopter_confirm <- data.frame(flight = 1:12, flight_time = round(rnorm(12, 2.4 + 0.15 + 0.02 - 0.05 + 0.20, 0.15), 2))
save(helicopter, file = "data/helicopter.rda", compress = "xz", version = 2)
save(helicopter_confirm, file = "data/helicopter_confirm.rda", compress = "xz", version = 2)
