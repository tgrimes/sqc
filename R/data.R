#' Four-head bottle filler: one bottle from each head per sample
#'
#' Simulated fills from a four-head filler, used throughout Week 1 (rational subgrouping).
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
