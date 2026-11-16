#' Hotelling's T-squared chart
#'
#' Charts several **correlated** measurements at once. Charting each variable separately defines a *box* of acceptable values; a point can lie inside the
#' box and still be unusual because the variables do not move together as they usually do. Hotelling's \eqn{T^2} measures how far a point is from the
#' center in units that account for the standard deviations **and the correlations**: for a vector \eqn{x} with mean \eqn{\mu} and covariance
#' matrix \eqn{\Sigma},
#' \deqn{T^2 = (x - \mu)^\top \Sigma^{-1} (x - \mu).}
#' For two variables with correlation \eqn{\rho} and standardized values \eqn{z_1, z_2}, \eqn{T^2 = (z_1^2 - 2\rho z_1 z_2 + z_2^2)/(1 - \rho^2)}. Points with
#' the same \eqn{T^2} lie on an *ellipse*.
#'
#' **Limits.** If the mean and covariance are known (`mean` and `cov` given), \eqn{T^2} follows a chi-square distribution with \eqn{p} degrees of
#' freedom for in-control data and the upper limit is `qchisq(1 - alpha, p)`. If they are estimated from the `m` rows of `phase1` (the default is all rows),
#' the limits are (Tracy, Young, and Mason, 1992):
#' * **Phase I** (the points used to estimate): \eqn{(m - 1)^2/m \cdot B_{1-\alpha}(p/2, (m - p - 1)/2)}, a beta quantile;
#' * **Phase II** (new points): \eqn{p(m + 1)(m - 1)/\{m(m - p)\} \cdot F_{1-\alpha}(p, m - p)}, an F quantile, which is larger than the chi-square limit
#'   when `m` is small.
#'
#' **Which variable drove a signal?** `contributions` has, for each point and variable \eqn{j}, the extra distance that variable adds,
#' \eqn{T^2 - T^2_{(-j)}}, where \eqn{T^2_{(-j)}} is the \eqn{T^2} of the other variables. It is the squared standardized difference between \eqn{x_j} and
#' the value predicted from the other variables, so a variable that is unusual **given the others** has a large contribution. `z` has the ordinary
#' standardized values (what separate charts see).
#'
#' @param data A data frame (or matrix) with one row per observation, in time order.
#' @param vars Names (or column numbers) of the variables to chart. Default: all numeric columns (at least two).
#' @param phase1 Rows used to estimate the mean and covariance: numbers or a logical vector. Default: all rows. Ignored if `mean` and `cov` are given.
#' @param mean,cov Known mean vector and covariance matrix (give both or neither).
#' @param alpha False alarm probability of a point. Default 0.0027 (the 3-sigma value).
#' @return An object of class `hotelling_chart`: a list with `statistics` (a tibble with `index`, the variables, `t2`, `phase` (`"I"`, `"II"`, or
#'   `"known"`), `ucl`, and `signal`), `contributions` and `z` (tibbles with one column per variable), `signals` (the indices of the points beyond the
#'   limit), `mean`, `cov`, `vars`, `p`, `m` (the number of estimation rows; `NA` if known), `alpha`, `ucl` (named limits), and `known`.
#' @references Tracy, N. D., Young, J. C., and Mason, R. L. (1992). Multivariate control charts for individual observations. *Journal of Quality
#'   Technology* 24, 88-95.
#' @seealso [t2_ellipse()], [sim_box_ellipse()], [sim_hotelling()], [qc_chart()] for one variable
#' @examples
#' h <- hotelling_chart(bp_home, c("systolic", "diastolic"), phase1 = 1:60)
#' h
#' autoplot(h)
#' autoplot(h, which = "ellipse")
#' @export
hotelling_chart <- function(data, vars = NULL, phase1 = NULL, mean = NULL, cov = NULL, alpha = 0.0027) {
  d <- as.data.frame(data)
  if (is.null(vars)) vars <- names(d)[vapply(d, is.numeric, logical(1))]
  if (is.numeric(vars)) vars <- names(d)[vars]
  if (!is.character(vars) || length(vars) < 2 || anyNA(vars) || !all(vars %in% names(d))) stop("`vars` must name at least two columns of `data`.", call. = FALSE)
  X <- as.matrix(d[vars]); storage.mode(X) <- "double"
  if (!all(is.finite(X))) stop("The variables must be numeric with no missing values.", call. = FALSE)
  check_probability(alpha, "alpha")
  n <- nrow(X); p <- ncol(X)
  known <- !is.null(mean) || !is.null(cov)
  if (known) {
    if (is.null(mean) || is.null(cov)) stop("Give both `mean` and `cov`, or neither.", call. = FALSE)
    mu <- as.numeric(mean); S <- as.matrix(cov)
    if (length(mu) != p || !all(is.finite(mu))) stop("`mean` must have one value for each variable.", call. = FALSE)
    if (!all(dim(S) == p) || !all(is.finite(S)) || max(abs(S - t(S))) > 1e-8 * max(1, max(abs(S)))) stop("`cov` must be a symmetric ", p, " by ", p, " matrix.", call. = FALSE)
    m <- NA_real_; p1 <- integer(0)
  } else {
    p1 <- if (is.null(phase1)) seq_len(n) else if (is.logical(phase1)) which(phase1) else as.integer(phase1)
    if (anyNA(p1) || any(p1 < 1 | p1 > n) || anyDuplicated(p1)) stop("`phase1` must select rows of `data`.", call. = FALSE)
    m <- length(p1)
    if (m < p + 2) stop("At least ", p + 2, " rows are needed in `phase1` to estimate the covariance of ", p, " variables.", call. = FALSE)
    mu <- colMeans(X[p1, , drop = FALSE]); S <- stats::cov(X[p1, , drop = FALSE])
  }
  if (qr(S)$rank < p) stop("The covariance matrix is singular: some variables are linearly dependent.", call. = FALSE)
  names(mu) <- vars; dimnames(S) <- list(vars, vars)
  t2 <- stats::mahalanobis(X, mu, S)
  if (known) {
    ucl <- c(known = stats::qchisq(1 - alpha, p)); phase <- rep("known", n); lim <- rep(ucl[[1]], n)
  } else {
    ucl <- c(I = (m - 1)^2 / m * stats::qbeta(1 - alpha, p / 2, (m - p - 1) / 2),
             II = p * (m + 1) * (m - 1) / (m * (m - p)) * stats::qf(1 - alpha, p, m - p))
    phase <- rep("II", n); phase[p1] <- "I"; lim <- unname(ucl[phase])
  }
  contrib <- vapply(seq_len(p), function(j) t2 - stats::mahalanobis(X[, -j, drop = FALSE], mu[-j], S[-j, -j, drop = FALSE]), numeric(n))
  if (n == 1) contrib <- matrix(contrib, nrow = 1)
  colnames(contrib) <- vars
  z <- sweep(sweep(X, 2, mu), 2, sqrt(diag(S)), "/")
  sig <- t2 > lim
  stat <- tibble::as_tibble(X)
  stat <- tibble::add_column(stat, index = seq_len(n), .before = 1)
  stat$t2 <- t2; stat$phase <- phase; stat$ucl <- lim; stat$signal <- sig
  structure(list(statistics = stat, contributions = tibble::as_tibble(contrib), z = tibble::as_tibble(z), signals = which(sig), mean = mu, cov = S,
                 vars = vars, p = p, m = m, alpha = alpha, ucl = ucl, known = known), class = "hotelling_chart")
}

#' @export
print.hotelling_chart <- function(x, ...) {
  cat("<hotelling_chart> T-squared chart of ", x$p, " variables (", paste(x$vars, collapse = ", "), ")\n", sep = "")
  if (x$known) cat("  known mean and covariance; limit ", format(round(x$ucl[[1]], 2)), " (chi-square, alpha = ", x$alpha, ")\n", sep = "") else
    cat("  estimated from m = ", x$m, " rows; limits ", format(round(x$ucl[["I"]], 2)), " (Phase I) and ", format(round(x$ucl[["II"]], 2)), " (Phase II), alpha = ", x$alpha, "\n", sep = "")
  if (length(x$signals) == 0) cat("  no point beyond the limit\n") else
    cat("  ", length(x$signals), " point(s) beyond the limit: observation(s) ", paste(utils::head(x$signals, 12), collapse = ", "), if (length(x$signals) > 12) ", ...", "\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.hotelling_chart <- function(x, ...) as.data.frame(x$statistics)

#' The ellipse of constant T-squared
#'
#' Points of the ellipse \eqn{(x - \mu)^\top \Sigma^{-1} (x - \mu) = c} for two variables. The ellipse with \eqn{c} equal to the control limit of a
#' [hotelling_chart()] is the **control region**: a point outside it signals.
#'
#' @param mean Center, a numeric vector of length 2.
#' @param cov A 2 by 2 covariance matrix.
#' @param c The value of \eqn{T^2} on the ellipse. Default: the 99.73% point of the chi-square distribution with 2 degrees of freedom.
#' @param n Number of points on the curve. Default 200.
#' @return A tibble with columns `x` and `y` (the closed curve).
#' @seealso [hotelling_chart()]
#' @examples
#' S <- matrix(c(36, 19.2, 19.2, 16), 2)
#' head(t2_ellipse(c(125, 80), S))
#' @export
t2_ellipse <- function(mean, cov, c = stats::qchisq(0.9973, 2), n = 200) {
  if (!(is.numeric(mean) && length(mean) == 2 && all(is.finite(mean)))) stop("`mean` must be a numeric vector of length 2.", call. = FALSE)
  cov <- as.matrix(cov)
  if (!all(dim(cov) == 2) || qr(cov)$rank < 2) stop("`cov` must be a nonsingular 2 by 2 matrix.", call. = FALSE)
  check_number(c, "c", positive = TRUE); check_count(n, "n", minimum = 8)
  th <- seq(0, 2 * pi, length.out = n)
  pts <- sqrt(c) * t(chol(cov)) %*% rbind(cos(th), sin(th))
  tibble::tibble(x = mean[1] + pts[1, ], y = mean[2] + pts[2, ])
}

#' Draw a T-squared chart with ggplot2
#'
#' `which = "t2"` draws the \eqn{T^2} statistic in time order with its limits (one for Phase I and one for Phase II when the parameters are
#' estimated). `which = "ellipse"` (two variables) draws the points in the plane of the variables with the control **ellipse** and, if `box = TRUE`,
#' the **box** of the two separate \eqn{\pm k}-sigma charts, so that points inside the box but outside the ellipse (and the reverse) can be seen.
#'
#' @param object A [hotelling_chart()] object.
#' @param which `"t2"` (default) or `"ellipse"`.
#' @param box For the ellipse plot, whether to draw the box of the separate charts. Default `TRUE`.
#' @param k Number of standard deviations of the box. Default 3.
#' @param title,subtitle Plot title and subtitle.
#' @param xlab,ylab Axis labels.
#' @param ... Not used.
#' @return A ggplot object.
#' @export
autoplot.hotelling_chart <- function(object, which = c("t2", "ellipse"), box = TRUE, k = 3, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, ...) {
  which <- match.arg(which)
  d <- object$statistics
  d$status <- ifelse(d$signal, "Beyond limit", "No signal")
  cols <- c(`No signal` = "black", `Beyond limit` = "red")
  if (which == "t2") {
    if (is.null(title)) title <- "Hotelling's T-squared chart"
    if (is.null(xlab)) xlab <- "Observation (time order)"
    if (is.null(ylab)) ylab <- expression(T^2)
    g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$index)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$ucl, group = .data$phase), linetype = "dashed") +
      ggplot2::geom_line(ggplot2::aes(y = .data$t2), color = "grey60") +
      ggplot2::geom_point(ggplot2::aes(y = .data$t2, color = .data$status), size = 1.6) +
      ggplot2::scale_color_manual(values = cols, breaks = "Beyond limit", name = NULL)
  } else {
    if (object$p != 2) stop("The ellipse plot needs exactly two variables.", call. = FALSE)
    v <- object$vars
    if (is.null(title)) title <- "Control ellipse"
    if (is.null(xlab)) xlab <- v[1]
    if (is.null(ylab)) ylab <- v[2]
    cu <- if (object$known) object$ucl[[1]] else object$ucl[["II"]]
    e <- t2_ellipse(object$mean, object$cov, cu)
    g <- ggplot2::ggplot()
    if (box) {
      sd <- sqrt(diag(object$cov)); lo <- object$mean - k * sd; hi <- object$mean + k * sd
      g <- g + ggplot2::annotate("rect", xmin = lo[1], xmax = hi[1], ymin = lo[2], ymax = hi[2], fill = NA, color = "#4c78a8", linetype = "dashed")
    }
    g <- g + ggplot2::geom_path(data = e, ggplot2::aes(x = .data$x, y = .data$y), color = "#e45756", linewidth = 0.9) +
      ggplot2::geom_point(data = d, ggplot2::aes(x = .data[[v[1]]], y = .data[[v[2]]], color = .data$status), size = 1.8) +
      ggplot2::scale_color_manual(values = cols, breaks = "Beyond limit", name = NULL)
  }
  g + ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}

#' Box and ellipse for two correlated variables: a simulation
#'
#' Simulates `n` points from a bivariate normal distribution (marginal means 0 and standard deviations 1, correlation `rho_true`, mean moved by `shift`)
#' and charts them two ways with parameters for correlation `rho` and no shift (the in-control state):
#' * the **box**: two separate charts, a point signals if either variable is beyond \eqn{\pm k} standard deviations;
#' * the **ellipse**: Hotelling's \eqn{T^2} with the chi-square limit `qchisq(1 - alpha, 2)`.
#'
#' With the default `rho_true = rho` and `shift = c(0, 0)` the process is in control and every signal is a false alarm: the box has about **twice** the
#' false alarm probability of one chart (two charts), and as `rho` grows most points beyond the ellipse are **inside the box**. Changing `rho_true`
#' (the variables stop moving together) or `shift` shows what each scheme detects.
#'
#' @param rho Correlation used by the charts, between -1 and 1.
#' @param n Number of points. Default 100000.
#' @param k Number of standard deviations of the box. Default 3.
#' @param alpha False alarm probability of the ellipse. Default 0.0027.
#' @param shift Shift of the two means, in standard deviations. Default none.
#' @param rho_true True correlation of the simulated data. Default `rho`.
#' @param points If `TRUE`, return the simulated points (with `t2`, `box`, `ellipse`, and `category`) instead of the summary. Default `FALSE`.
#' @return By default a one-row tibble with `rho`, `rho_true`, and the shares of points beyond the `box`, beyond the `ellipse`, in `box_only`,
#'   `ellipse_only`, `both`, and `neither`. With `points = TRUE`, a tibble with one row per point (`x`, `y`, `t2`, `box`, `ellipse`, `category`).
#' @seealso [hotelling_chart()], [t2_ellipse()]
#' @examples
#' set.seed(1)
#' sim_box_ellipse(0.8)
#' sim_box_ellipse(0.8, rho_true = 0, shift = c(0, 0))
#' @export
sim_box_ellipse <- function(rho, n = 100000, k = 3, alpha = 0.0027, shift = c(0, 0), rho_true = rho, points = FALSE) {
  if (!(is.numeric(rho) && length(rho) == 1 && is.finite(rho) && abs(rho) < 1)) stop("`rho` must be a single number between -1 and 1.", call. = FALSE)
  if (!(is.numeric(rho_true) && length(rho_true) == 1 && is.finite(rho_true) && abs(rho_true) < 1)) stop("`rho_true` must be a single number between -1 and 1.", call. = FALSE)
  if (!(is.numeric(shift) && length(shift) == 2 && all(is.finite(shift)))) stop("`shift` must be two numbers.", call. = FALSE)
  check_count(n, "n", minimum = 10); check_number(k, "k", positive = TRUE); check_probability(alpha, "alpha")
  z1 <- stats::rnorm(n); z2 <- stats::rnorm(n)
  x <- z1 + shift[1]; y <- rho_true * z1 + sqrt(1 - rho_true^2) * z2 + shift[2]
  t2 <- (x^2 - 2 * rho * x * y + y^2) / (1 - rho^2)
  box <- abs(x) > k | abs(y) > k
  ell <- t2 > stats::qchisq(1 - alpha, 2)
  cat <- ifelse(box & ell, "both", ifelse(box, "box only", ifelse(ell, "ellipse only", "neither")))
  if (points) return(tibble::tibble(x = x, y = y, t2 = t2, box = box, ellipse = ell, category = cat))
  sh <- c(box = mean(box), ellipse = mean(ell), box_only = mean(box & !ell), ellipse_only = mean(ell & !box), both = mean(box & ell), neither = mean(!box & !ell))
  tibble::tibble(rho = rho, rho_true = rho_true, box = sh[["box"]], ellipse = sh[["ellipse"]], box_only = sh[["box_only"]], ellipse_only = sh[["ellipse_only"]],
                 both = sh[["both"]], neither = sh[["neither"]])
}

#' Simulate the false alarm probability of T-squared limits when the covariance is estimated
#'
#' Hotelling's \eqn{T^2} is unchanged by a linear change of the variables, so in-control data can be simulated as `p` independent standard normal
#' variables. Each of `reps` replications draws `m` Phase I points, estimates the mean and covariance, and then
#' * computes \eqn{T^2} of the `m` Phase I points and compares it with the chi-square limit and with the beta limit;
#' * draws `new` Phase II points and compares their \eqn{T^2} with the chi-square limit and with the F limit.
#'
#' The chi-square limit is only right when the parameters are known: for small `m` it gives a Phase II false alarm probability that is **too large** and
#' a Phase I one that is **too small**, while the beta and F limits have the nominal `alpha` (see [hotelling_chart()]).
#'
#' @param m Number of Phase I points (must exceed `p + 1`).
#' @param p Number of variables. Default 2.
#' @param reps Number of replications. Default 2000.
#' @param new Number of Phase II points per replication. Default 100.
#' @param alpha Nominal false alarm probability. Default 0.0027.
#' @return A tibble with `phase`, `limit` (`"chisq"`, `"beta"`, `"F"`), `ucl`, and `rate` (the simulated false alarm probability of a point).
#' @examples
#' set.seed(1)
#' sim_hotelling(m = 30, reps = 500)
#' @export
sim_hotelling <- function(m, p = 2, reps = 2000, new = 100, alpha = 0.0027) {
  check_count(p, "p", minimum = 1); check_count(m, "m", minimum = p + 2); check_count(reps, "reps", minimum = 1); check_count(new, "new", minimum = 1)
  check_probability(alpha, "alpha")
  u_chi <- stats::qchisq(1 - alpha, p)
  u_beta <- (m - 1)^2 / m * stats::qbeta(1 - alpha, p / 2, (m - p - 1) / 2)
  u_f <- p * (m + 1) * (m - 1) / (m * (m - p)) * stats::qf(1 - alpha, p, m - p)
  a <- numeric(4)
  for (r in seq_len(reps)) {
    X <- matrix(stats::rnorm(m * p), m, p)
    mu <- colMeans(X); S <- stats::cov(X)
    t1 <- stats::mahalanobis(X, mu, S)
    t2 <- stats::mahalanobis(matrix(stats::rnorm(new * p), new, p), mu, S)
    a <- a + c(mean(t1 > u_chi), mean(t1 > u_beta), mean(t2 > u_chi), mean(t2 > u_f))
  }
  tibble::tibble(phase = c("I", "I", "II", "II"), limit = c("chisq", "beta", "chisq", "F"), ucl = c(u_chi, u_beta, u_chi, u_f), rate = a / reps)
}
