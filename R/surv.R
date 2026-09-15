#' Cumulative calibration assessment for time-to-event outcomes
#'
#' Assesses the calibration of a time-to-event risk prediction model at a fixed
#' horizon \code{tau}, on the cumulative domain, using the exact
#' counting-process formulation. Subjects are ordered by predicted risk at
#' \code{tau} and the martingale residuals
#' \eqn{M_i(\tau) = A_i - \Lambda_i(\tilde{T}_i \wedge \tau)} are accumulated in
#' that order. Under the null hypothesis that the cumulative hazard of the model
#' is correct on \eqn{[0, \tau]}, this partial-sum process converges weakly to
#' Brownian motion once time-changed by its predictable variation
#' \eqn{\sum_i \Lambda_i(\tilde{T}_i \wedge \tau)}.
#'
#' Because the compensator is supplied by the model rather than estimated from
#' the data, the clock is observable and exact: it involves no nuisance
#' parameter, no estimate of the censoring distribution and no variance
#' correction. It is also automatically non-decreasing, since every increment
#' \eqn{\Lambda_i \ge 0}.
#'
#' The method is agnostic to how the predictions were produced. Any model that
#' can report predicted survival \eqn{S_i(t)} at two time points per subject
#' supplies everything that is needed, via \eqn{\Lambda_i(t) = -\log S_i(t)}:
#' Cox, accelerated failure time, parametric and flexible-parametric models,
#' random survival forests and discrete-time neural networks all qualify. Note
#' that \code{p} is evaluated at \code{tau} whereas \code{cumhaz} is evaluated
#' at each subject's own truncated follow-up time \code{min(time, tau)}. These
#' differ for anyone not followed to \code{tau}, so both are genuinely required
#' and neither can be recovered from the other.
#'
#' @section Null hypothesis:
#' This test addresses the whole hazard trajectory on \eqn{[0, \tau]}, not only
#' the risk at \eqn{\tau}. A model whose shape is wrong but whose value at
#' \eqn{\tau} is right will be rejected here. This is a stronger null than the
#' one addressed by a test based only on predicted risks at the horizon.
#'
#' @return an object of class cumulcalib that can be printed or plotted
#' @param y a right-censored \code{Surv} object, as returned by
#'   \code{survival::Surv(time, status)}. Follow-up is truncated at \code{tau}
#'   internally. Counting-process, left-truncated, interval-censored and
#'   multi-state (competing risks) objects are not supported.
#' @param p vector of predicted risks at \code{tau} (a risk, not a survival
#'   probability). Used to order the observations, and as the secondary axis
#'   when the result is plotted.
#' @param tau the prediction horizon; a single finite positive value.
#' @param cumhaz vector of predicted cumulative hazards evaluated at
#'   \code{min(time, tau)}, one value per subject. For a model reporting
#'   predicted survival, \code{cumhaz = -log(S_i(min(time, tau)))}; for a model
#'   reporting a predicted risk \code{r} at that same time,
#'   \code{cumhaz = -log(1 - r)}. If \code{NULL} (the default) the marginal
#'   approach would be used, which is not yet implemented.
#' @param method string with either BB (Brownian bridge test, default method), BM (Brownian motion test), BM2p (two-part BM test - experimental), BB1p (one-part BB test with only the 'bridge' component). Multiple methods can be specified. The first one will be the 'main' method (e.g., when submitting the resulting object to plot()). Default is c("BB","BM")
#' @param ordered if TRUE, the data are already ordered based on ascending values of \code{p}. This is to speed up simulations.
#' @param ties how to handle tied values of \code{p}. \code{"group"} (default) averages the location and time increments within each tied group. \code{"random"} randomly reorders tied observations; this is valid, but results will vary across runs unless a seed is set. \code{"ignore"} keeps tied observations in their input order.
#' @param n_sim if >0, indicates a simulation-based test is requested for inference.
#' @param tol relative tolerance for the internal consistency check between
#'   \code{p} and \code{cumhaz} among subjects followed to \code{tau}, for whom
#'   \code{cumhaz} must equal \code{-log(1 - p)}. The default, \code{0.01}, is
#'   loose enough to accommodate predictions that have been rounded for
#'   reporting, while still far tighter than the discrepancy produced by the
#'   error it guards against (supplying the cumulative hazard at the
#'   untruncated follow-up time, which is typically wrong by O(1) or more).
#' @examples
#' # A correctly specified exponential model, so the null holds by construction
#' set.seed(1)
#' n <- 2000
#' tau <- 1
#' x <- rnorm(n)
#' lambda <- exp(-1.2 + 0.8 * x)
#' Tt <- rexp(n, lambda)
#' Cc <- rexp(n, 0.35)
#' time <- pmin(Tt, Cc)
#' status <- as.numeric(Tt <= Cc)
#' p <- 1 - exp(-lambda * tau)        # predicted risk at tau
#' cumhaz <- lambda * pmin(time, tau) # compensator at min(time, tau)
#' res <- cumulcalibSurv(survival::Surv(time, status), p, tau, cumhaz)
#' summary(res)
#' @export
cumulcalibSurv <- function(
  y,
  p,
  tau,
  cumhaz = NULL,
  method = c("BB", "BM"),
  ordered = FALSE,
  ties = c("group", "random", "ignore"),
  n_sim = 0,
  tol = 0.01
) {
  ties <- match.arg(ties)

  if (is.null(cumhaz)) {
    stop(
      "cumulcalibSurv: 'cumhaz' is required. The marginal approach, which needs only predicted risks at tau, is not yet implemented."
    )
  }

  #---- outcome --------------------------------------------------------------
  if (!inherits(y, "Surv")) {
    stop(
      "cumulcalibSurv: 'y' must be a Surv object, e.g. survival::Surv(time, status)."
    )
  }
  y_type <- attr(y, "type")
  if (is.null(y_type) || y_type != "right") {
    stop(sprintf(
      "cumulcalibSurv: only right-censored data are supported, but 'y' has type \"%s\". Counting-process/left-truncated, interval-censored and multi-state (competing risks) outcomes require a different construction.",
      if (is.null(y_type)) "unknown" else y_type
    ))
  }
  time <- as.numeric(y[, 1])
  status <- as.numeric(y[, 2])

  #---- horizon --------------------------------------------------------------
  if (length(tau) != 1 || !is.finite(tau) || tau <= 0) {
    stop("cumulcalibSurv: 'tau' must be a single finite positive value.")
  }

  #---- lengths and basic validity -------------------------------------------
  n <- length(time)
  if (length(p) != n || length(cumhaz) != n) {
    stop(sprintf(
      "cumulcalibSurv: 'p' (length %d) and 'cumhaz' (length %d) must both have the same length as 'y' (length %d).",
      length(p),
      length(cumhaz),
      n
    ))
  }
  if (anyNA(time) || anyNA(status) || anyNA(p) || anyNA(cumhaz)) {
    stop(
      "cumulcalibSurv: missing values are not allowed in 'y', 'p' or 'cumhaz'."
    )
  }
  if (any(p < 0) || any(p >= 1)) {
    stop(
      "cumulcalibSurv: 'p' must be predicted risks in [0, 1). Note that 'p' is a risk, not a survival probability."
    )
  }
  if (any(cumhaz < 0) || any(!is.finite(cumhaz))) {
    stop("cumulcalibSurv: 'cumhaz' must be finite and non-negative.")
  }

  #---- the compensator must be evaluated at min(time, tau) ------------------
  #For anyone followed to tau we have min(time, tau) = tau, so cumhaz must equal
  #-log(1 - p) exactly. This catches the most common input error, namely
  #supplying Lambda at the subject's untruncated follow-up time (for instance
  #predict(fit, type = "expected") from a model fitted without truncating at tau
  #first). It is an identity of the definitions, so it holds whatever produced
  #the predictions, and it does not look at the observed events.
  reached <- which(time >= tau)
  if (length(reached) > 0) {
    implied <- -log(1 - p[reached])
    bad <- which(abs(cumhaz[reached] - implied) > tol * pmax(1, abs(implied)))
    if (length(bad) > 0) {
      i <- reached[bad[1]]
      stop(sprintf(
        "cumulcalibSurv: 'cumhaz' is inconsistent with 'p' for %d of the %d subjects followed to tau. For those subjects min(time, tau) = tau, so cumhaz must equal -log(1 - p). First mismatch at index %d: cumhaz = %g but -log(1 - p) = %g. The usual cause is supplying the cumulative hazard at the untruncated follow-up time; truncate the data at tau before fitting and predicting. If instead the predictions were merely rounded, raise 'tol' (currently %g).",
        length(bad),
        length(reached),
        i,
        cumhaz[i],
        implied[bad[1]],
        tol
      ))
    }
  }

  #---- event indicator at the horizon ---------------------------------------
  A <- as.numeric(status == 1 & time <= tau)
  if (sum(A) == 0) {
    warning(
      "cumulcalibSurv: no events are observed before tau; the process is deterministic given the predictions."
    )
  }

  #---- order by predicted risk ----------------------------------------------
  if (!ordered) {
    o <- order(p)
    p <- p[o]
    A <- A[o]
    cumhaz <- cumhaz[o]
  }

  #---- ties in p ------------------------------------------------------------
  #Both the location increment (A - cumhaz) and the time increment (cumhaz) vary
  #within a tied group, because subjects sharing a predicted risk at tau may
  #have been followed for different lengths of time. Averaging both within the
  #group preserves the group totals and makes the interior of the path
  #independent of the arbitrary input order of tied rows; this is equivalent to
  #linearly interpolating the cumulative path between the (already correct)
  #group boundaries.
  run_lengths <- rle(p)$lengths
  tied_lengths <- run_lengths[run_lengths > 1]
  if (length(tied_lengths) > 0) {
    n_groups <- length(tied_lengths)
    n_affected <- sum(tied_lengths)
    if (ties == "group") {
      grp <- rep(seq_along(run_lengths), times = run_lengths)
      ends <- cumsum(run_lengths)
      group_mean <- function(v) (diff(c(0, cumsum(v)[ends])) / run_lengths)[grp]
      A <- group_mean(A)
      cumhaz <- group_mean(cumhaz)
      message(sprintf(
        "cumulcalibSurv: %d groups of tied predicted risk values (%d of %d observations) detected; averaging the location and time increments within tied groups (ties = \"group\").",
        n_groups,
        n_affected,
        n
      ))
    } else if (ties == "random") {
      o2 <- order(p, stats::runif(n))
      p <- p[o2]
      A <- A[o2]
      cumhaz <- cumhaz[o2]
      warning(sprintf(
        "cumulcalibSurv: %d groups of tied predicted risk values (%d of %d observations) detected; using ties = \"random\" (tied observations randomly reordered). Results are stochastic and may differ across runs unless a random seed is set.",
        n_groups,
        n_affected,
        n
      ))
    } else {
      warning(sprintf(
        "cumulcalibSurv: %d groups of tied predicted risk values (%d of %d observations) detected, but ties = \"ignore\" was requested; results may depend on the arbitrary order of tied observations in the input data and may not be reproducible.",
        n_groups,
        n_affected,
        n
      ))
    }
  }

  #---- the process ----------------------------------------------------------
  #Clock: the predictable variation of the martingale, which is observable here
  #rather than estimated. Non-decreasing by construction, every increment >= 0.
  s2 <- cumsum(cumhaz)
  T_ <- s2[n]
  if (T_ <= 0) {
    stop(
      "cumulcalibSurv: the total predicted cumulative hazard is zero; the process is degenerate."
    )
  }
  if (T_ < 30) {
    warning(
      "Total observed time (the accumulated predictable variation) is less than 30; the data might be too small for reliable inference."
    )
  }

  C <- cumsum(A - cumhaz) / n #Scaled partial sum of martingale residuals
  t <- s2 / T_
  S <- C * n / sqrt(T_)

  out <- inference(t, S, method)
  out$T <- T_
  out$C_n <- C[n] #Mean calibration error: observed minus expected, divided by n
  out$C_star <- max(abs(C))
  out$tau <- tau
  out$data <- cbind(t = t, S = S, X = p, C = C)
  out$approach <- "exact"

  class(out) <- c("cumulcalib", "cumulcalibSurv")
  return(out)
}
