## Tests for the exact counting-process arm.
##
## Throughout, the data are generated from an exponential model whose rate is
## also handed to cumulcalibSurv() as the prediction. The null therefore holds
## by construction, with no estimation anywhere, so any failure is a failure of
## the construction rather than of a fitted model.

skip_if_not_installed("survival")

sim_exact <- function(n = 2000, tau = 1, crate = 0.35, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- stats::rnorm(n)
  lambda <- exp(-1.2 + 0.8 * x)
  Tt <- stats::rexp(n, lambda)
  Cc <- if (crate > 0) stats::rexp(n, crate) else rep(Inf, n)
  time <- pmin(Tt, Cc)
  status <- as.numeric(Tt <= Cc)
  list(
    y = survival::Surv(time, status),
    p = 1 - exp(-lambda * tau),
    cumhaz = lambda * pmin(time, tau),
    time = time,
    status = status,
    tau = tau
  )
}

test_that("returns a well-formed cumulcalib object", {
  d <- sim_exact(seed = 1)
  res <- cumulcalibSurv(d$y, d$p, d$tau, d$cumhaz)

  expect_s3_class(res, "cumulcalib")
  expect_s3_class(res, "cumulcalibSurv")
  expect_identical(res$approach, "exact")
  expect_identical(res$tau, d$tau)
  expect_true(all(c("t", "S", "X", "C") %in% colnames(res$data)))
  expect_equal(nrow(res$data), length(d$p))
  expect_true(all(is.finite(res$data)))
  expect_gte(min(res$data[, "t"]), 0)
  expect_equal(max(res$data[, "t"]), 1)
  expect_true(res$pval >= 0 && res$pval <= 1)
})

test_that("the clock is non-decreasing by construction", {
  ## This is the defining advantage of the exact arm: every time increment is a
  ## predicted cumulative hazard and so cannot be negative. No monotonisation.
  for (crate in c(0, 0.35, 1.5)) {
    d <- sim_exact(n = 500, crate = crate, seed = 42)
    res <- cumulcalibSurv(d$y, d$p, d$tau, d$cumhaz)
    expect_true(all(diff(res$data[, "t"]) >= 0))
  }
})

test_that("the terminal value is the observed-minus-expected statistic", {
  d <- sim_exact(seed = 7)
  res <- cumulcalibSurv(d$y, d$p, d$tau, d$cumhaz)
  A <- as.numeric(d$status == 1 & d$time <= d$tau)
  expect_equal(res$C_n, mean(A) - mean(d$cumhaz))
  expect_equal(res$T, sum(d$cumhaz))
})

test_that("observations are ordered by predicted risk", {
  d <- sim_exact(n = 300, seed = 3)
  res <- cumulcalibSurv(d$y, d$p, d$tau, d$cumhaz)
  expect_equal(res$data[, "X"], sort(d$p), ignore_attr = TRUE)
})

test_that("ordered = TRUE reproduces ordered = FALSE on pre-sorted input", {
  d <- sim_exact(n = 300, seed = 5)
  o <- order(d$p)
  pre <- cumulcalibSurv(
    survival::Surv(d$time[o], d$status[o]),
    d$p[o],
    d$tau,
    d$cumhaz[o],
    ordered = TRUE
  )
  post <- cumulcalibSurv(d$y, d$p, d$tau, d$cumhaz)
  expect_equal(pre$data, post$data)
  expect_equal(pre$stat, post$stat)
})

test_that("a correctly specified model is not rejected", {
  ## A loose sanity check at a fixed seed, not a calibration study: under the
  ## null the statistic should be unremarkable.
  d <- sim_exact(n = 4000, seed = 11)
  res <- cumulcalibSurv(d$y, d$p, d$tau, d$cumhaz)
  expect_gt(res$pval, 0.05)
  expect_lt(abs(res$S_n), 3)
})

test_that("miscalibration is detected", {
  ## Inflate every predicted hazard by 30%. Note that p and cumhaz must be
  ## inflated together: a miscalibrated model is still internally consistent
  ## about its own predictions, which is exactly what the p/cumhaz gate checks.
  d <- sim_exact(n = 4000, seed = 11)
  lambda <- d$cumhaz / pmin(d$time, d$tau)
  p_mis <- 1 - exp(-1.3 * lambda * d$tau)
  cumhaz_mis <- 1.3 * lambda * pmin(d$time, d$tau)
  res <- cumulcalibSurv(d$y, p_mis, d$tau, cumhaz_mis)
  expect_lt(res$pval, 0.05)
  expect_lt(res$C_n, 0)
})

test_that("cumhaz is required", {
  d <- sim_exact(n = 100, seed = 1)
  expect_error(cumulcalibSurv(d$y, d$p, d$tau), "not yet implemented")
})

test_that("the outcome must be a right-censored Surv object", {
  d <- sim_exact(n = 100, seed = 1)
  expect_error(
    cumulcalibSurv(d$time, d$p, d$tau, d$cumhaz),
    "must be a Surv object"
  )
  ## Counting-process (left-truncated) input must be refused, not mishandled.
  cp <- survival::Surv(rep(0, 100), d$time, d$status)
  expect_error(
    cumulcalibSurv(cp, d$p, d$tau, d$cumhaz),
    "only right-censored data are supported"
  )
})

test_that("tau, lengths and ranges are validated", {
  d <- sim_exact(n = 100, seed = 1)
  expect_error(cumulcalibSurv(d$y, d$p, c(1, 2), d$cumhaz), "single finite positive")
  expect_error(cumulcalibSurv(d$y, d$p, -1, d$cumhaz), "single finite positive")
  expect_error(cumulcalibSurv(d$y, d$p[-1], d$tau, d$cumhaz), "same length")
  expect_error(
    cumulcalibSurv(d$y, rep(1, 100), d$tau, d$cumhaz),
    "risk, not a survival probability"
  )
  expect_error(
    cumulcalibSurv(d$y, d$p, d$tau, -d$cumhaz),
    "finite and non-negative"
  )
})

test_that("cumhaz supplied at the untruncated follow-up time is caught", {
  ## The most likely user error: Lambda at the subject's own follow-up time
  ## rather than at min(time, tau). For subjects followed to tau these must
  ## agree with -log(1 - p), so the mistake is detectable without using the
  ## observed events.
  d <- sim_exact(n = 500, seed = 2)
  expect_gt(sum(d$time >= d$tau), 0)

  bad <- d$cumhaz
  i <- which(d$time >= d$tau)[1]
  bad[i] <- bad[i] * 1.5
  expect_error(cumulcalibSurv(d$y, d$p, d$tau, bad), "inconsistent with 'p'")

  ## And the realistic version of that mistake, untruncated throughout.
  lambda <- d$cumhaz / pmin(d$time, d$tau)
  expect_error(
    cumulcalibSurv(d$y, d$p, d$tau, lambda * d$time),
    "truncate the data at tau"
  )
})

test_that("tied predicted risks are handled and reported", {
  ## Ties are created by rounding the predicted risk. The compensator is then
  ## rebuilt from the rounded risk so that the two remain internally
  ## consistent; otherwise this would exercise the p/cumhaz gate rather than
  ## the tie logic.
  d <- sim_exact(n = 400, seed = 9)
  p_tied <- round(d$p, 2)
  ch_tied <- (-log(1 - p_tied) / d$tau) * pmin(d$time, d$tau)

  expect_message(
    res <- cumulcalibSurv(d$y, p_tied, d$tau, ch_tied),
    "tied predicted risk values"
  )
  ## Averaging within groups preserves the totals, so the endpoint is unchanged.
  A <- as.numeric(d$status == 1 & d$time <= d$tau)
  expect_equal(res$C_n, mean(A) - mean(ch_tied))
  expect_true(all(diff(res$data[, "t"]) >= 0))

  expect_warning(
    cumulcalibSurv(d$y, p_tied, d$tau, ch_tied, ties = "ignore"),
    "ties = \"ignore\""
  )
})

test_that("the consistency gate tolerates rounding but not the truncation bug", {
  ## Rounding the reported risk perturbs the identity slightly; that must pass
  ## at the default tolerance, while the untruncated-cumhaz mistake, which is
  ## wrong by O(1), must still fail.
  ## (rounding also creates ties, hence the tie message; suppressed here so the
  ## expectation is about the gate alone)
  d <- sim_exact(n = 500, seed = 4)
  expect_no_error(
    suppressMessages(cumulcalibSurv(d$y, round(d$p, 3), d$tau, d$cumhaz))
  )

  lambda <- d$cumhaz / pmin(d$time, d$tau)
  expect_error(
    cumulcalibSurv(d$y, d$p, d$tau, lambda * d$time, tol = 0.5),
    "inconsistent with 'p'"
  )
})
