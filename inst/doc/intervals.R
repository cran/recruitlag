## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 3.6)

## ----message = FALSE----------------------------------------------------------
library(recruitlag)

## ----message = FALSE, warning = FALSE-----------------------------------------
set.seed(20260913)
k <- lag_kernels(4, bin = 6, unit = "mo")
b  <- lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                       K = 4, kernels = k, R = 2000)
bb <- lag_ceiling_bayesboot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                            K = 4, kernels = k, draws = 4000)
knitr::kable(ceiling_draws_table(`cluster bootstrap` = b, `Bayesian bootstrap` = bb)[
  , c("record", "quantity", "estimate", "lower", "upper", "interval", "probability", "probability_is")],
  digits = 2)
# the numbers quoted in the sentences below, taken from the two objects
ex <- function(o) o$table[o$table$quantity == "exceedance Gini", ]
fb <- ex(b); fbb <- ex(bb)

## -----------------------------------------------------------------------------
ex <- bb$draws[, "exceedance_gini"]
eti <- quantile(ex, c(0.05, 0.95))
hdi <- function(v, level = 0.9) { v <- sort(v); n <- length(v); k <- floor(level * n)
  w <- v[(k + 1):n] - v[1:(n - k)]; i <- which.min(w); c(v[i], v[i + k]) }
rbind(`equal-tailed, on the ratio`       = eti,
      `equal-tailed, on log, back`       = exp(quantile(log(ex), c(0.05, 0.95))),
      `highest-density, on the ratio`    = hdi(ex),
      `highest-density, on log, back`    = exp(hdi(log(ex))))

## ----fig.height = 2.6---------------------------------------------------------
rain_forest(`23 host trees` = b, `23 host trees` = bb, index = "both")

