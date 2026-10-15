## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4.5,
                      eval = requireNamespace("dplyr", quietly = TRUE))

## ----message = FALSE----------------------------------------------------------
library(recruitlag)
library(dplyr)

## -----------------------------------------------------------------------------
lepanthes_monthly
lepanthes_census |> select(period, census, adults, inflorescences, recruits)
lepanthes_hosts

## -----------------------------------------------------------------------------
k_m <- lag_kernels(6, bin = 1, unit = "mo")
k_6 <- lag_kernels(4, bin = 6, unit = "mo")
c_m <- lag_ceiling(lepanthes_monthly, reproduction = "reproductive_adults", K = 6, lag0 = 0, kernels = k_m)
c_6 <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4, lag0 = 1, kernels = k_6)
c_m
c_6

## -----------------------------------------------------------------------------
rain_strip(`monthly 1994` = c_m, `six-monthly 1999-2004` = c_6)

## -----------------------------------------------------------------------------
profiles <- list(
  `one interval only (projection-matrix null)` = c(1, 0, 0, 0, 0),
  `concentrated early (0 to 6 months)` = c(.7, .2, .07, .03, 0),
  `evenly spread (0 to 30 months)` = rep(.2, 5),
  `concentrated late (24 to 30 months)` = c(0, .03, .07, .2, .7),
  `fitted posterior mean` = lepanthes_profile$weight)
rain_profiles(profiles, bin = 6, unit = "months")
rain_expected(c_6, profiles)

## -----------------------------------------------------------------------------
set.seed(20260913)
k4 <- lag_kernels(4, bin = 6, unit = "mo")
lag_ceiling_boot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                 K = 4, kernels = k4, R = 2000)

## -----------------------------------------------------------------------------
lag_ceiling_bayesboot(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                      K = 4, kernels = k4, draws = 4000)

## ----eval = FALSE-------------------------------------------------------------
# lag_ceiling_stan(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
#                  K = 4, kernels = k4)

## -----------------------------------------------------------------------------
set.seed(20260903)
cal <- ceiling_calibration(c_6)
cal

## -----------------------------------------------------------------------------
lepanthes_hosts |>
  filter(!is.na(recruits)) |>
  group_by(period) |>
  summarise(recruits = sum(recruits), hosts_recruiting = sum(recruits > 0))

## -----------------------------------------------------------------------------
set.seed(20260910)
ht <- host_lag_test(lepanthes_hosts, unit = "host", reproduction = "inflorescences",
                    K = 4, nsim = 4000,
                    kernels = lag_kernels(4, bin = 6, unit = "mo",
                                          extra = list(`fitted posterior mean` = lepanthes_profile$weight)),
                    phi = 0.089)
ht

## -----------------------------------------------------------------------------
plot(ht, statistic = "silent")
plot(ht, statistic = "gini")

## -----------------------------------------------------------------------------
ht$verdict
ht$table |>
  filter(phi_source == "moment",
         label %in% c("flat reference", "pure delay 24-30 mo", "uniform 0-30 mo", "fitted posterior mean")) |>
  select(label, phi, gini_expected, p_silent, p_gini, p_cv)

