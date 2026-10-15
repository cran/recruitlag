## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4,
                      eval = requireNamespace("dplyr", quietly = TRUE) &&
                             requireNamespace("tidyr", quietly = TRUE))

## ----message = FALSE----------------------------------------------------------
library(recruitlag)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
regimes

## -----------------------------------------------------------------------------
r <- regimes |> filter(record == "aseasonal_gated")
r

## -----------------------------------------------------------------------------
totals <- r |>
  group_by(period) |>
  summarise(reproduction = sum(reproduction),
            recruits = if (all(is.na(recruits))) NA else sum(recruits))
totals

## -----------------------------------------------------------------------------
totals |>
  tidyr::pivot_longer(c(reproduction, recruits), names_to = "series", values_to = "count") |>
  filter(!is.na(count)) |>
  ggplot(aes(period, count)) +
  geom_col(fill = "grey40", na.rm = TRUE) +
  facet_wrap(~ series, ncol = 1, scales = "free_y") +
  theme_minimal()

## -----------------------------------------------------------------------------
rain_gini(totals$reproduction)
rain_gini(totals$recruits)

## -----------------------------------------------------------------------------
ce <- lag_ceiling(totals, K = 4, kernels = lag_kernels(4))
ce

## -----------------------------------------------------------------------------
rain_strip(`aseasonal, gated` = ce)

## -----------------------------------------------------------------------------
delays <- list(
  `one interval only` = c(1, 0, 0, 0, 0),
  `early`             = c(.7, .2, .07, .03, 0),
  `evenly spread`     = rep(.2, 5),
  `late`              = c(0, .03, .07, .2, .7))
rain_profiles(delays)

## -----------------------------------------------------------------------------
rain_expected(ce, delays)

## -----------------------------------------------------------------------------
set.seed(1)
ht <- host_lag_test(r, K = 4, nsim = 2000, kernels = lag_kernels(4))
ht

## -----------------------------------------------------------------------------
ht$table |> select(label, phi, p_silent, p_gini, p_cv) |> arrange(desc(p_gini)) |> head()

## -----------------------------------------------------------------------------
r2 <- regimes |> filter(record == "aseasonal_delayed")
ce2 <- lag_ceiling(r2, K = 4, kernels = lag_kernels(4))   # a unit-by-period table is summed for you
ce2
set.seed(1)
ht2 <- host_lag_test(r2, K = 4, nsim = 2000, kernels = lag_kernels(4))
ht2

## -----------------------------------------------------------------------------
r3 <- regimes |> filter(record == "seasonal_delayed")
totals3 <- r3 |> group_by(period) |>
  summarise(reproduction = sum(reproduction),
            recruits = if (all(is.na(recruits))) NA else sum(recruits))
totals3$reproduction

## -----------------------------------------------------------------------------
rain_gini(totals3$reproduction)
rain_gini(totals3$recruits)

## -----------------------------------------------------------------------------
ce3 <- lag_ceiling(totals3, K = 4, kernels = lag_kernels(4))
ce3
rain_strip(`aseasonal, gated` = ce, `seasonal, true delay` = ce3)

## -----------------------------------------------------------------------------
r4 <- regimes |> filter(record == "seasonal_gated")
ce4 <- lag_ceiling(r4, K = 4, kernels = lag_kernels(4))
ce4
rain_strip(`seasonal, true delay` = ce3, `seasonal, gated` = ce4)

## -----------------------------------------------------------------------------
set.seed(1)
host_lag_test(r4, K = 4, nsim = 2000, kernels = lag_kernels(4))

