## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4.5,
                      eval = requireNamespace("dplyr", quietly = TRUE) &&
                             requireNamespace("tidyr", quietly = TRUE))

## ----message = FALSE----------------------------------------------------------
library(recruitlag)
library(dplyr)
library(tidyr)
library(ggplot2)

## -----------------------------------------------------------------------------
regimes
regimes_truth

## -----------------------------------------------------------------------------
totals <- regimes |>
  group_by(record, label, period) |>
  summarise(reproduction = sum(reproduction),
            recruits = if (all(is.na(recruits))) NA_real_ else sum(recruits),
            .groups = "drop")
totals

## ----fig.height = 8-----------------------------------------------------------
totals |>
  pivot_longer(c(reproduction, recruits), names_to = "series", values_to = "count") |>
  filter(!is.na(count)) |>                      # unscored periods are not drawn
  mutate(series = factor(series, c("reproduction", "recruits")),
         label = factor(label, unique(regimes$label))) |>
  mutate(panel = factor(paste(label, series, sep = ": "),
                        paste(rep(unique(regimes$label), each = 2), c("reproduction", "recruits"), sep = ": "))) |>
  ggplot(aes(period, count, fill = series)) +
  geom_col(width = 0.7, na.rm = TRUE) +
  facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c(reproduction = "#0072B2", recruits = "#D55E00"), guide = "none") +
  theme_minimal(base_size = 8.5) +
  theme(strip.text = element_text(hjust = 0))

## -----------------------------------------------------------------------------
k4 <- lag_kernels(4)
ceilings <- regimes |>
  group_split(record) |>
  setNames(sort(unique(regimes$record))) |>
  lapply(lag_ceiling, K = 4, kernels = k4)

tibble(record = names(ceilings),
       ceiling = sapply(ceilings, function(x) x$theorem["gini"]),
       observed = sapply(ceilings, function(x) x$obs["gini"]),
       exceedance = sapply(ceilings, function(x) x$exceedance["gini"])) |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))

## ----fig.height = 5-----------------------------------------------------------
rain_strip(ceilings)

## -----------------------------------------------------------------------------
set.seed(1)
cal <- lapply(ceilings, ceiling_calibration, nsim = 3000)
bind_rows(cal, .id = "record") |> select(record, null, phi, p_gini, p_cv)

## -----------------------------------------------------------------------------
set.seed(2)
tests <- regimes |>
  group_split(record) |>
  setNames(sort(unique(regimes$record))) |>
  lapply(host_lag_test, K = 4, nsim = 2000, kernels = k4)

verdict <- bind_rows(lapply(tests, `[[`, "verdict"), .id = "record") |>
  select(record, statistic, max_p) |>
  pivot_wider(names_from = statistic, values_from = max_p, names_prefix = "P_")
verdict

## -----------------------------------------------------------------------------
plot(tests$aseasonal_gated, statistic = "gini")

## -----------------------------------------------------------------------------
set.seed(11)
H <- 20; Tn <- 24; K <- 4; t_R <- 6:24
w_true <- regimes_truth |> filter(mechanism == "delayed") |> pull(value)
base <- round(exp(rnorm(H, log(12), 0.6)))
n_rep <- 12
p_true <- replicate(n_rep, {
  d <- tibble(unit = rep(sprintf("u%02d", 1:H), each = Tn), period = rep(1:Tn, H),
              reproduction = as.vector(t(sapply(base, function(b) rpois(Tn, b)))))
  d <- d |> group_by(unit) |>
    mutate(mu = c(rep(NA, 5), 0.08 * convolve_lag(reproduction, w_true, t_R, 1)),
           recruits = ifelse(is.na(mu), NA, rnbinom(n(), size = 2, mu = ifelse(is.na(mu), 0, mu)))) |>
    ungroup()
  ht <- host_lag_test(d, K = K, nsim = 300, flat = FALSE,
                      kernels = lag_kernels(K, families = "delay", extra = list(truth = w_true)))
  ht$table$p_gini[ht$table$label == "truth"]
})
round(sort(p_true), 2)

