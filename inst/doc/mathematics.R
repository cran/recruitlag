## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4)

## -----------------------------------------------------------------------------
library(recruitlag)

## -----------------------------------------------------------------------------
X <- lepanthes_census$inflorescences
convolve_lag(X, c(0, 0, 0, 0, 1), t_R = 3:12)          # pure 24 to 30 month delay
convolve_lag(X, flat_kernel(4), t_R = 3:12)             # equal weights

## -----------------------------------------------------------------------------
ce <- lag_ceiling(lepanthes_census, reproduction = "inflorescences", K = 4,
                  kernels = lag_kernels(4, n_dirichlet = 2000, bin = 6, unit = "mo"))
ce$theorem; ce$ceiling; ce$best$label

## -----------------------------------------------------------------------------
s <- ce$search
plot(s$gini, s$cv, pch = 16, cex = 0.5, col = ifelse(s$family == "pure delay", "#D55E00", "#56B4E9"),
     xlab = "Gini of the expected series", ylab = "CV of the expected series",
     xlim = c(0, 0.7), ylim = c(0, 1.5))
abline(v = ce$theorem["gini"], h = ce$theorem["cv"], lty = 2)
points(ce$obs["gini"], ce$obs["cv"], pch = 17, cex = 1.6)
legend("topleft", bty = "n", pch = c(16, 16, 17), col = c("#56B4E9", "#D55E00", "black"),
       legend = c("searched kernels", "pure delays", "observed recruitment"))

## -----------------------------------------------------------------------------
h <- host_matrices(lepanthes_hosts, unit = "host", reproduction = "inflorescences")
mu <- expected_recruits(h$R, h$X, w = c(0, 0, 0, 0, 1))
rowSums(mu)[1:5]; rowSums(h$R)[1:5]     # scaled to each host's total
phi_moment(h$R, mu)

