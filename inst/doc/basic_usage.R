## ---- include = FALSE---------------------------------------------------------
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
}

## ----setup.libraries----------------------------------------------------------
library(diseq)
library(magrittr)
library(Formula)

## ----setup.data---------------------------------------------------------------
nobs <- 1000
tobs <- 10

alpha_d <- -0.3
beta_d0 <- 6.8
beta_d <- c(0.3, -0.02)
eta_d <- c(0.6, -0.1)

alpha_s <- 0.6
beta_s0 <- 4.1
beta_s <- c(0.9)
eta_s <- c(-0.5, 0.2)

gamma <- 1.2
beta_p0 <- 0.9
beta_p <- c(-0.1)

sigma_d <- 1
sigma_s <- 1
sigma_p <- 1
rho_ds <- 0.0
rho_dp <- 0.0
rho_sp <- 0.0

seed <- 4430

stochastic_adjustment_data <- simulate_data(
  "diseq_stochastic_adjustment", nobs, tobs,
  alpha_d, beta_d0, beta_d, eta_d,
  alpha_s, beta_s0, beta_s, eta_s,
  gamma, beta_p0, beta_p,
  sigma_d = sigma_d, sigma_s = sigma_s, sigma_p = sigma_p,
  rho_ds = rho_ds, rho_dp = rho_dp, rho_sp = rho_sp,
  seed = seed
)

## ----model.parameters---------------------------------------------------------
market_spec <- Q | P | id | date ~ P + Xd1 + Xd2 + X1 + X2 | P + Xs1 + X1 + X2

## ----model.constructor--------------------------------------------------------
eq_reg <- equilibrium_model(
  market_spec, stochastic_adjustment_data,
  estimation_options = list(method = "2SLS")
)
eq_fit <- equilibrium_model(
  market_spec, stochastic_adjustment_data,
  estimation_options = list(standard_errors = c("id"))
)
bs_fit <- diseq_basic(
  market_spec, stochastic_adjustment_data,
  estimation_options = list(
    method = "Nelder-Mead", control = list(maxit = 1e+5),
    standard_errors = "heteroscedastic"
  )
)
dr_fit <- diseq_directional(
  formula(update(Formula(market_spec), . ~ . | . - P)),
  stochastic_adjustment_data,
  estimation_options = list(standard_errors = "heteroscedastic")
)
da_fit <- diseq_deterministic_adjustment(
  market_spec, stochastic_adjustment_data,
  estimation_options = list(standard_errors = c("id"))
)
sa_fit <- diseq_stochastic_adjustment(
  formula(update(Formula(market_spec), . ~ . | . | Xp1)),
  stochastic_adjustment_data,
  estimation_options = list(control = list(maxit = 1e+5))
)

## ----analysis.summaries-------------------------------------------------------
summary(eq_reg@fit[[1]]$first_stage_model)
summary(eq_reg)
summary(eq_fit)
summary(bs_fit)
summary(da_fit)
summary(sa_fit)

## ----analysis.effects---------------------------------------------------------
diseq_abbrs <- c("bs", "dr", "da", "sa")
diseq_fits <- c(bs_fit, dr_fit, da_fit, sa_fit)
variables <- c("P", "Xd1", "Xd2", "X1", "X2", "Xs1")

apply_marginal <- function(fnc, ...) {
  function(fit) {
    sapply(variables, function(v) fnc(fit, v, ...), USE.NAMES = FALSE)
  }
}

mspm <- sapply(diseq_fits, apply_marginal(shortage_probability_marginal))
colnames(mspm) <- diseq_abbrs
# Mean Shortage Probabilities' Marginal Effects
mspm

spmm <- sapply(
  diseq_fits,
  apply_marginal(shortage_probability_marginal, aggregate = "at_the_mean")
)
colnames(spmm) <- diseq_abbrs
# Shortage Probabilities' Marginal Effects at the Mean
spmm

## ----analysis.estimates-------------------------------------------------------
fit <- sa_fit
mdt <- tibble::add_column(
  fit@model_tibble,
  shortage_indicators = c(shortage_indicators(fit)),
  normalized_shortages = c(normalized_shortages(fit)),
  shortage_probabilities = c(shortage_probabilities(fit)),
  relative_shortages = c(relative_shortages(fit))
)

## ----analysis.shortages-------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  pdt <- tibble::tibble(
    Shortage = c(mdt$normalized_shortages, mdt$relative_shortages),
    Type = c(rep("Normalized", nrow(mdt)), rep("Relative", nrow(mdt))),
    xpos = c(rep(-1.0, nrow(mdt)), rep(1.0, nrow(mdt))),
    ypos = c(
      rep(0.8 * max(mdt$normalized_shortages), nrow(mdt)),
      rep(0.8 * max(mdt$relative_shortages), nrow(mdt))
    )
  )
  ggplot2::ggplot(pdt) +
    ggplot2::stat_density(ggplot2::aes(Shortage,
      linetype = Type,
      color = Type
    ), geom = "line") +
    ggplot2::ggtitle(paste0("Estimated shortages densities (", model_name(fit), ")")) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "transparent"),
      plot.background = ggplot2::element_rect(
        fill = "transparent",
        color = NA
      ),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.box.background = ggplot2::element_rect(
        fill = "transparent",
        color = NA
      ),
      legend.position = c(0.8, 0.8)
    )
} else {
    summary(mdt[, grep("shortage", colnames(mdt))])
}

## ----analysis.market_forces---------------------------------------------------
market <- cbind(
  demand = demanded_quantities(fit)[, 1],
  supply = supplied_quantities(fit)[, 1]
)
summary(market)

## ----analysis.aggregation-----------------------------------------------------
aggregates <- aggregate_demand(fit) %>%
  dplyr::left_join(aggregate_supply(fit), by = "date") %>%
  dplyr::mutate(date = as.numeric(date)) %>%
  dplyr::rename(demand = D_Q, supply = S_Q)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  pdt <- tibble::tibble(
    Date = c(aggregates$date, aggregates$date),
    Quantity = c(aggregates$demand, aggregates$supply),
    Side = c(rep("Demand", nrow(aggregates)), rep("Supply", nrow(aggregates)))
  )
  ggplot2::ggplot(pdt, ggplot2::aes(x = Date)) +
    ggplot2::geom_line(ggplot2::aes(y = Quantity, linetype = Side, color = Side)) +
    ggplot2::ggtitle(paste0(
      "Aggregate estimated demand and supply  (", model_name(fit), ")"
    )) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "transparent"),
      plot.background = ggplot2::element_rect(
        fill = "transparent", color = NA
      ),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.box.background = ggplot2::element_rect(
        fill = "transparent", color = NA
      ),
      legend.position = c(0.8, 0.5)
    )
} else {

    aggregates
}

