## ---- include = FALSE---------------------------------------------------------
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
}

## ----setup.libraries----------------------------------------------------------
library(diseq)
library(Formula)

## ----setup.data---------------------------------------------------------------
nobs <- 2000
tobs <- 5

alpha_d <- -0.3
beta_d0 <- 6.8
beta_d <- c(0.3, -0.02)
eta_d <- c(0.6, -0.1)

alpha_s <- 0.6
beta_s0 <- 2.1
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

seed <- 443

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
market_spec <-   Q | P | id | date ~ P + Xd1 + Xd2 + X1 + X2 | P + Xs1 + X1 + X2

## ----model.constructor--------------------------------------------------------
eqmdl_reg <- equilibrium_model(
  market_spec, stochastic_adjustment_data,
  estimation_options = list(method = "2SLS")
)
eqmdl_fit <- equilibrium_model(
    market_spec, stochastic_adjustment_data,
    estimation_options = list(standard_errors = c("id"))
)
bsmdl_fit <- diseq_basic(
    market_spec, stochastic_adjustment_data,
    estimation_options = list(
        method = "Nelder-Mead", control = list(maxit = 1e+5),
        standard_errors = "heteroscedastic"
    )
)
drmdl_fit <- diseq_directional(
    formula(update(Formula(market_spec), . ~ . | . - P)),
    stochastic_adjustment_data,
    estimation_options = list(standard_errors = "heteroscedastic")
)
damdl_fit <- diseq_deterministic_adjustment(
    market_spec, stochastic_adjustment_data,
    estimation_options = list(standard_errors = c("id"))
)
samdl_fit <- diseq_stochastic_adjustment(
    formula(update(Formula(market_spec), . ~ . | . | Xp1)),
    stochastic_adjustment_data,
    estimation_options = list(control = list(maxit = 1e+5))
)

## ----analysis.effects---------------------------------------------------------
variables <- c("P", "Xd1", "Xd2", "X1", "X2", "Xs1")

bsmdl_mme <- sapply(variables,
  function(v) shortage_probability_marginal(bsmdl_fit, v),
  USE.NAMES = FALSE
)
drmdl_mme <- sapply(variables,
  function(v) shortage_probability_marginal(drmdl_fit, v),
  USE.NAMES = FALSE
)
damdl_mme <- sapply(variables,
  function(v) shortage_probability_marginal(damdl_fit, v),
  USE.NAMES = FALSE
)
samdl_mme <- sapply(variables,
  function(v) shortage_probability_marginal(samdl_fit, v),
  USE.NAMES = FALSE
)
bsmdl_mem <- sapply(variables,
  function(v) {
    shortage_probability_marginal(bsmdl_fit, v, aggregate = "at_the_mean")
  },
  USE.NAMES = FALSE
)
drmdl_mem <- sapply(variables,
  function(v) {
    shortage_probability_marginal(drmdl_fit, v, aggregate = "at_the_mean")
  },
  USE.NAMES = FALSE
)
damdl_mem <- sapply(variables,
  function(v) {
    shortage_probability_marginal(damdl_fit, v, aggregate = "at_the_mean")
  },
  USE.NAMES = FALSE
)
samdl_mem <- sapply(variables,
  function(v) {
    shortage_probability_marginal(samdl_fit, v, aggregate = "at_the_mean")
  },
  USE.NAMES = FALSE
)

cbind(
  bsmdl_mme, drmdl_mme, damdl_mme, samdl_mme, bsmdl_mem, drmdl_mem, damdl_mem,
  samdl_mem
)

## ----analysis.estimates-------------------------------------------------------
mdt <- tibble::add_column(
  bsmdl_fit@model_tibble,
  normalized_shortages = c(normalized_shortages(bsmdl_fit)),
  shortage_probabilities = c(shortage_probabilities(bsmdl_fit)),
  relative_shortages = c(relative_shortages(bsmdl_fit))
)

## ----analysis.shortages-------------------------------------------------------
abs_fitsep <- c(
  nobs = length(shortage_indicators(bsmdl_fit)),
  nshortages = sum(shortage_indicators(bsmdl_fit)),
  nsurpluses = sum(!shortage_indicators(bsmdl_fit))
)
print(abs_fitsep)

rel_fitsep <- abs_fitsep / abs_fitsep["nobs"]
names(rel_fitsep) <- c("total", "shortages_share", "surpluses_share")
print(rel_fitsep)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(mdt, ggplot2::aes(shortage_probabilities)) +
    ggplot2::geom_density() +
    ggplot2::ggtitle(paste0(
      "Normalized shortages density (",
      model_name(bsmdl_fit), ")"
    )) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "transparent"),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
    )
}

## ----analysis.summaries-------------------------------------------------------
summary(eqmdl_reg@fit[[1]]$first_stage_model)
summary(eqmdl_reg)
summary(eqmdl_fit)
summary(bsmdl_fit)
summary(damdl_fit)
summary(samdl_fit)

## ----analysis.market_forces---------------------------------------------------
market <- cbind(
  demand = demanded_quantities(bsmdl_fit)[, 1],
  supply = supplied_quantities(bsmdl_fit)[, 1]
)
summary(market)

## ----analysis.aggregation-----------------------------------------------------
aggregates <- c(
  demand = aggregate_demand(bsmdl_fit),
  supply = aggregate_supply(bsmdl_fit)
)
aggregates

