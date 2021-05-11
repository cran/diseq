## ---- include = FALSE---------------------------------------------------------
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
}

## ----setup.libraries----------------------------------------------------------
library(diseq)
library(magrittr)

## ----setup.data---------------------------------------------------------------
nobs <- 2000
tobs <- 5

alpha_d <- -0.1
beta_d0 <- 9.8
beta_d <- c(0.3, -0.2)
eta_d <- c(0.6, -0.1)

alpha_s <- 0.1
beta_s0 <- 5.1
beta_s <- c(0.9)
eta_s <- c(-0.5, 0.2)

gamma <- 1.2
beta_p0 <- 3.1
beta_p <- c(0.8)

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

## ----model.parameters.key-----------------------------------------------------
key_columns <- c("id", "date")

## ----model.parameters.time----------------------------------------------------
time_column <- c("date")

## ----model.parameters.quantity------------------------------------------------
quantity_column <- "Q"

## ----model.parameters.price---------------------------------------------------
price_column <- "P"

## ----model.parameters.specifications------------------------------------------
demand_specification <- paste0(price_column, " + Xd1 + Xd2 + X1 + X2")
supply_specification <- "Xs1 + X1 + X2"
price_specification <- "Xp1"

## ----model.parameters.verbose-------------------------------------------------
verbose <- 2

## ----model.parameters.correlated_shocks---------------------------------------
correlated_shocks <- TRUE

## ----model.constructor--------------------------------------------------------
eqmdl <- new(
  "equilibrium_model",
  key_columns,
  quantity_column, price_column,
  demand_specification, paste0(price_column, " + ", supply_specification),
  stochastic_adjustment_data,
  correlated_shocks = correlated_shocks, verbose = verbose
)
bsmdl <- new(
  "diseq_basic",
  key_columns,
  quantity_column, price_column,
  demand_specification, paste0(price_column, " + ", supply_specification),
  stochastic_adjustment_data,
  correlated_shocks = correlated_shocks, verbose = verbose
)
drmdl <- new(
  "diseq_directional",
  key_columns, time_column,
  quantity_column, price_column,
  demand_specification, supply_specification,
  stochastic_adjustment_data,
  correlated_shocks = correlated_shocks, verbose = verbose
)
damdl <- new(
  "diseq_deterministic_adjustment",
  key_columns, time_column,
  quantity_column, price_column,
  demand_specification, paste0(price_column, " + ", supply_specification),
  stochastic_adjustment_data,
  correlated_shocks = correlated_shocks, verbose = verbose
)
samdl <- new(
  "diseq_stochastic_adjustment",
  key_columns, time_column,
  quantity_column, price_column,
  demand_specification, paste0(price_column, " + ", supply_specification),
  price_specification,
  stochastic_adjustment_data,
  correlated_shocks = correlated_shocks, verbose = verbose
)

## ----estimation.parameters.method---------------------------------------------
optimization_method <- "BFGS"
optimization_controls <- list(REPORT = 10, maxit = 10000, reltol = 1e-6)

## ----estimation.execution-----------------------------------------------------
eqmdl_reg <- estimate(eqmdl, method = "2SLS")
eqmdl_est <- estimate(eqmdl,
  control = optimization_controls, method = optimization_method,
  standard_errors = c("id")
)
bsmdl_est <- estimate(bsmdl,
  control = optimization_controls, method = optimization_method,
  standard_errors = "heteroscedastic"
)
drmdl_est <- estimate(drmdl,
  control = optimization_controls, method = optimization_method,
  standard_errors = "heteroscedastic"
)
damdl_est <- estimate(damdl,
  control = optimization_controls, method = optimization_method,
  standard_errors = c("id")
)
samdl_est <- estimate(samdl,
  control = optimization_controls, method = optimization_method,
  standard_errors = c("id")
)

## ----analysis.effects---------------------------------------------------------
variables <- c(price_column, "Xd1", "Xd2", "X1", "X2", "Xs1")

bsmdl_mme <- sapply(variables,
  function(v) shortage_probability_marginal(bsmdl, bsmdl_est@coef, v),
  USE.NAMES = FALSE
)
drmdl_mme <- sapply(variables,
  function(v) shortage_probability_marginal(drmdl, drmdl_est@coef, v),
  USE.NAMES = FALSE
)
damdl_mme <- sapply(variables,
  function(v) shortage_probability_marginal(damdl, damdl_est@coef, v),
  USE.NAMES = FALSE
)
samdl_mme <- sapply(variables,
  function(v) shortage_probability_marginal(samdl, samdl_est@coef, v),
  USE.NAMES = FALSE
)
bsmdl_mem <- sapply(variables,
  function(v) {
    shortage_probability_marginal(bsmdl, bsmdl_est@coef, v,
      aggregate = "at_the_mean"
    )
  },
  USE.NAMES = FALSE
)
drmdl_mem <- sapply(variables,
  function(v) {
    shortage_probability_marginal(drmdl, drmdl_est@coef, v,
      aggregate = "at_the_mean"
    )
  },
  USE.NAMES = FALSE
)
damdl_mem <- sapply(variables,
  function(v) {
    shortage_probability_marginal(damdl, damdl_est@coef, v,
      aggregate = "at_the_mean"
    )
  },
  USE.NAMES = FALSE
)
samdl_mem <- sapply(variables,
  function(v) {
    shortage_probability_marginal(samdl, samdl_est@coef, v,
      aggregate = "at_the_mean"
    )
  },
  USE.NAMES = FALSE
)

cbind(
  bsmdl_mme, drmdl_mme, damdl_mme, samdl_mme, bsmdl_mem, drmdl_mem, damdl_mem,
  samdl_mem
)

## ----analysis.estimates-------------------------------------------------------
mdt <- tibble::add_column(
  bsmdl@model_tibble,
  normalized_shortages = c(normalized_shortages(bsmdl, bsmdl_est@coef)),
  shortage_probabilities = c(shortage_probabilities(bsmdl, bsmdl_est@coef)),
  relative_shortages = c(relative_shortages(bsmdl, bsmdl_est@coef))
)

## ----analysis.shortages-------------------------------------------------------
abs_estsep <- c(
  nobs = length(shortage_indicators(bsmdl, bsmdl_est@coef)),
  nshortages = sum(shortage_indicators(bsmdl, bsmdl_est@coef)),
  nsurpluses = sum(!shortage_indicators(bsmdl, bsmdl_est@coef))
)
print(abs_estsep)

rel_estsep <- abs_estsep / abs_estsep["nobs"]
names(rel_estsep) <- c("total", "shortages_share", "surpluses_share")
print(rel_estsep)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(mdt, ggplot2::aes(normalized_shortages)) +
    ggplot2::geom_density() +
    ggplot2::ggtitle(paste0(
      "Normalized shortages density (",
      model_name(bsmdl), ")"
    ))
}

## ----analysis.summaries-------------------------------------------------------
summary(eqmdl_reg$first_stage_model)
summary(eqmdl_reg$system_model)
bbmle::summary(eqmdl_est)
bbmle::summary(bsmdl_est)
bbmle::summary(damdl_est)
bbmle::summary(samdl_est)

## ----analysis.market_forces---------------------------------------------------
market <- cbind(
  demand = demanded_quantities(bsmdl, bsmdl_est@coef)[, 1],
  supply = supplied_quantities(bsmdl, bsmdl_est@coef)[, 1]
)
summary(market)

## ----analysis.aggregation-----------------------------------------------------
aggregates <- c(
  demand = aggregate_demand(bsmdl, bsmdl_est@coef),
  supply = aggregate_supply(bsmdl, bsmdl_est@coef)
)
aggregates

