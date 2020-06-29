## ---- include = FALSE---------------------------------------------------------
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
}

## ----setup.libraries----------------------------------------------------------
library(diseq)
library(magrittr)

## ----setup.data---------------------------------------------------------------
nobs <- 5000
tobs <- 5

alpha_d <- -1.9
beta_d0 <- 4.9
beta_d <- c(2.1, -0.7)
eta_d <- c(3.5, 6.25)

alpha_s <- 2.8
beta_s0 <- 1.2
beta_s <- c(0.65)
eta_s <- c(1.15, 4.2)

sigma_d <- 1
sigma_s <- 1
rho_ds <- 0.5

seed <- 443

eq_data <- simulate_model_data(
  "eq_fiml", nobs, tobs,
  alpha_d, beta_d0, beta_d, eta_d,
  alpha_s, beta_s0, beta_s, eta_s,
  NA, NA, c(NA),
  sigma_d = sigma_d, sigma_s = sigma_s, rho_ds = rho_ds,
  seed = seed
)

## ----model.parameters---------------------------------------------------------
key_columns <- c("id", "date")
time_column <- c("date")
quantity_column <- "Q"
price_column <- "P"
demand_specification <- paste0(price_column, " + Xd1 + Xd2 + X1 + X2")
supply_specification <- "Xs1 + X1 + X2"
price_specification <- "Xp1"
verbose <- 2
use_correlated_shocks <- TRUE

## ----model.constructor--------------------------------------------------------
eq2sls <- new(
  "eq_2sls",
  key_columns,
  quantity_column, price_column,
  demand_specification, paste0(price_column, " + ", supply_specification),
  eq_data[eq_data$date != 1, ],
  verbose = verbose
)
eqfiml <- new(
  "eq_fiml",
  key_columns,
  quantity_column, price_column,
  demand_specification, paste0(price_column, " + ", supply_specification),
  eq_data[eq_data$date != 1, ],
  use_correlated_shocks = use_correlated_shocks, verbose = verbose
)
bsmdl <- new(
  "diseq_basic",
  key_columns,
  quantity_column, price_column,
  demand_specification, paste0(price_column, " + ", supply_specification),
  eq_data[eq_data$date != 1, ],
  use_correlated_shocks = use_correlated_shocks, verbose = verbose
)
damdl <- new(
  "diseq_deterministic_adjustment",
  key_columns, time_column,
  quantity_column, price_column,
  demand_specification, paste0(price_column, " + ", supply_specification),
  eq_data,
  use_correlated_shocks = use_correlated_shocks, verbose = verbose
)

## ----estimation.parameters----------------------------------------------------
optimization_method <- "BFGS"
optimization_controls <- list(maxit = 10000, reltol = 1e-8)

## ----estimation.execution-----------------------------------------------------
eq2sls <- estimate(eq2sls)
eqfiml_est <- estimate(eqfiml,
  control = optimization_controls,
  method = optimization_method
)
bsmdl_est <- estimate(bsmdl,
  control = optimization_controls,
  method = optimization_method
)
damdl_est <- estimate(damdl,
  control = optimization_controls,
  method = optimization_method
)

## ----analysis.summaries-------------------------------------------------------
summary(eq2sls@first_stage_model)
summary(eq2sls@system_model)
bbmle::summary(eqfiml_est)
bbmle::summary(bsmdl_est)
bbmle::summary(damdl_est)

## ----analysis.estimates-------------------------------------------------------
sim_coef <- c(
  alpha_d, beta_d0, beta_d, eta_d,
  alpha_s, beta_s0, beta_s, eta_s,
  NA,
  sigma_d, sigma_s,
  rho_ds
)
names(sim_coef) <- names(damdl_est@coef)

dm_inc <- eq2sls@system_model$coefficients[
  grep(
    "demand",
    names(eq2sls@system_model$coefficients)
  )
]
sp_inc <- eq2sls@system_model$coefficients[
  grep(
    "supply",
    names(eq2sls@system_model$coefficients)
  )
]
lm_coef <- c(
  dm_inc[2], dm_inc[-2], sp_inc[2], sp_inc[-2],
  NA,
  NA, NA,
  NA
)

eqfiml_coef <- append(
  eqfiml_est@coef, c(NA),
  after = which(names(eqfiml_est@coef) ==
    get_prefixed_variance_variable(eqfiml@system@demand)) - 1
)

bsmdl_coef <- append(
  bsmdl_est@coef, c(NA),
  after = which(names(bsmdl_est@coef) ==
    get_prefixed_variance_variable(bsmdl@system@demand)) - 1
)

damdl_coef <- damdl_est@coef

comp <- tibble::tibble(
  parameter = names(sim_coef),
  sim = sim_coef, lm = lm_coef, fi = eqfiml_coef,
  bm = bsmdl_coef, da = damdl_coef,
  lmerr = abs(lm_coef - sim_coef), fierr = abs(eqfiml_coef - sim_coef),
  bmerr = abs(bsmdl_coef - sim_coef), daerr = abs(damdl_coef - sim_coef)
)
comp

## ----analysis.averages--------------------------------------------------------
comp_means <- colMeans(comp[, grep("err", colnames(comp))], na.rm = TRUE)
comp_means

## ----analysis.model.selection-------------------------------------------------
model_names <- c(
  eqfiml@model_type_string,
  bsmdl@model_type_string, damdl@model_type_string
)
model_obs <- c(
  get_number_of_observations(eqfiml),
  get_number_of_observations(bsmdl),
  get_number_of_observations(damdl)
)
model_errors <- c(
  comp_means["fierr"],
  comp_means["bmerr"],
  comp_means["daerr"]
)
seltbl <- AIC(eqfiml_est, bsmdl_est, damdl_est) %>%
  tibble::add_column(Model = model_names, .before = 1) %>%
  tibble::add_column(Obs. = model_obs, `Mean Error` = model_errors) %>%
  dplyr::rename(D.F. = df) %>%
  dplyr::arrange(AIC)
seltbl

