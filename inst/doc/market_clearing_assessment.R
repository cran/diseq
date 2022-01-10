## ---- include = FALSE---------------------------------------------------------
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
}

## ----setup.libraries----------------------------------------------------------
library(diseq)
library(magrittr)

## ----setup.data---------------------------------------------------------------
nobs <- 1000
tobs <- 5

alpha_d <- -3.9
beta_d0 <- 18.9
beta_d <- c(2.1, -0.7)
eta_d <- c(3.5, 6.25)

alpha_s <- 2.8
beta_s0 <- 3.2
beta_s <- c(2.65)
eta_s <- c(1.15, 4.2)

sigma_d <- 0.8
sigma_s <- 1.1
rho_ds <- 0.0

seed <- 42

eq_data <- simulate_data(
  "equilibrium_model", nobs, tobs,
  alpha_d, beta_d0, beta_d, eta_d,
  alpha_s, beta_s0, beta_s, eta_s,
  NA, NA, c(NA),
  sigma_d = sigma_d, sigma_s = sigma_s, rho_ds = rho_ds,
  seed = seed
)

## ----model.parameters---------------------------------------------------------
verbose <- 2
correlated_shocks <- TRUE
formula <-   Q | P | id | date ~ P + Xd1 + Xd2 + X1 + X2 | P + Xs1 + X1 + X2

## ----estimation.parameters----------------------------------------------------
optimization_method <- "BFGS"
optimization_options <- list(maxit = 10000, reltol = 1e-8)

## ----model.constructor--------------------------------------------------------
eqmdl_reg <- equilibrium_model(
  formula, eq_data[eq_data$date != 1, ],
  correlated_shocks = correlated_shocks, verbose = verbose,
  estimation_options = list(method = "2SLS")
)
eqmdl_fit <- equilibrium_model(
  formula, eq_data[eq_data$date != 1, ],
  correlated_shocks = correlated_shocks, verbose = verbose,
  estimation_options = list(
    control = optimization_options, method = optimization_method
  )
)
bsmdl_fit <- diseq_basic(
  formula, eq_data[eq_data$date != 1, ],
  correlated_shocks = correlated_shocks, verbose = verbose,
  estimation_options = list(
    control = optimization_options, method = optimization_method
  )
)
damdl_fit <- diseq_deterministic_adjustment(
  formula, eq_data,
  correlated_shocks = correlated_shocks, verbose = verbose,
  estimation_options = list(
    control = optimization_options, method = optimization_method
  )
)

## ----analysis.summaries-------------------------------------------------------
summary(eqmdl_reg@fit[[1]]$first_stage_model)
summary(eqmdl_reg)
summary(eqmdl_fit)
summary(bsmdl_fit)
summary(damdl_fit)

## ----analysis.estimates-------------------------------------------------------
sim_coef <- c(
  alpha_d, beta_d0, beta_d, eta_d,
  alpha_s, beta_s0, beta_s, eta_s,
  NA,
  sigma_d, sigma_s,
  rho_ds
)
names(sim_coef) <- names(coef(damdl_fit))

dm_inc <- coef(eqmdl_reg)[
  grep("demand", names(coef(eqmdl_reg)))
]
sp_inc <- coef(eqmdl_reg)[
  grep("supply", names(coef(eqmdl_reg)))
]
lm_coef <- c(
  dm_inc[2], dm_inc[-2], sp_inc[2], sp_inc[-2],
  NA,
  NA, NA,
  NA
)

eqmdl_coef <- append(
  coef(eqmdl_fit), c(NA),
  after = which(names(coef(eqmdl_fit)) ==
    prefixed_variance_variable(eqmdl_fit@system@demand)) - 1
)

bsmdl_coef <- append(
  coef(bsmdl_fit), c(NA),
  after = which(names(coef(bsmdl_fit)) ==
    prefixed_variance_variable(bsmdl_fit@system@demand)) - 1
)

damdl_coef <- coef(damdl_fit)

comp <- tibble::tibble(
  parameter = names(sim_coef),
  sim = sim_coef, lm = lm_coef, fi = eqmdl_coef,
  bm = bsmdl_coef, da = damdl_coef,
  lmerr = abs(lm_coef - sim_coef), fierr = abs(eqmdl_coef - sim_coef),
  bmerr = abs(bsmdl_coef - sim_coef), daerr = abs(damdl_coef - sim_coef)
)
comp

## ----analysis.averages--------------------------------------------------------
comp_means <- colMeans(comp[, grep("err", colnames(comp))], na.rm = TRUE)
comp_means

## ----analysis.model.selection-------------------------------------------------
model_names <- c(
  eqmdl_fit@model_type_string,
  bsmdl_fit@model_type_string, damdl_fit@model_type_string
)
model_obs <- c(nobs(eqmdl_fit), nobs(bsmdl_fit), nobs(damdl_fit))
model_errors <- c(
  comp_means["fierr"],
  comp_means["bmerr"],
  comp_means["daerr"]
)
seltbl <- AIC(eqmdl_fit, bsmdl_fit, damdl_fit) %>%
  tibble::add_column(Model = model_names, .before = 1) %>%
  tibble::add_column(Obs. = model_obs, `Mean Error` = model_errors) %>%
  dplyr::rename(D.F. = df) %>%
  dplyr::arrange(AIC)
seltbl

