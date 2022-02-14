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
beta_d0 <- 28.9
beta_d <- c(2.1, -0.7)
eta_d <- c(3.5, 6.25)

alpha_s <- 2.8
beta_s0 <- 26.2
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
eq_reg <- equilibrium_model(
  formula, eq_data[eq_data$date != 1, ],
  correlated_shocks = correlated_shocks, verbose = verbose,
  estimation_options = list(method = "2SLS")
)
eq_fit <- equilibrium_model(
  formula, eq_data[eq_data$date != 1, ],
  correlated_shocks = correlated_shocks, verbose = verbose,
  estimation_options = list(
    control = optimization_options, method = optimization_method
  )
)
bs_fit <- diseq_basic(
  formula, eq_data[eq_data$date != 1, ],
  correlated_shocks = correlated_shocks, verbose = verbose,
  estimation_options = list(
    control = optimization_options, method = optimization_method
  )
)
da_fit <- diseq_deterministic_adjustment(
  formula, eq_data,
  correlated_shocks = correlated_shocks, verbose = verbose,
  estimation_options = list(
    control = optimization_options, method = optimization_method
  )
)

## ----analysis.summaries-------------------------------------------------------
summary(eq_reg@fit[[1]]$first_stage_model)
summary(eq_reg)
summary(eq_fit)
summary(bs_fit)
summary(da_fit)

## ----analysis.estimates-------------------------------------------------------
da_coef <- coef(da_fit)
coef_names <- names(da_coef)

sim_coef <- c(
  alpha_d, beta_d0, beta_d, eta_d,
  alpha_s, beta_s0, beta_s, eta_s,
  NA,
  sigma_d, sigma_s,
  rho_ds
)

coef_tbl <- function(fit) {
    tibble::tibble(coef = names(coef(fit)),!!substitute(fit) := coef(fit))
}

comp <- coef_tbl(da_fit) %>% 
    dplyr::left_join(coef_tbl(bs_fit), by = "coef") %>%
    dplyr::left_join(coef_tbl(eq_reg), by = "coef") %>%
    dplyr::left_join(coef_tbl(eq_fit), by = "coef") %>%
    dplyr::mutate(sim = sim_coef) %>%
    dplyr::mutate(sim = sim_coef,
                  da_fit_err = abs(da_fit - sim),
                  bs_fit_err = abs(bs_fit - sim),
                  eq_fit_err = abs(eq_fit - sim)) 

comp

## ----analysis.averages--------------------------------------------------------
model_errors <- colMeans(comp[, grep("err", colnames(comp))], na.rm = TRUE) 
model_errors

## ----analysis.model.selection-------------------------------------------------
fits <- c(da_fit, bs_fit, eq_fit)
model_names <- sapply(fits, function(m) m@model_type_string)
model_obs <- sapply(fits, nobs)
aic <- sapply(fits, AIC)
df <- sapply(fits, function(m) attr(logLik(m), "df"))
seltbl <- tibble::tibble(Model = model_names, AIC = aic,
                         D.F = df, Obs. = model_obs,
                         Abs.Error = model_errors) %>%
  dplyr::arrange(AIC)
seltbl

