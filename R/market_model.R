#' @include model_logger.R
#' @include system_base.R
#' @importFrom bbmle parnames mle2
#' @importFrom grid grid.raster
#' @importFrom png readPNG
#' @importFrom rlang :=
#' @importFrom stats formula lm model.matrix na.omit median qnorm sd var
#' @import dplyr magrittr tibble

setClassUnion("characterOrNULL", c("character", "NULL"))
setOldClass(c("spec_tbl_df", "tbl_df", "tbl", "data.frame"))
utils::globalVariables("where")

#' @title Market model classes
#'
#' @slot logger Logger object.
#' @slot key_columns Vector of column names that uniquely identify data records. For
#' panel data this vector should contain an entity and a time point identifier.
#' @slot time_column Column name for the time point data.
#' @slot explanatory_columns Vector of explanatory column names for all model's
#' equations.
#' @slot data_columns Vector of model's data column names. This is the union of the
#' quantity, price and explanatory columns.
#' @slot columns Vector of primary key and data column names for all model's equations.
#' @slot model_tibble Model data \code{tibble}.
#' @slot model_type_string Model type string description.
#' @slot system Model's system of equations.
#' @name market_models
#' @seealso initialize_market_model
NULL

#' @describeIn market_models Base class for market models
setClass(
  "market_model",
  representation(
    ## Logging
    logger = "model_logger",

    ## Column fields
    key_columns = "vector",
    time_column = "characterOrNULL",
    explanatory_columns = "vector",
    data_columns = "vector",
    columns = "vector",

    ## Model data
    model_tibble = "tbl_df",
    model_type_string = "character",
    market_type_string = "character",
    system = "system_base"
  )
)

#' @title Model initialization
#'
#' @details
#' The following two subsections describe the common initialization steps of all market
#' model classes.
#'
#' \subsection{Variable construction}{
#' The constructor prepares the model's variables using the passed specifications. The
#' specification strings are expected to follow the syntax of
#' \code{\link[stats]{formula}}. The construction of the model's data uses the variables
#' that are extracted by these specification. The demand variables are extracted by a
#' formula that uses the \code{quantity_column} on the left hand side and the
#' \code{demand_specification} on the right hand side of the formula. The supply
#' variables are constructed by the the\code{quantity_column} and the
#' \code{supply_specification}. In the case of the
#' \code{\linkS4class{diseq_stochastic_adjustment}} model, the price dynamics'
#' variables are extracted using the \code{price_specification}. The
#' \code{price_specification} for the \code{\linkS4class{diseq_stochastic_adjustment}}
#' should contain only terms other than that of excess demand. The excess demand term of
#' the price equation is automatically added by the constructor.
#' }
#'
#' \subsection{Data preparation}{
#'   1. If the passed data set contains rows with NA values, they are dropped. If the
#' verbosity level allows warnings, a warning is emitted reporting how many rows were
#' dropped.
#'
#'   2. After dropping the rows, factor levels may be invalidated. If needed, the
#' constructor readjusts the factor variables by removing the unobserved levels. Factor
#' indicators and interaction terms are automatically created.
#'
#'   3. The primary key column is constructed by pasting the values of the key_columns.
#'
#'   4. In the cases of the \code{\linkS4class{diseq_directional}},
#' \code{\linkS4class{diseq_deterministic_adjustment}}, and
#' the \code{\linkS4class{diseq_stochastic_adjustment}} models, a column with lagged
#' prices is constructed. Since lagged prices are unavailable for the observations of
#' the first time point, these observations are dropped. If the verbosity level allows
#' the emission of information messages, the constructor prints the number of dropped
#' observations.
#'
#'   5. In the cases of the \code{\linkS4class{diseq_directional}}
#' and the \code{\linkS4class{diseq_stochastic_adjustment}} models, a column with price
#' differences is created.
#' }
#'
#' @param .Object The object to be Constructed.
#' @param verbose Verbosity level.
#' @param key_columns Key columns of the data set.
#' @param time_column The time column of the data set.
#' @param quantity_column The quantity variable of the data set.
#' @param price_column The price variable of the data set.
#' @param demand_specification A formula representation of the right hand side of the
#'   demand equation.
#' @param supply_specification A formula representation of the right hand side of the
#'   supply equation.
#' @param price_specification A formula representation of the price equation.
#' @param correlated_shocks Should the model be estimated using correlated shocks?
#' @param data The data set.
#' @return The initialized model.
#' @name initialize_market_model
NULL

setMethod(
  "initialize", "market_model",
  function(.Object, model_type_string, verbose,
           key_columns, time_column, quantity_column, price_column,
           demand_specification, supply_specification, price_specification,
           correlated_shocks,
           data,
           system_initializer) {

    ## Model assignments
    .Object@model_type_string <- model_type_string
    .Object@logger <- new("model_logger", verbose)
    .Object@system@correlated_shocks <- correlated_shocks
    print_info(.Object@logger, "This is '", model_name(.Object), "' model")

    .Object@key_columns <- key_columns
    .Object@time_column <- time_column

    .Object@explanatory_columns <- unique(c(
      all.vars(formula(paste0(quantity_column, " ~ ", demand_specification))),
      all.vars(formula(paste0(quantity_column, " ~ ", supply_specification)))
    ))
    if (.Object@model_type_string %in% c("Stochastic Adjustment")) {
      .Object@explanatory_columns <- unique(c(
        .Object@explanatory_columns,
        all.vars(formula(paste0(price_column, " ~ ", price_specification)))
      ))
    }
    .Object@data_columns <- unique(c(quantity_column, price_column,
                                     .Object@explanatory_columns))
    .Object@columns <- unique(c(.Object@key_columns, .Object@data_columns))

    ## Data assignment
    .Object@model_tibble <- data

    ## Create model tibble
    len <- nrow(.Object@model_tibble)
    .Object@model_tibble <- .Object@model_tibble %>%
      dplyr::select(!!!.Object@columns) %>%
      na.omit()
    drops <- len - nrow(.Object@model_tibble)
    if (drops) {
      print_warning(.Object@logger, "Dropping ", drops, " rows due to omitted values.")
    }

    remove_unused_levels <- function(x) {
      initial_levels <- levels(x)
      x <- factor(x)
      remaining_levels <- levels(x)
      removed_levels <- initial_levels[!(initial_levels %in% remaining_levels)]
      if (length(removed_levels)) {
        print_warning(
          .Object@logger, "Removing unobserved '",
          paste0(removed_levels, collapse = ", "), "' level(s)."
        )
      }
      x
    }
    .Object@model_tibble <- .Object@model_tibble %>%
      dplyr::mutate(dplyr::across(
        where(is.factor),
        remove_unused_levels
      ))

    ## Create primary key column
    key_columns_syms <- rlang::syms(.Object@key_columns)
    .Object@model_tibble <- .Object@model_tibble %>%
      dplyr::mutate(pk = as.integer(paste0(!!!key_columns_syms)))

    ## Do we need to use lags?
    if (.Object@model_type_string %in% c(
      "Directional", "Deterministic Adjustment", "Stochastic Adjustment"
    )) {
      ## Generate lags
      key_syms <- rlang::syms(.Object@key_columns[.Object@key_columns !=
                                                  .Object@time_column])
      price_sym <- rlang::sym(price_column)
      time_sym <- rlang::sym(.Object@time_column)
      lagged_price_column <- paste0("LAGGED_", price_column)
      lagged_price_sym <- rlang::sym(lagged_price_column)

      .Object@model_tibble <- .Object@model_tibble %>%
        dplyr::group_by(!!!key_syms) %>%
        dplyr::mutate(
          !!lagged_price_sym := dplyr::lag(!!price_sym, order_by = !!time_sym)
        ) %>%
        dplyr::ungroup()

      drop_rows <- .Object@model_tibble %>%
        dplyr::select(!!lagged_price_sym) %>%
        is.na() %>%
        c()
      .Object@model_tibble <- .Object@model_tibble[!drop_rows, ]
      print_info(
        .Object@logger, "Dropping ",
        sum(drop_rows), " rows by generating '", lagged_price_column, "'."
      )

      ## Do we need to use first differences?
      if (.Object@model_type_string %in% c("Directional", "Deterministic Adjustment")) {
        ## Generate first differences
        diff_column <- paste0(price_column, "_DIFF")
        diff_sym <- rlang::sym(diff_column)

        .Object@model_tibble <- .Object@model_tibble %>%
          dplyr::group_by(!!!key_syms) %>%
          dplyr::mutate(!!diff_sym := !!price_sym - !!lagged_price_sym) %>%
          dplyr::ungroup()
      }
    }

    if (.Object@model_type_string %in% c("Stochastic Adjustment")) {
      .Object@system <- system_initializer(
        quantity_column, price_column,
        demand_specification, supply_specification, price_specification,
        .Object@model_tibble, correlated_shocks
      )
    }
    else {
      .Object@system <- system_initializer(
        quantity_column, price_column, demand_specification, supply_specification,
        .Object@model_tibble, correlated_shocks
      )
    }

    print_verbose(.Object@logger, "Using columns ",
                  paste0(.Object@columns, collapse = ", "), ".")

    .Object
  }
)

#' Prints a short description of the model.
#'
#' Sends basic information about the model to standard output.
#' @param object A model object.
#' @examples
#' \donttest{
#' model <- simulate_model(
#'   "diseq_stochastic_adjustment", list(
#'     # observed entities, observed time points
#'     nobs = 500, tobs = 3,
#'     # demand coefficients
#'     alpha_d = -0.1, beta_d0 = 9.8, beta_d = c(0.3, -0.2), eta_d = c(0.6, -0.1),
#'     # supply coefficients
#'     alpha_s = 0.1, beta_s0 = 5.1, beta_s = c(0.9), eta_s = c(-0.5, 0.2),
#'     # price equation coefficients
#'     gamma = 1.2, beta_p0 = 3.1, beta_p = c(0.8)
#'   ),
#'   seed = 31
#' )
#'
#' # print short model information
#' show(model)
#' }
#' @rdname show
#' @export
setGeneric("show")

#' @rdname show
setMethod("show", signature(object = "market_model"), function(object) {
  cat(sprintf(
    "\n%s Model for Markets in %s\n",
    object@model_type_string, object@market_type_string
  ))
  show_implementation(object@system)
  cat(sprintf(
    "  %-18s: %s\n", "Shocks",
    ifelse(object@system@correlated_shocks, "Correlated", "Independent")
  ))
})

#' Summarizes the model.
#'
#' Prints basic information about the passed model object. In addition to the output of
#' the \code{\link{show}} method, \code{summary} prints
#' - the number of observations,
#' - the number of observations in each equation for models with sample separation, and
#' - various categories of variables.
#' @param object A model object.
#' @examples
#' \donttest{
#' model <- simulate_model(
#'   "diseq_stochastic_adjustment", list(
#'     # observed entities, observed time points
#'     nobs = 500, tobs = 3,
#'     # demand coefficients
#'     alpha_d = -0.1, beta_d0 = 9.8, beta_d = c(0.3, -0.2), eta_d = c(0.6, -0.1),
#'     # supply coefficients
#'     alpha_s = 0.1, beta_s0 = 5.1, beta_s = c(0.9), eta_s = c(-0.5, 0.2),
#'     # price equation coefficients
#'     gamma = 1.2, beta_p0 = 3.1, beta_p = c(0.8)
#'   ),
#'   seed = 556
#' )
#'
#' # print model summary
#' summary(model)
#' }
#' @rdname summary
#' @export
setMethod("summary", signature(object = "market_model"), function(object) {
  show(object)
  cat(sprintf("  %-18s: %d\n", "Nobs", nrow(object@model_tibble)))
  summary_implementation(object@system)
  cat(sprintf(
    "  %-18s: %s\n", "Key Var(s)",
    paste0(object@key_columns, collapse = ", ")
  ))
  if (!is.null(object@time_column)) {
    cat(sprintf(
      "  %-18s: %s\n", "Time Var",
      paste0(object@time_column, collapse = ", ")
    ))
  }
})

#' Plots the model.
#'
#' Displays a graphical illustration of the passed model object.
#' @param x A model object.
#' @examples
#' \donttest{
#' model <- simulate_model(
#'   "diseq_basic", list(
#'     # observed entities, observed time points
#'     nobs = 500, tobs = 3,
#'     # demand coefficients
#'     alpha_d = -0.9, beta_d0 = 8.9, beta_d = c(0.3, -0.2), eta_d = c(-0.03, -0.01),
#'     # supply coefficients
#'     alpha_s = 0.9, beta_s0 = 4.2, beta_s = c(0.03), eta_s = c(0.05, 0.02)
#'   ),
#'   seed = 44
#' )
#'
#' # show model's illustration plot
#' plot(model)
#' }
#' @rdname plot
#' @export
setMethod("plot", signature(x = "market_model"), function(x) {
  filename <- paste0(class(x)[1], ".png")
  path <- system.file("help", "figures", filename, package = "diseq")
  if (path == "") {
    path <- system.file("man", "figures", filename, package = "diseq")
  }
  grid::grid.raster(png::readPNG(path))
})

#' Minus log-likelihood.
#'
#' Returns the opposite of the log-likelihood. The likelihood functions are based on
#' Maddala and Nelson (1974) \doi{10.2307/1914215}. The likelihoods expressions
#' that the function uses are derived in
#' Karapanagiotis (2020) \doi{10.2139/ssrn.3525622}. The function calculates
#' the model's log likelihood by evaluating the log likelihood of each observation in
#' the sample and summing the evaluation results.
#' @param object A model object.
#' @param parameters A vector of parameters at which the function is to be evaluated.
#' @return The opposite of the sum of the likelihoods evaluated for each observation.
#' @rdname minus_log_likelihood
#' @export
setGeneric("minus_log_likelihood", function(object, parameters) {
  standardGeneric("minus_log_likelihood")
})

#' Gradient
#'
#' Returns the gradient of the opposite of the log-likelihood evaluated at the passed
#' parameters.
#' @param object A model object.
#' @param parameters A vector of parameters at which the gradient is to be evaluated.
#' @return The opposite of the model log likelihood's gradient.
#' @rdname gradient
#' @export
setGeneric("gradient", function(object, parameters) {
  standardGeneric("gradient")
})

#' Hessian
#'
#' Returns the hessian of the opposite of the log-likelihood evaluated at the passed
#' parameters.
#' @param object A model object.
#' @param parameters A vector of parameters at which the hessian is to be evaluated.
#' @return The opposite of the model log likelihood's hessian.
#' @rdname hessian
#' @export
setGeneric("hessian", function(object, parameters) {
  standardGeneric("hessian")
})

validate_gradient_option <- function(object, option) {
  allowed <- c("calculated", "numerical")
  if (!(option %in% allowed)) {
    print_error(
      object@logger,
      paste0(
        "Invalid `gradient` option '", option, "'. Valid options are ('",
        paste0(allowed, collapse = "', '"), "')."
      )
    )
  }
}

validate_hessian_option <- function(object, option) {
  allowed <- c("skip", "calculated", "numerical")
  if (!(option %in% allowed)) {
    print_error(
      object@logger,
      paste0(
        "Invalid `hessian` option '", option, "'. Valid options are ('",
        paste0(allowed, collapse = "', '"), "')."
      )
    )
  }
}

validate_standard_error_option <- function(object, option) {
  allowed <- c("homoscedastic", "heteroscedastic")
  if (!(option %in% allowed || all(option %in% object@columns))) {
    print_error(
      object@logger,
      paste0(
        "Invalid `standard_error` option '", option, "'. Valid options are ('",
        paste0(allowed, collapse = "', '"), "') or a vector of model variable names."
      )
    )
  }
}

#' Model estimation.
#'
#' All models are estimated using full information maximum likelihood. The
#' \code{\linkS4class{equilibrium_model}} can also be estimated using two-stage
#' least squares. The maximum likelihood estimation is based on
#' \code{\link[bbmle]{mle2}}. If no starting values are provided, the function uses
#' linear regression estimates as initializing values. The default optimization method is
#' BFGS. For other alternatives see \code{\link[bbmle]{mle2}}. The implementation of
#' the two-stage least square estimation of the \code{\linkS4class{equilibrium_model}}
#' is based on \code{\link[systemfit]{systemfit}}.
#' @param object A model object.
#' @param ... Named parameter used in the model's estimation. These are passed further
#' down to the estimation call. For the \code{\linkS4class{equilibrium_model}} model, the
#' parameters are passed to \code{\link[systemfit]{systemfit}}, if the method is set to
#' \code{2SLS}, or to \code{\link[bbmle]{mle2}} for any other method. For the rest of
#' the models, the parameters are passed to \code{\link[bbmle]{mle2}}.
#' @return The object that holds the estimation result.
#' @rdname estimate
#' @examples
#' \donttest{
#' # initialize the model using the houses dataset
#' model <- new(
#'   "diseq_deterministic_adjustment", # model type
#'   c("ID", "TREND"), "TREND", "HS", "RM", # keys, time, quantity, and price variables
#'   "RM + TREND + W + CSHS + L1RM + L2RM + MONTH", # demand specification
#'   "RM + TREND + W + L1RM + MA6DSF + MA3DHF + MONTH", # supply specification
#'   fair_houses(), # data
#'   correlated_shocks = FALSE # allow shocks to be correlated
#' )
#'
#' # estimate the model object (BFGS is used by default)
#' est <- estimate(model)
#'
#' # estimate the model by specifying the optimization details passed to the optimizer.
#' est <- estimate(model, control = list(maxit = 1e+5), method = "BFGS")
#'
#' # summarize results
#' bbmle::summary(est)
#' }
#' @export
setGeneric("estimate", function(object, ...) {
  standardGeneric("estimate")
})

#' @describeIn estimate Full information maximum likelihood estimation.
#' @param gradient One of two potential options: `numerical` and `calculated`. By
#' default, all the models are estimated using the analytic expressions of their
#' likelihoods' gradients.
#' @param hessian One of three potential options: `skip`, `numerical`, and `calculated`.
#' The default is to use the `calculated` Hessian for the model that expressions are
#' available and the `numerical` Hessian in other cases. Calculated Hessian expressions
#' are available for the basic and directional models.
#' @param standard_errors One of three potential options: `homoscedastic`,
#' `heteroscedastic`, or a vector with variables names for which standard error
#' clusters are to be created. The default value is `homoscedastic`. If the option
#' `heteroscedastic` is passed, the variance-covariance matrix is calculated using
#' heteroscedasticity adjusted (Huber-White) standard errors. If the vector is
#' supplied, the variance-covariance matrix is calculated by grouping the score matrix
#' based on the passed variables.
setMethod(
  "estimate", signature(object = "market_model"),
  function(object, gradient = "calculated", hessian = "calculated",
           standard_errors = "homoscedastic", ...) {
    validate_gradient_option(object, gradient)
    validate_hessian_option(object, hessian)
    validate_standard_error_option(object, standard_errors)

    va_args <- list(...)

    if (hessian == "skip" ||
      ((object@model_type_string %in% c("Basic", "Directional")) &&
        hessian == "calculated")) {
      va_args$skip.hessian <- TRUE
    } else {
      hessian <- "numerical"
    }

    va_args$start <- prepare_initializing_values(object, va_args$start)

    if (is.null(va_args$method)) {
      va_args$method <- "BFGS"
    }

    va_args$minuslogl <- function(...) minus_log_likelihood(object, ...)
    bbmle::parnames(va_args$minuslogl) <- likelihood_variables(object@system)
    if (gradient == "calculated") {
      va_args$gr <- function(...) gradient(object, ...)
      bbmle::parnames(va_args$gr) <- likelihood_variables(object@system)
    }

    est <- do.call(bbmle::mle2, va_args)
    est@call.orig <- call("bbmle::mle2", va_args)

    if (hessian == "calculated") {
      print_verbose(object@logger, "Calculating hessian and variance-covariance matrix.")
      est@details$hessian <- hessian(object, est@coef)
      tryCatch(
        est@vcov <- MASS::ginv(est@details$hessian),
        error = function(e) print_warning(object@logger, e$message)
      )
    }

    if (length(standard_errors) == 1) {
      if (standard_errors == "heteroscedastic") {
        est <- set_heteroscedasticity_consistent_errors(object, est)
      } else if (standard_errors != "homoscedastic") {
        est <- set_clustered_errors(object, est, standard_errors)
      }
    } else {
      est <- set_clustered_errors(object, est, standard_errors)
    }

    est
  }
)


#' Maximize the log-likelihood.
#'
#' Maximizes the log-likelihood using the
#' \href{https://www.gnu.org/software/gsl/doc/html/multimin.html}{\code{GSL}}
#' implementation of the BFGS algorithm. This function is primarily intended for
#' advanced usage. The \code{\link{estimate}} functionality is a fast,
#' analysis-oriented alternative. If the
#' \href{https://www.gnu.org/software/gsl/doc/html/multimin.html}{\code{GSL}} is not
#' available, the function returns a trivial result list with status set equal to -1.
#' If the
#' \href{https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag_t}{C++17
#' execution policies}
#' are available, the implementation of the optimization is parallelized.
#' @param object A model object.
#' @param start Initializing vector.
#' @param step Optimization step.
#' @param objective_tolerance Objective optimization tolerance.
#' @param gradient_tolerance Gradient optimization tolerance.
#' @return A list with the optimization output.
#' @rdname maximize_log_likelihood
#' @seealso estimate
#' @examples
#' \donttest{
#' model <- simulate_model(
#'   "equilibrium_model", list(
#'     # observed entities, observed time points
#'     nobs = 500, tobs = 3,
#'     # demand coefficients
#'     alpha_d = -0.9, beta_d0 = 14.9, beta_d = c(0.3, -0.2), eta_d = c(-0.03, -0.01),
#'     # supply coefficients
#'     alpha_s = 0.9, beta_s0 = 3.2, beta_s = c(0.03), eta_s = c(0.05, 0.02)
#'   ),
#'   seed = 99
#' )
#'
#' # maximize the model's log-likelihood
#' mll <- maximize_log_likelihood(
#'   model,
#'   start = NULL, step = 1e-5,
#'   objective_tolerance = 1e-4, gradient_tolerance = 1e-3
#' )
#' mll
#' }
#' @export
setGeneric("maximize_log_likelihood", function(object, start, step, objective_tolerance,
                                               gradient_tolerance) {
  standardGeneric("maximize_log_likelihood")
})

#' Likelihood scores.
#'
#' It calculates the gradient of the likelihood at the given parameter point for each
#' observation in the sample. It, therefore, returns an n x k matrix, where n denotes
#' the number of observations in the sample and k the number of estimated parameters.
#' There order of the parameters is the same as the one that is used in the summary
#' of the results.
#' @param object A model object.
#' @param parameters A vector with model parameters.
#' @return The score matrix.
#' @rdname scores
#' @examples
#' \donttest{
#' model <- simulate_model(
#'   "diseq_basic", list(
#'     # observed entities, observed time points
#'     nobs = 500, tobs = 3,
#'     # demand coefficients
#'     alpha_d = -0.9, beta_d0 = 8.9, beta_d = c(0.6), eta_d = c(-0.2),
#'     # supply coefficients
#'     alpha_s = 0.9, beta_s0 = 4.2, beta_s = c(0.03, 1.2), eta_s = c(0.1)
#'   ),
#'   seed = 7523
#' )
#'
#' # estimate the model object (BFGS is used by default)
#' est <- estimate(model)
#'
#' # Calculate the score matrix
#' head(scores(model, est@coef))
#' }
#' @export
setGeneric("scores", function(object, parameters) {
  standardGeneric("scores")
})

setGeneric("set_heteroscedasticity_consistent_errors", function(object, ...) {
  standardGeneric("set_heteroscedasticity_consistent_errors")
})

setGeneric("set_clustered_errors", function(object, ...) {
  standardGeneric("set_clustered_errors")
})

#' Model description.
#'
#' A unique identifying string for the model.
#' @param object A model object.
#' @return A string representation of the model.
#' @rdname model_name
#' @export
setGeneric("model_name", function(object) {
  standardGeneric("model_name")
})

#' Number of observations.
#'
#' Returns the number of observations that are used by an initialized model. The number
#' of used observations may differ from the numbers of observations of the data set
#' that was passed to the model's initialization.
#' @param object A model object.
#' @return The number of used observations.
#' @rdname number_of_observations
#' @export
setGeneric("number_of_observations", function(object) {
  standardGeneric("number_of_observations")
})

#' @title Market side descriptive statistics
#' @details Calculates and returns basic descriptive statistics for the model's demand
#' or supply side data. Factor variables are excluded from the calculations. The function
#' calculates and returns:
#' \itemize{
#' \item \code{nobs} Number of observations.
#' \item \code{nmval} Number of missing values.
#' \item \code{min} Minimum observation.
#' \item \code{max} Maximum observation.
#' \item \code{range} Observations' range.
#' \item \code{sum} Sum of observations.
#' \item \code{median} Median observation.
#' \item \code{mean} Mean observation.
#' \item \code{mean_se} Mean squared error.
#' \item \code{mean_ce} Confidence interval bound.
#' \item \code{var} Variance.
#' \item \code{sd} Standard deviation.
#' \item \code{coef_var} Coefficient of variation.
#' }
#' @param object A model object.
#' @return A data \code{tibble} containing descriptive statistics.
#' @examples
#' # initialize the basic model using the houses dataset
#' model <- new(
#'   "diseq_basic", # model type
#'   c("ID", "TREND"), "HS", "RM", # keys, quantity, and price variables
#'   "RM + TREND + W + CSHS + L1RM + L2RM + MONTH", # demand specification
#'   "RM + TREND + W + L1RM + MA6DSF + MA3DHF + MONTH", # supply specification
#'   fair_houses(), # data
#'   correlated_shocks = FALSE # allow shocks to be correlated
#' )
#'
#' # get descriptive statistics of demand side variables
#' demand_descriptives(model)
#'
#' # get descriptive statistics of supply side variables
#' supply_descriptives(model)
#' @name market_descriptives
NULL

setGeneric("descriptives", function(object, variables) {
  standardGeneric("descriptives")
})

#' @describeIn market_descriptives Demand descriptive statistics.
#' @export
setGeneric("demand_descriptives", function(object) {
  standardGeneric("demand_descriptives")
})

#' @describeIn market_descriptives Supply descriptive statistics.
#' @export
setGeneric("supply_descriptives", function(object) {
  standardGeneric("supply_descriptives")
})

setMethod(
  "set_heteroscedasticity_consistent_errors", signature(object = "market_model"),
  function(object, est) {
    est@details$original_hessian <- est@details$hessian
    scores <- scores(object, est@coef)
    adjustment <- MASS::ginv(t(scores) %*% scores)
    est@details$hessian <- est@details$hessian %*% adjustment %*% est@details$hessian
    est@vcov <- MASS::ginv(est@details$hessian)
    est
  }
)

setMethod(
  "set_clustered_errors", signature(object = "market_model"),
  function(object, est, cluster_errors_by) {
    if (!(cluster_errors_by %in% names(object@model_tibble))) {
      print_error(
        object@logger, "Cluster variable is not among model data variables."
      )
    }
    est@details$original_hessian <- est@details$hessian
    cluster_var <- rlang::syms(cluster_errors_by)
    clustered_scores <- tibble::tibble(
      object@model_tibble %>% dplyr::select(!!!cluster_var),
      tibble::as_tibble(scores(object, est@coef))
    ) %>%
      dplyr::group_by(!!!cluster_var) %>%
      dplyr::group_map(~t(as.matrix(.)) %*% (as.matrix(.)))
    est@details$number_of_clusters <- length(clustered_scores)
    adjustment <- MASS::ginv(Reduce("+", clustered_scores))
    est@details$hessian <- est@details$hessian %*% adjustment %*% est@details$hessian
    est@vcov <- MASS::ginv(est@details$hessian) * (
      est@details$number_of_clusters / (est@details$number_of_clusters - 1)
    )
    est
  }
)

#' @rdname model_name
setMethod("model_name", signature(object = "market_model"), function(object) {
  paste0(
    object@model_type_string, " with ",
    ifelse(object@system@correlated_shocks, "correlated", "independent"), " shocks"
  )
})

#' @rdname number_of_observations
setMethod("number_of_observations", signature(object = "market_model"),
          function(object) {
  nrow(object@model_tibble)
})

setMethod(
  "descriptives", signature(object = "market_model"),
  function(object, variables = NULL) {
    if (is.null(variables)) {
      variables <- object@columns
    }
    variables <- variables[sapply(
      variables,
      function(c) !is.factor(object@model_tibble[[c]])
    )]

    tibble::as_tibble(apply(
      object@model_tibble[, variables], 2,
      function(x) {
        c(
          nobs = length(x), nmval = sum(is.na(x)),
          min = min(x), max = max(x), range = max(x) - min(x),
          sum = sum(x), median = median(x), mean = mean(x),
          mean_se = sqrt(var(x) / length(x)),
          mean_ce = qnorm(0.975) * sqrt(var(x) / length(x)),
          var = var(x), sd = sd(x), coef_var = sd(x) / mean(x)
        )
      }
    ), rownames = "col")
  }
)

#' @rdname market_descriptives
setMethod("demand_descriptives", signature(object = "market_model"), function(object) {
  descriptives(object, object@system@demand@independent_variables)
})

#' @rdname market_descriptives
setMethod("supply_descriptives", signature(object = "market_model"), function(object) {
  descriptives(object, object@system@supply@independent_variables)
})

setGeneric("calculate_initializing_values", function(object) {
  standardGeneric("calculate_initializing_values")
})

setMethod("calculate_initializing_values", signature(object = "market_model"),
          function(object) {
  dlm <- object@system@demand@linear_model

  slm <- object@system@supply@linear_model

  ## Set demand initializing values
  varloc <-
    !(prefixed_independent_variables(object@system@demand) %in% names(dlm$coefficients))
  if (sum(varloc) > 0) {
    print_error(
      object@logger,
      "Misspecified model matrix. ",
      "The matrix should contain all the variables except the variance."
    )
  }
  if (any(is.na(dlm$coefficients))) {
    print_warning(
      object@logger,
      "Setting demand side NA initial values to zero: ",
      paste0(names(dlm$coefficients[is.na(dlm$coefficients)]), collapse = ", "), "."
    )
    dlm$coefficients[is.na(dlm$coefficients)] <- 0
  }
  start_names <- c(
    prefixed_price_variable(object@system@demand),
    prefixed_control_variables(object@system@demand)
  )
  start <- c(dlm$coefficients[start_names])

  ## Set supply initializing values
  varloc <-
    !(prefixed_independent_variables(object@system@supply) %in% names(slm$coefficients))
  if (sum(varloc) > 0) {
    print_error(
      object@logger,
      "Misspecified model matrix. ",
      "The matrix should contain all the variables except the variance."
    )
  }
  if (any(is.na(slm$coefficients))) {
    print_warning(
      object@logger,
      "Setting supply side NA initial values to zero: ",
      paste0(names(slm$coefficients[is.na(slm$coefficients)]), collapse = ", ")
    )
    slm$coefficients[is.na(slm$coefficients)] <- 0
  }
  start_names <- c(
    prefixed_price_variable(object@system@supply),
    prefixed_control_variables(object@system@supply)
  )
  start <- c(start, slm$coefficients[start_names])

  if (object@model_type_string %in% c("Deterministic Adjustment",
                                      "Stochastic Adjustment")) {
    start <- c(start, gamma = 1)
    names(start)[length(start)] <- price_differences_variable(object@system)
  }

  start <- c(start, 1, 1)
  names(start)[(length(start) - 1):length(start)] <- c(
    prefixed_variance_variable(object@system@demand),
    prefixed_variance_variable(object@system@supply)
  )

  if (object@system@correlated_shocks) {
    start <- c(start, rho = 0)
    names(start)[length(start)] <- correlation_variable(object@system)
  }

  start
})

setGeneric("prepare_initializing_values", function(object, initializing_vector) {
  standardGeneric("prepare_initializing_values")
})

setMethod(
  "prepare_initializing_values", signature(object = "market_model"),
  function(object, initializing_vector) {
    if (is.null(initializing_vector)) {
      print_verbose(object@logger, "Initializing using linear regression estimations.")
      initializing_vector <- calculate_initializing_values(object)
    }
    names(initializing_vector) <- likelihood_variables(object@system)
    print_debug(
      object@logger, "Using starting values: ",
      paste(names(initializing_vector), initializing_vector, sep = " = ",
            collapse = ", ")
    )

    initializing_vector
  }
)

#' @title Market side aggregation.
#'
#' @details Calculates the sample's aggregate demand or supply at the passed set of
#' parameters. If the model is static, as is for example the case of
#' \code{\linkS4class{equilibrium_model}}, then all observations are aggregated. If the
#' used data have a time dimension and aggregation per date is required, it can be
#' manually performed using the \code{\link{demanded_quantities}} and
#' \code{\link{supplied_quantities}} functions. If the model has a dynamic component,
#' such as the \code{\linkS4class{diseq_deterministic_adjustment}}, then demanded
#' and supplied quantities are automatically aggregated for each time point.
#' @param object A model object.
#' @param parameters A vector of model's parameters.
#' @return The sum of the estimated demanded or supplied quantities evaluated at the
#' given parameters.
#' @name market_aggregation
#' @examples
#' \donttest{
#' # initialize the basic model using the houses dataset
#' model <- new(
#'   "diseq_basic", # model type
#'   c("ID", "TREND"), "HS", "RM", # keys, quantity, and price variables
#'   "RM + TREND + W + CSHS + L1RM + L2RM + MONTH", # demand specification
#'   "RM + TREND + W + L1RM + MA6DSF + MA3DHF + MONTH", # supply specification
#'   fair_houses(), # data
#'   correlated_shocks = FALSE # allow shocks to be correlated
#' )
#'
#' # estimate the model object (BFGS is used by default)
#' est <- estimate(model)
#'
#' # get estimated aggregate demand
#' aggregate_demand(model, est@coef)
#'
#' # simulate the deterministic adjustment model
#' model <- simulate_model(
#'   "diseq_deterministic_adjustment", list(
#'     # observed entities, observed time points
#'     nobs = 500, tobs = 3,
#'     # demand coefficients
#'     alpha_d = -0.6, beta_d0 = 9.8, beta_d = c(0.3, -0.2), eta_d = c(0.6, -0.1),
#'     # supply coefficients
#'     alpha_s = 0.2, beta_s0 = 4.1, beta_s = c(0.9), eta_s = c(-0.5, 0.2),
#'     # price equation coefficients
#'     gamma = 0.9
#'   ), seed = 1356
#' )
#'
#' # estimate the model object
#' est <- estimate(model)
#'
#' # get estimated aggregate demand
#' aggregate_demand(model, est@coef)
#'
#' # get estimated aggregate demand
#' aggregate_supply(model, est@coef)
#' }
#' @seealso demanded_quantities, supplied_quantities
NULL

setGeneric("aggregate_equation", function(object, parameters, equation) {
  standardGeneric("aggregate_equation")
})

setMethod("aggregate_equation", signature(object = "market_model"),
          function(object, parameters, equation) {
  object@system <- set_parameters(object@system, parameters)
  quantities <- quantities(slot(object@system, equation))
  result <- NULL
  if (!is.null(object@time_column)) {
    time_symbol <- rlang::sym(object@time_column)
    aggregate_symbol <- rlang::sym(colnames(quantities))
    result <- object@model_tibble[, object@time_column] %>%
      dplyr::mutate(!!aggregate_symbol := quantities) %>%
      dplyr::group_by(!!time_symbol) %>%
      dplyr::summarise(!!aggregate_symbol := sum(!!aggregate_symbol))
  } else {
    result <- sum(quantities)
  }
  result
})


#' @describeIn market_aggregation Demand aggregation.
#' @export
setGeneric("aggregate_demand", function(object, parameters) {
  standardGeneric("aggregate_demand")
})

#' @rdname market_aggregation
setMethod("aggregate_demand", signature(object = "market_model"),
          function(object, parameters) {
  aggregate_equation(object, parameters, "demand")
})

#' @title Estimated market quantities.
#'
#' @details Calculates and returns the estimated demanded or supplied quantities for
#' each observation at the passed vector of parameters.
#' @param object A model object.
#' @param parameters A vector of model's parameters.
#' @return A vector with the demanded quantities evaluated at the given parameter
#' vector.
#' @examples
#' \donttest{
#' # initialize the model using the houses dataset
#' model <- new(
#'   "diseq_basic", # model type
#'   c("ID", "TREND"), "HS", "RM", # keys, quantity, and price variables
#'   "RM + TREND + W + CSHS + L1RM + L2RM + MONTH", # demand specification
#'   "RM + TREND + W + L1RM + MA6DSF + MA3DHF + MONTH", # supply specification
#'   fair_houses(), # data
#'   correlated_shocks = FALSE # allow shocks to be correlated
#' )
#'
#' # estimate the model object (BFGS is used by default)
#' est <- estimate(model)
#'
#' # get estimated demanded and supplied quantities
#' head(cbind(demanded_quantities(model, est@coef),
#'      supplied_quantities(model, est@coef)))
#' }
#' @name market_quantities
NULL

#' @describeIn market_quantities Estimated demanded quantities.
#' @export
setGeneric("demanded_quantities", function(object, parameters) {
  standardGeneric("demanded_quantities")
})

#' @rdname market_quantities
setMethod(
  "demanded_quantities", signature(object = "market_model"),
  function(object, parameters) {
    object@system <- set_parameters(object@system, parameters)
    quantities(object@system@demand)
  }
)

#' @describeIn market_aggregation Supply aggregation.
#' @export
setGeneric("aggregate_supply", function(object, parameters) {
  standardGeneric("aggregate_supply")
})

#' @rdname market_aggregation
setMethod("aggregate_supply", signature(object = "market_model"),
          function(object, parameters) {
  aggregate_equation(object, parameters, "supply")
})

#' @describeIn market_quantities Estimated supplied quantities.
#' @export
setGeneric("supplied_quantities", function(object, parameters) {
  standardGeneric("supplied_quantities")
})

#' @rdname market_quantities
setMethod(
  "supplied_quantities", signature(object = "market_model"),
  function(object, parameters) {
    object@system <- set_parameters(object@system, parameters)
    quantities(object@system@supply)
  }
)
