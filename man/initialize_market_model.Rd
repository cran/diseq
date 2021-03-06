% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/market_model.R, R/diseq_basic.R,
%   R/diseq_deterministic_adjustment.R, R/diseq_directional.R,
%   R/diseq_stochastic_adjustment.R, R/equilibrium_model.R
\name{initialize_market_model}
\alias{initialize_market_model}
\alias{initialize,diseq_basic-method}
\alias{initialize,diseq_deterministic_adjustment-method}
\alias{initialize,diseq_directional-method}
\alias{initialize,diseq_stochastic_adjustment-method}
\alias{initialize,equilibrium_model-method}
\title{Model initialization}
\usage{
\S4method{initialize}{diseq_basic}(
  .Object,
  key_columns,
  quantity_column,
  price_column,
  demand_specification,
  supply_specification,
  data,
  correlated_shocks = TRUE,
  verbose = 0
)

\S4method{initialize}{diseq_deterministic_adjustment}(
  .Object,
  key_columns,
  time_column,
  quantity_column,
  price_column,
  demand_specification,
  supply_specification,
  data,
  correlated_shocks = TRUE,
  verbose = 0
)

\S4method{initialize}{diseq_directional}(
  .Object,
  key_columns,
  time_column,
  quantity_column,
  price_column,
  demand_specification,
  supply_specification,
  data,
  correlated_shocks = TRUE,
  verbose = 0
)

\S4method{initialize}{diseq_stochastic_adjustment}(
  .Object,
  key_columns,
  time_column,
  quantity_column,
  price_column,
  demand_specification,
  supply_specification,
  price_specification,
  data,
  correlated_shocks = TRUE,
  verbose = 0
)

\S4method{initialize}{equilibrium_model}(
  .Object,
  key_columns,
  quantity_column,
  price_column,
  demand_specification,
  supply_specification,
  data,
  correlated_shocks = TRUE,
  verbose = 0
)
}
\arguments{
\item{.Object}{The object to be Constructed.}

\item{key_columns}{Key columns of the data set.}

\item{quantity_column}{The quantity variable of the data set.}

\item{price_column}{The price variable of the data set.}

\item{demand_specification}{A formula representation of the right hand side of the
demand equation.}

\item{supply_specification}{A formula representation of the right hand side of the
supply equation.}

\item{data}{The data set.}

\item{correlated_shocks}{Should the model be estimated using correlated shocks?}

\item{verbose}{Verbosity level.}

\item{time_column}{The time column of the data set.}

\item{price_specification}{A formula representation of the price equation.}
}
\value{
The initialized model.
}
\description{
Model initialization
}
\details{
The following two subsections describe the common initialization steps of all market
model classes.

\subsection{Variable construction}{
The constructor prepares the model's variables using the passed specifications. The
specification strings are expected to follow the syntax of
\code{\link[stats]{formula}}. The construction of the model's data uses the variables
that are extracted by these specification. The demand variables are extracted by a
formula that uses the \code{quantity_column} on the left hand side and the
\code{demand_specification} on the right hand side of the formula. The supply
variables are constructed by the the\code{quantity_column} and the
\code{supply_specification}. In the case of the
\code{\linkS4class{diseq_stochastic_adjustment}} model, the price dynamics'
variables are extracted using the \code{price_specification}. The
\code{price_specification} for the \code{\linkS4class{diseq_stochastic_adjustment}}
should contain only terms other than that of excess demand. The excess demand term of
the price equation is automatically added by the constructor.
}

\subsection{Data preparation}{
  1. If the passed data set contains rows with NA values, they are dropped. If the
verbosity level allows warnings, a warning is emitted reporting how many rows were
dropped.

  2. After dropping the rows, factor levels may be invalidated. If needed, the
constructor readjusts the factor variables by removing the unobserved levels. Factor
indicators and interaction terms are automatically created.

  3. The primary key column is constructed by pasting the values of the key_columns.

  4. In the cases of the \code{\linkS4class{diseq_directional}},
\code{\linkS4class{diseq_deterministic_adjustment}}, and
the \code{\linkS4class{diseq_stochastic_adjustment}} models, a column with lagged
prices is constructed. Since lagged prices are unavailable for the observations of
the first time point, these observations are dropped. If the verbosity level allows
the emission of information messages, the constructor prints the number of dropped
observations.

  5. In the cases of the \code{\linkS4class{diseq_directional}}
and the \code{\linkS4class{diseq_stochastic_adjustment}} models, a column with price
differences is created.
}
}
\section{Functions}{
\itemize{
\item \code{initialize,diseq_basic-method}: Basic disequilibrium model base constructor

\item \code{initialize,diseq_deterministic_adjustment-method}: Disequilibrium model with deterministic price
adjustment constructor

\item \code{initialize,diseq_directional-method}: Directional disequilibrium model base constructor

\item \code{initialize,diseq_stochastic_adjustment-method}: Disequilibrium model with stochastic price
adjustment constructor

\item \code{initialize,equilibrium_model-method}: Equilibrium model constructor
}}

\examples{
simulated_data <- simulate_data(
  "diseq_basic", 500, 3, # model type, observed entities, observed time points
  -0.9, 8.9, c(0.3, -0.2), c(-0.03, -0.01), # demand coefficients
  0.9, 4.2, c(0.03), c(-0.05, 0.02) # supply coefficients
)

# initialize the model
model <- new(
  "diseq_basic", # model type
  c("id", "date"), "Q", "P", # keys, quantity, and price variables
  "P + Xd1 + Xd2 + X1 + X2", "P + Xs1 + X1 + X2", # equation specifications
  simulated_data, # data
  correlated_shocks = FALSE # use independent shocks
)

show(model)
simulated_data <- simulate_data(
  # model type, observed entities and time points
  "diseq_deterministic_adjustment", 500, 3,
  # demand coefficients
  -0.9, 8.9, c(0.03, -0.02), c(-0.03, -0.01),
  # supply coefficients
  0.9, 4.2, c(0.03), c(0.05, 0.02),
  # price adjustment coefficient
  1.4
)

# initialize the model
model <- new(
  "diseq_deterministic_adjustment", # model type
  c("id", "date"), "date", "Q", "P", # keys, quantity, and price variables
  "P + Xd1 + Xd2 + X1 + X2", "P + Xs1 + X1 + X2", # equation specifications
  simulated_data, # data
  correlated_shocks = TRUE # allow shocks to be correlated
)

show(model)
\donttest{
simulated_data <- simulate_data(
  "diseq_directional", 500, 3, # model type, observed entities, observed time points
  -0.2, 4.3, c(0.03, 0.02), c(0.03, 0.01), # demand coefficients
  0.0, 4.0, c(0.03), c(0.05, 0.02) # supply coefficients
)

# in the directional model prices cannot be included in both demand and supply
model <- new(
  "diseq_directional", # model type
  c("id", "date"), "date", "Q", "P", # keys, time point, quantity, and price variables
  "P + Xd1 + Xd2 + X1 + X2", "Xs1 + X1 + X2", # equation specifications
  simulated_data, # data
  correlated_shocks = TRUE # allow shocks to be correlated
)

show(model)
}
simulated_data <- simulate_data(
  # model type, observed entities and time points
  "diseq_stochastic_adjustment", 500, 3,
  # demand coefficients
  -0.1, 9.8, c(0.3, -0.2), c(0.6, 0.1),
  # supply coefficients
  0.1, 5.1, c(0.9), c(-0.5, 0.2),
  # price adjustment coefficient
  1.4, 3.1, c(0.8)
)

# initialize the model
model <- new(
  "diseq_stochastic_adjustment", # model type
  c("id", "date"), "date", "Q", "P", # keys, quantity, and price variables
  "P + Xd1 + Xd2 + X1 + X2", "P + Xs1 + X1 + X2", "Xp1", # equation specifications
  simulated_data, # data
  correlated_shocks = TRUE # allow shocks to be correlated
)

show(model)
simulated_data <- simulate_data(
  "equilibrium_model", 500, 3, # model type, observed entities and time points
  -0.9, 14.9, c(0.3, -0.2), c(-0.03, -0.01), # demand coefficients
  0.9, 3.2, c(0.3), c(0.5, 0.02) # supply coefficients
)

# initialize the model
model <- new(
  "equilibrium_model", # model type
  c("id", "date"), "Q", "P", # keys, quantity, and price variables
  "P + Xd1 + Xd2 + X1 + X2", "P + Xs1 + X1 + X2", # equation specifications
  simulated_data, # data
  correlated_shocks = TRUE # allow shocks to be correlated
)

show(model)
}
