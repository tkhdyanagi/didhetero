#' Data generating function
#'
#' @param n The number of cross-sectional units
#' @param tau The length of time series such that tau > 2
#' @param HC Boolean for whether or not to consider heteroscedastic error terms.
#' If it is TRUE (resp. FALSE), heteroscedastic (resp. homoscedastic) error terms are generated.
#' Default is FALSE.
#' @param dimX The dimension of covariates X.
#' It should be a positive integer.
#' Default is 1.
#' @param DGPY The data generating process for the treated potential outcome equations.
#' Options are 1 and 2.
#' Default is 1.
#' @param continuous Boolean for whether or not to consider a continuous covariate Z.
#' If it is TRUE (resp. FALSE), a continuous (resp. discrete) covariate Z is generated.
#' Default is TRUE.
#'
#' @returns A data.frame that contains the following elements.
#' \item{id}{The unit index}
#' \item{period}{The time period index}
#' \item{Y}{The outcome}
#' \item{G}{The group}
#' \item{Z}{The scalar continuous covariate}
#' \item{X}{The other covariates}
#'
#' @importFrom purrr %>%
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' data <- datageneration(n = 500,
#'                        tau = 4,
#'                        HC = FALSE,
#'                        dimX = 1,
#'                        DGPY = 1,
#'                        continuous = TRUE)
#'
datageneration <- function(n, tau, HC = FALSE, dimX = 1, DGPY = 1, continuous = TRUE) {

  # Pre-treatment covariate Z
  if (continuous) {

    Z <- stats::rnorm(n)

  } else {

    Z <- sample(x = c(-1, 0, 1),
                size = n,
                replace = TRUE,
                prob = c(1/3, 1/3, 1/3))

  }

  # All pre-treatment covariates X
  X <- Z

  if (dimX > 1) {
    for (k in 1:(dimX - 1)) {

      X <- cbind(X, stats::rnorm(n))

    }
  }

  # Support of "groups" where 0 indicates the never treated group
  group_supp <- c(0, 2:tau)

  # The n * tau matrix of group choice probability
  group_P <- NULL
  for (g in group_supp) {
    gamma <- 0.5 * g / tau
    group_P <- cbind(group_P, exp(Z * gamma))
  }
  group_P <- group_P / rowSums(group_P)

  # Group choice
  G <- rep(NA, n)
  for (i in 1:n) {
    G[i] <- sample(x = group_supp, size = 1, prob = group_P[i, ])
  }

  # Individual effect
  eta <- stats::rnorm(n, mean = G, sd = 1)

  # Heteroscedasticity or homoscedasticity for error terms
  if (HC) {
    u_sd <- 0.5 + stats::pnorm(Z)
    v_sd <- (G / tau) + stats::pnorm(Z)
  } else {
    u_sd <- 1
    v_sd <- 1
  }

  # Data generation
  data <- NULL

  for (t in 1:tau) {

    # Error term for the untreated potential outcome
    u <- stats::rnorm(n, mean = 0, sd = u_sd)

    # Untreated potential outcome
    delta_t <- t
    beta_t_0 <- rep(t, dimX) / seq(1, dimX, by = 1)
    if (dimX == 1) {

      Y_0 <- delta_t + eta + X * beta_t_0 + u

    } else if (dimX > 1) {

      Y_0 <- delta_t + eta + X %*% beta_t_0 + u

    }

    # Error term for the treated potential outcome
    v <- stats::rnorm(n, mean = 0, sd = v_sd)

    # Coefficient in the treated potential outcome
    delta_e <- t - G + 1

    # A term in the potential outcome equation
    if (DGPY == 1) {

      mgt <- (G / t) * sin(pi * Z)

    } else if (DGPY == 2) {

      mgt <- Z * G / t

    }

    # Potential outcome given G = g
    Y_g <- Y_0 + mgt + delta_e + v - u

    # Observed outcome
    Y <- (G != 0 & G <= t) * Y_g + (G == 0 | G > t) * Y_0

    # Observed data
    data <- rbind(data,
                  cbind(1:n, t, Y, G, X))
  }

  # Column names
  name <- c("id", "period", "Y", "G", "Z")

  if (dimX > 1) {
    for (k in 1:(dimX - 1)) {

      name <- c(name, paste0("X", k))

    }
  }

  # Data arrangement
  colnames(data) <- name
  data <- as.data.frame(data) %>%
    dplyr::arrange(id, period)

  # Return
  return(data)

}
