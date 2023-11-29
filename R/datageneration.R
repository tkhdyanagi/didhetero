#' Data generating function
#'
#' @param n The number of cross-sectional units
#' @param tau The length of time series such that tau > 2
#' @param continuous Boolean for whether or not to consider a continuous covariate.
#' If it is TRUE (resp. FALSE), a continuous (resp. discrete) covariate is generated.
#' Default is TRUE.
#'
#' @returns A data.frame that contains the following elements.
#' \item{id}{The unit index}
#' \item{period}{The time period index}
#' \item{Y}{The outcome}
#' \item{G}{The group}
#' \item{Z}{The scalar continuous covariate}
#'
#' @importFrom purrr %>%
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' data <- datageneration(n = 1000, tau = 4, continuous = TRUE)
#'
datageneration <- function(n, tau, continuous = TRUE) {

  # Pre-treatment covariate
  if (continuous) {

    Z <- stats::rnorm(n)

  } else {

    Z <- sample(x = c(-1, 0, 1),
                size = n,
                replace = TRUE,
                prob = c(1/3, 1/3, 1/3))

  }

  # Support of "groups"
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

  # Data generation
  data <- NULL
  for (t in 1:tau) {

    # Untreated potential outcome
    delta_t <- t
    beta_t_0 <- t
    u <- stats::rnorm(n)
    Y_0 <- delta_t + eta + Z * beta_t_0 + u

    # Potential outcome given G = g
    beta_t_g <- beta_t_0 + (G / tau)
    delta_e <- t - G + 1
    v <- stats::rnorm(n)
    Y_g <- delta_t + eta + Z * beta_t_g + delta_e + v

    # Observed outcome
    Y <- (G != 0 & G <= t) * Y_g + (G == 0 | G > t) * Y_0

    # Observed data
    data <- rbind(data,
                  cbind(1:n, t, Y, G, Z))
  }

  # Data arrangement
  colnames(data) <- c("id", "period", "Y", "G", "Z")
  data <- as.data.frame(data) %>%
    dplyr::arrange(id, period)

  # Return
  return(data)

}
