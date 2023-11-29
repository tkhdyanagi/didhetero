#' Kernel density estimation
#'
#' @param x An n-dimensional vector
#' @param eval A vector of evaluation points
#' @param kernel The kernel function. Options are `epa` and `gau`.
#' @param h A vector or scalar of bandwidth.
#' If h is a vector, the length of h should equal to that of eval.
#'
#' @return A vector of kernel density estimates
#'
#' @noRd
#'
kde_func <- function(x, eval, kernel, h) {

  # Sample size
  n <- length(x)

  # Number of evaluation points
  num_eval <- length(eval)

  # Number of bandwidths
  num_h <- length(h)

  if (num_h == 1) {
    h <- rep(h, num_eval)
  }

  # A vector to save kernel density estimates
  Estimate <- NULL

  # Kernel density estimation
  for (j in 1:num_eval) {

    Estimate <-
      c(Estimate,
        sum(kernel_func(u = (x - eval[j])/h[j], kernel = kernel)) / (n * h[j]))

  }

  # Return
  return(Estimate)

}

#' Local polynomial regression estimation
#'
#' @param y An n-dimensional vector of dependent variables
#' @param x An n-dimensional vector of regressors
#' @param eval A vector of evaluation points
#' @param p A polynomial order
#' @param deriv A derivative order to be estimated
#' @param kernel The kernel function. Options are `epa` and `gau`.
#' @param h A vector or scalar of bandwidth.
#' If h is a vector, the length of h should equal to that of eval.
#' @param weight
#' A scalar weight or an n-dimensional vector of multiplier bootstrap weights.
#' Default is 1, which corresponds to the original estimator without bootstrapping.
#'
#' @return A vector of local polynomial regression estimates
#'
#' @importFrom purrr %>%
#'
#' @noRd
#'
lpr_func <- function(y, x, eval, p, deriv, kernel, h, weight = 1) {

  # Number of evaluation points
  num_eval <- length(eval)

  # Number of bandwidths
  num_h <- length(h)

  if (num_h == 1) {
    h <- rep(h, num_eval)
  }

  # A vector to save LPR estimates
  Estimate <- NULL

  # Local polynomial regression estimation
  for (j in 1:num_eval) {

    # An n times (p + 1) matrix of polynomials
    Z_mat <- polynomial_func(u = x - eval[j], p = p)

    # An n-dimensional vector of (kernel values * weights)
    K_h <- weight * kernel_func(u = (x - eval[j])/h[j], kernel = kernel)

    # Observations with non-zero kernel values
    index   <- (K_h != 0)
    Z_mat_1 <- Z_mat[index, ]
    K_h_1   <- K_h[index]
    y_1     <- y[index]

    if (sum(index) == 0) {

      Estimate <- c(Estimate, NA)

    } else if (sum(index) == 1) {

      Estimate <- c(Estimate, y_1)

    } else {

      # Inverse matrix
      Gamma_inv <- crossprod(K_h_1 * Z_mat_1, Z_mat_1) %>%
        chol() %>%
        chol2inv()

      # LPR coefficient
      beta <- Gamma_inv %*% crossprod(K_h_1 * Z_mat_1, y_1)

      # Unit vector
      e_vec <- numeric(p + 1)
      e_vec[deriv + 1] <- 1

      # LPR estimate
      Estimate <- c(Estimate,
                    gamma(deriv + 1) * crossprod(e_vec, beta))

    }
  }

  # Return
  return(Estimate)

}

#' Compute a vector of kernel values
#'
#' @param u An n-dimensional vector
#' @param kernel The kernel function. Options are `epa` and `gau`.
#'
#' @returns An n-dimensional vector of kernel values
#'
#' @noRd
#'
kernel_func <- function(u, kernel) {
  if (kernel == "epa") {

    return((abs(u) <= 1) * 0.75 * (1 - u^2))

  } else if (kernel == "gau") {

    return(stats::dnorm(u))

  }
}

#' Compute a matrix of polynomials
#'
#' @param u An n-dimensional vector
#' @param p A polynomial order
#'
#' @return An n times (p + 1) matrix of polynomials
#'
#' @noRd
#'
polynomial_func <- function(u, p) {
  r_p <- NULL
  for (j in 0:p) {
    if (length(u) == 1) {
      r_p <- c(r_p, u^j)
    } else if (length(u) > 1) {
      r_p <- cbind(r_p, u^j)
    }
  }
  return(r_p)
}

#' Compute some integrals related to the kernel function
#'
#' @param l A non-negative integer
#' @param m A non-negative integer
#' @param kernel The kernel function. Options are `epa` and `gau`.
#'
#' @returns The value of integral
#'
#' @noRd
#'
integrate_func <- function(l, m, kernel) {

  # Kernel function
  K <- function(u) {
    kernel_func(u, kernel = kernel)
  }

  # Compute the integral
  integral <- stats::integrate(f = function(u){u^l * (K(u))^m},
                               lower = -Inf,
                               upper =  Inf)

  # Return
  return(integral)

}
