#' First-stage parametric estimation
#'
#' @param xformla
#' A formula for the covariates to include in the model.
#' It should be of the form `~ X1 + X2`.
#' @param data The name of the data.frame that contains the data
#' @param gteval The matrix of the evaluation points g and t
#' @param control_group Which units to use the control group.
#' Options are "nevertreated" and "notyettreated".
#' @param anticipation The number of time periods before participating in the
#' treatment where units can anticipate participating in the treatment and
#' therefore it can affect their untreated potential outcomes.
#' @param period1
#' The first time period
#'
#' @return A list that contains the following elements:
#' GPS: A data.frame of unit index 'id', group index 'g', time period 't', and estimated GPS 'est'.
#' OR: A data.frame of unit index 'id', group index 'g', time period 't', and estimated OR 'est'.
#' GPS_coef: A data.frame of group index 'g', time period 't', and coefficients' estimates 'coef1', 'coef2' and so on.
#' OR_coef: A data.frame of group index 'g', time period 't', and coefficients' estimates 'coef1', 'coef2' and so on.
#'
#' @importFrom dplyr filter mutate pull select
#' @importFrom purrr %>%
#'
#' @noRd
#'
parametric_func <- function(xformla,
                            data,
                            gteval,
                            control_group,
                            anticipation,
                            period1) {

  # Variable definitions -------------------------------------------------------

  # Number of evaluation points (g, t) and vector of evaluation points g's
  if (is.vector(gteval)) {

    num_gteval <- 1

    geval <- gteval[1]

  } else if (is.matrix(gteval)) {

    num_gteval <- nrow(gteval)

    geval <- gteval[, 1] %>%
      unique() %>%
      sort()

  }

  # Number of individuals
  n <- length(unique(data$id))

  # Estimation of generalized propensity score ---------------------------------

  GPS <- GPS_coef <- NULL

  if (control_group == "nevertreated") {

    for (g in geval) {

      # Subset of data with G == g | G == 0
      data1 <- data %>%
        filter(period == period1) %>%
        filter(G == g | G == 0) %>%
        mutate(G_g = ifelse(G == g, 1, 0))

      # Indicator 1{G = g}
      G_g <- data1$G_g

      # Covariates
      covariates <- stats::model.matrix(xformla, data = data1)

      # Dimension of covariates
      dim_X <- ncol(covariates)

      # Logit estimation
      logit <- stats::glm(G_g ~ covariates - 1,
                          family = stats::binomial(link = logit))

      # Estimated coefficients in logit estimation
      pi_coef <- logit$coefficients

      GPS_coef <- rbind(GPS_coef, c(g, pi_coef))

      # Estimated GPS for all units
      data2 <- data %>%
        filter(period == period1)

      covariates <- stats::model.matrix(xformla, data = data2)

      p <- exp(covariates %*% pi_coef)

      GPS <- rbind(GPS,
                   cbind(data2$id, g, p / (1 + p)))

    }

  } else if (control_group == "notyettreated") {

    for (id_gt in 1:num_gteval) {

      # Evaluation points (g, t)
      if (is.vector(gteval)) {

        g <- gteval[1]
        t <- gteval[2]

      } else if (is.matrix(gteval)) {

        g <- gteval[id_gt, 1]
        t <- gteval[id_gt, 2]

      }

      # Subset of data with G == g | G == 0 | G > max(g, t) + delta
      data1 <- data %>%
        filter(period == period1) %>%
        filter(G == g | G == 0 | G > max(g, t) + anticipation) %>%
        mutate(G_g = ifelse(G == g, 1, 0))

      # Indicator 1{G = g}
      G_g <- data1$G_g

      # Covariates
      covariates <- stats::model.matrix(xformla, data = data1)

      # Dimension of covariates
      dim_X <- ncol(covariates)

      # Logit estimation
      logit <- stats::glm(G_g ~ covariates - 1,
                          family = stats::binomial(link = logit))

      # Estimated coefficients in logit estimation
      pi_coef <- logit$coefficients

      GPS_coef <- rbind(GPS_coef, c(g, t, pi_coef))

      # Estimated GPS
      data2 <- data %>%
        filter(period == period1)

      covariates <- stats::model.matrix(xformla, data = data2)

      p <- exp(covariates %*% pi_coef)

      GPS <- rbind(GPS, cbind(data2$id, g, t, p / (1 + p)))

    }
  }

  # Estimation of outcome regression function ----------------------------------

  OR <- OR_coef <- NULL

  if (control_group == "nevertreated") {

    data3 <- data %>%
      filter(G == 0)

    for (id_gt in 1:num_gteval) {

      # Evaluation points (g, t)
      if (is.vector(gteval)) {

        g <- gteval[1]
        t <- gteval[2]

      } else if (is.matrix(gteval)) {

        g <- gteval[id_gt, 1]
        t <- gteval[id_gt, 2]

      }

      # Variables for outcome regression
      Y_t <- data3 %>%
        filter(period == t) %>%
        pull(Y)

      Y_g <- data3 %>%
        filter(period == g - anticipation - 1) %>%
        pull(Y)

      # Covariates
      data3c <- data3 %>%
        filter(period == period1)

      covariates <- stats::model.matrix(xformla, data = data3c)

      # Dimension of covariates
      dim_X <- ncol(covariates)

      # Outcome regression
      OLS <- stats::lm((Y_t - Y_g) ~ covariates - 1)

      beta <- OLS$coefficients

      OR_coef <- rbind(OR_coef, c(g, t, beta))

      # Estimated OR
      data4 <- data %>%
        filter(period == period1)

      covariates <- stats::model.matrix(xformla, data = data4)

      OR <- rbind(OR,
                  cbind(data4$id, g, t, covariates %*% beta))

    }

  } else if (control_group == "notyettreated") {

    for (id_gt in 1:num_gteval) {

      # Evaluation points (g, t)
      if (is.vector(gteval)) {

        g <- gteval[1]
        t <- gteval[2]

      } else if (is.matrix(gteval)) {

        g <- gteval[id_gt, 1]
        t <- gteval[id_gt, 2]

      }

      data3 <- data %>%
        filter(G == 0 | G > max(g, t) + anticipation)

      # Variables for outcome regression
      Y_t <- data3 %>%
        filter(period == t) %>%
        pull(Y)

      Y_g <- data3 %>%
        filter(period == g - anticipation - 1) %>%
        pull(Y)

      # Covariates
      data3c <- data3 %>%
        filter(period == period1)

      covariates <- stats::model.matrix(xformla, data = data3c)

      # Dimension of covariates
      dim_X <- ncol(covariates)

      # Outcome regression
      OLS <- stats::lm((Y_t - Y_g) ~ covariates - 1)

      beta <- OLS$coefficients

      OR_coef <- rbind(OR_coef, c(g, t, beta))

      # Estimated OR
      data4 <- data %>%
        filter(period == period1)

      covariates <- stats::model.matrix(xformla, data = data4)

      OR <- rbind(OR,
                  cbind(data4$id, g, t, covariates %*% beta))

    }
  }

  # Names
  if (control_group == "nevertreated") {

    colnames(GPS)         <- c("id", "g",      "est")
    colnames(OR)          <- c("id", "g", "t", "est")

    if (num_gteval == 1) {

      names(GPS_coef) <- c("g",      paste0("coef", 1:dim_X))
      names(OR_coef)  <- c("g", "t", paste0("coef", 1:dim_X))

    } else {

      colnames(GPS_coef) <- c("g",      paste0("coef", 1:dim_X))
      colnames(OR_coef)  <- c("g", "t", paste0("coef", 1:dim_X))

    }

  } else if (control_group == "notyettreated") {

    colnames(GPS) <- colnames(OR) <- c("id", "g", "t", "est")

    if (num_gteval == 1) {

      names(GPS_coef) <- c("g", "t", paste0("coef", 1:dim_X))
      names(OR_coef)  <- c("g", "t", paste0("coef", 1:dim_X))

    } else {

      colnames(GPS_coef) <- c("g", "t", paste0("coef", 1:dim_X))
      colnames(OR_coef)  <- c("g", "t", paste0("coef", 1:dim_X))

    }
  }

  # Convert to data.frame
  GPS         <- as.data.frame(GPS)
  OR          <- as.data.frame(OR)
  GPS_coef    <- as.data.frame(GPS_coef)
  OR_coef     <- as.data.frame(OR_coef)

  # Return ---------------------------------------------------------------------

  return(list(GPS = GPS,
              GPS_coef = GPS_coef,
              OR = OR,
              OR_coef = OR_coef))

}
