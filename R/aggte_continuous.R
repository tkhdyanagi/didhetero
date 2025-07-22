#' Uniform Inference for Summary Parameters that aggregate CATTs
#'
#' `aggte_continuous` computes doubly robust estimates and uniform confidence
#' bands for summary parameters that aggregate group-time conditional
#' average treatment effects (CATT) given a continuous pre-treatment covariate
#' in the staggered difference-in-differences setup of
#' Callaway and Sant'Anna (2021) for balanced panel data.
#' See Imai, Qin, and Yanagi (2025) for more details.
#'
#' @param output The output of the `catt_gt_continuous` function.
#' In doing so, several arguments and the uniform inference results for CATT from
#' the `catt_gt_continuous` function can be used in this function.
#' In particular, the following are passed down from output:
#' xformla, zeval, pretrend, control_group, anticipation, alp, and kernel.
#' @param type Which type of the summary parameter is of interest.
#' Options are "simple", "dynamic", "group", and "calendar".
#' Default is "dynamic".
#' @param eval The vector of the evaluation point specific to the chosen summary parameter.
#' If type is set to "dynamic", it is the evaluation point e.
#' If type is set to "group", it is the evaluation point g'.
#' If type is set to "calendar", it is the evaluation point t'.
#' If type is set to "simple", there is no evaluation point specific to this summary parameter, and eval should be `NULL`.
#' Default is `NULL` and eval is constructed automatically.
#' @param bstrap
#' Boolean for whether or not to perform multiplier bootstrapping.
#' Default is `TRUE`.
#' If bstrap is `FALSE`, only the analytical critical value is used.
#' @param biters
#' The number of bootstrap iterations.
#' This parameter is only applicable if `bstrap = TRUE`.
#' Default is 1000.
#' @param porder
#' The polynomial order used for the second- and third-stage estimation.
#' Options are 1 and 2,
#' which correspond to the local linear and quadratic regressions, respectively.
#' Default is 2.
#' @param bwselect
#' The bandwidth selection method used for the aggregation.
#' Options are "IMSE1", "IMSE2", "US1", and "manual".
#' "IMSE1" and "IMSE2" mean the IMSE-optimal bandwidths for the local linear and quadratic regressions, respectively.
#' "US1" means the rule-of-thumb undersmoothing for the local linear regression.
#' "manual" means the manual selection and bw should be specified in this case.
#' Default is "IMSE1", which is recommended for use with `porder = 2`.
#' @param bw
#' The bandwidth used for the aggregation.
#' Default is `NULL` and the bandwidth is chosen automatically.
#' This parameter is only applicable if bwselect is "manual", and
#' should be a scalar or a vector whose length equals to the number of rows of eval.
#' @param uniformall
#' Boolean for whether or not to perform the uniform inference over both eval and z.
#' Default is `TRUE` and the uniform inference over eval and z is performed.
#' If `FALSE`, the uniform inference only over z is performed.
#'
#' @return A list that contains the following elements.
#' \item{Estimate}{A data.frame that contains the following elements. \cr
#' eval: The evaluation point specific to the chosen summary parameter. \cr
#' z: The covariate value. \cr
#' est: The doubly robust estimate of the chosen summary parameter. \cr
#' se: The standard error. \cr
#' ci1_lower: The lower bound of the analytical UCB. \cr
#' ci1_upper: The upper bound of the analytical UCB. \cr
#' ci2_lower: The lower bound of the bootstrap UCB. \cr
#' ci2_upper: The upper bound of the bootstrap UCB. \cr
#' bw: The bandwidth.}
#' \item{Figure1}{A list that contains the ggplot elements for the analytical UCB}
#' \item{Figure2}{A list that contains the ggplot elements for the bootstrap UCB}
#'
#' @importFrom dplyr arrange filter mutate pull sym
#' @importFrom purrr %>%
#' @importFrom rlang UQ
#' @importFrom utils globalVariables
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' data <- datageneration(
#'   n = 500,
#'   tau = 4,
#'   continuous = TRUE
#' )
#'
#' output1 <- catt_gt_continuous(
#'   yname = "Y",
#'   tname = "period",
#'   idname = "id",
#'   gname = "G",
#'   zname = "Z",
#'   xformla = ~ Z,
#'   data = data,
#'   zeval = seq(-1, 1, by = 0.1),
#'   gteval = NULL,
#'   pretrend = FALSE,
#'   control_group = "notyettreated",
#'   anticipation = 0,
#'   alp = 0.05,
#'   bstrap = TRUE,
#'   biters = 1000,
#'   porder = 2,
#'   kernel = "gau",
#'   bwselect = "IMSE1",
#'   bw = NULL,
#'   uniformall = TRUE
#' )
#'
#' output2 <- aggte_continuous(
#'   output = output1,
#'   type = "dynamic",
#'   eval = NULL,
#'   bstrap = TRUE,
#'   biters = 1000,
#'   porder = 2,
#'   bwselect = "IMSE1",
#'   bw = NULL,
#'   uniformall = TRUE
#' )
#'
#' }
#'
#' @references
#' Callaway, B., & Sant’Anna, P. H. (2021).
#' Difference-in-differences with multiple time periods.
#' Journal of Econometrics, 225(2), 200-230.
#'
#' Imai, S., Qin, L., & Yanagi, T. (2025).
#' Doubly robust uniform confidence bands
#' for group-time conditional average treatment effects
#' in difference-in-differences.
#' arXiv preprint arXiv:2305.02185.
#'
aggte_continuous <- function(output,
                             type = "dynamic",
                             eval = NULL,
                             bstrap = TRUE,
                             biters = 1000,
                             porder = 2,
                             bwselect = "IMSE1",
                             bw = NULL,
                             uniformall = TRUE) {

  #-----------------------------------------------------------------------------
  # Error handling
  #-----------------------------------------------------------------------------

  error <- errorhandling_aggte_continuous(output = output,
                                          type = type,
                                          eval = eval,
                                          bstrap = bstrap,
                                          biters = biters,
                                          bwselect = bwselect,
                                          bw = bw,
                                          uniformall = uniformall)

  #-----------------------------------------------------------------------------
  # Variable definitions
  #-----------------------------------------------------------------------------

  # Arguments of catt_gt_continuous
  yname         <- output$Arguments$yname
  tname         <- output$Arguments$tname
  idname        <- output$Arguments$idname
  gname         <- output$Arguments$gname
  zname         <- output$Arguments$zname
  xformla       <- output$Arguments$xformla
  data          <- output$Arguments$data
  zeval         <- output$Arguments$zeval
  control_group <- output$Arguments$control_group
  anticipation  <- output$Arguments$anticipation
  alp           <- output$Arguments$alp
  kernel        <- output$Arguments$kernel

  # Other variable definitions
  gteval0    <- output$gteval
  all_B_g_t  <- output$B_g_t
  all_G_g    <- output$G_g
  all_mu_G_g <- output$mu_G_g
  gbar       <- output$gbar
  kd0_Z      <- output$kd0_Z
  kd1_Z      <- output$kd1_Z
  Z          <- output$Z
  Z_supp     <- seq(min(Z), max(Z), length = 100)

  # Exclude the pre-treatment periods for "group", "calendar", and "simple"
  if (type != "dynamic") {

    gteval0 <- gteval0[gteval0[, 2] - gteval0[, 1] >= 0, ]

  }

  # Construct evaluation points (g, t, .) specific to the chosen summary parameter
  # That is:
  # (g, t, e)  for "dynamic"
  # (g, t, g') for "group"
  # (g, t, t') for "calendar"
  # (g, t, NA) for "simple"
  if (is.null(eval)) {
    if (type == "dynamic") {
      eval <- sort(unique(gteval0[, 2] - gteval0[, 1]))
    } else if (type == "group") {
      eval <- sort(unique(gteval0[, 1]))
    } else if (type == "calendar") {
      eval <- sort(unique(gteval0[, 2]))
    }
  }

  if (type == "simple") {
    eval <- NA
  }

  # Number of units
  n <- nrow(all_G_g)

  # Number of evaluation points z
  num_zeval <- length(zeval)

  # Number of evaluation points in eval
  num_eval <- length(eval)

  # Specify I_l_Km = \int u^l K^m(u) du and lambda
  if (kernel == "epa") {

    I_2_K1 <- 0.2
    I_4_K1 <- 0.08571429
    I_6_K1 <- 0.04761905
    I_0_K2 <- 0.6
    I_2_K2 <- 0.08571429
    I_4_K2 <- 0.02857143
    I_6_K2 <- 0.01298701

    lambda <- 2.5

  } else if (kernel == "gau") {

    I_2_K1 <- 1
    I_4_K1 <- 3
    I_6_K1 <- 15
    I_0_K2 <- 0.2820948
    I_2_K2 <- 0.1410474
    I_4_K2 <- 0.2115711
    I_6_K2 <- 0.5289277

    lambda <- 0.5

  }

  # Constants in \mathcal{V}_{\theta} for the LLR and LQR estimators
  const_V1 <- I_0_K2

  const_V2 <-
    (I_4_K1^2 * I_0_K2 - 2 * I_2_K1 * I_4_K1 * I_2_K2 + I_2_K1^2 * I_4_K2) /
    (I_4_K1 - I_2_K1^2)^2

  const_V <- (porder == 1) * const_V1 + (porder == 2) * const_V2

  # Construct a matrix of evaluation points (g, t, .)
  all_gteeval <- NULL

  for (id_gt in 1:nrow(gteval0)) {

    g0 <- gteval0[id_gt, 1]
    t0 <- gteval0[id_gt, 2]

    for (id_eval in 1:num_eval) {

      e0 <- eval[id_eval]

      if (type == "dynamic" & e0 == t0 - g0 & g0 + e0 < gbar) {

        all_gteeval <- rbind(all_gteeval, c(g0, t0, e0))

      } else if (type == "group" & g0 == e0 & t0 >= e0 & t0 < gbar) {

        all_gteeval <- rbind(all_gteeval, c(g0, t0, e0))

      } else if (type == "calendar" & g0 <= e0 & t0 == e0) {

        all_gteeval <- rbind(all_gteeval, c(g0, t0, e0))

      } else if (type == "simple" & g0 <= t0 & t0 < gbar) {

        all_gteeval <- rbind(all_gteeval, c(g0, t0, e0))

      }
    }
  }

  #-----------------------------------------------------------------------------
  # Aggregation with pilot bandwidth
  #-----------------------------------------------------------------------------

  # Empty lists
  list_aggte_est    <- list()
  list_aggte_weight <- list()
  list_B_g_t        <- list()
  list_catt         <- list()
  list_ci1_lower    <- list()
  list_ci1_upper    <- list()
  list_ci2_lower    <- list()
  list_ci2_upper    <- list()
  list_G_g          <- list()
  list_J            <- list()
  list_mu_G_g       <- list()
  list_se           <- list()
  list_xi_g_t       <- list()

  for (id_eval in 1:num_eval) {

    # Last element of (g, t, .)
    e1 <- eval[id_eval]

    if (type %in% c("dynamic", "group", "calendar")) {

      # Number and index of evaluation points (g, t, .) such that last element = this e
      num_gte   <-   sum(all_gteeval[, 3] == e1)
      index_gte <- which(all_gteeval[, 3] == e1)

      # Evaluation points (g, t, .) such that last element = this e
      gteeval <- all_gteeval[index_gte, ]

    } else if (type == "simple") {

      # Same as above but need modifications when "simple" because e = NA
      num_gte   <- nrow(all_gteeval)
      index_gte <- 1:num_gte
      gteeval   <- all_gteeval

    }

    if (num_gte == 1) {

      # Specify g and t
      g1 <- gteeval[1]
      t1 <- gteeval[2]

      # CATT inference results
      output1 <- output$Estimate %>%
        filter(g == g1, t == t1)

      # Save as lists
      list_aggte_est    <- append(list_aggte_est,    list(output1$est))
      list_aggte_weight <- append(list_aggte_weight, list(rep(1, num_zeval)))
      list_B_g_t        <- append(list_B_g_t,        list(all_B_g_t[, , index_gte]))
      list_catt         <- append(list_catt,         list(output1$est))
      list_ci1_lower    <- append(list_ci1_lower,    list(output1$ci1_lower))
      list_ci1_upper    <- append(list_ci1_upper,    list(output1$ci1_upper))
      list_ci2_lower    <- append(list_ci2_lower,    list(output1$ci2_lower))
      list_ci2_upper    <- append(list_ci2_upper,    list(output1$ci2_upper))
      list_G_g          <- append(list_G_g,          list(all_G_g[, index_gte]))
      list_J            <- append(list_J,            list(matrix(NA, nrow = n, ncol = num_zeval)))
      list_mu_G_g       <- append(list_mu_G_g,       list(all_mu_G_g[, index_gte]))
      list_se           <- append(list_se,           list(output1$se))
      list_xi_g_t       <- append(list_xi_g_t,       list(matrix(NA, nrow = n, ncol = num_zeval)))

    } else if (num_gte > 1) {

      # Subsets used for aggregation
      B_g_t  <- all_B_g_t[, , index_gte]
      G_g    <- all_G_g[, index_gte]
      mu_G_g <- all_mu_G_g[, index_gte]

      # Aggregation weights
      aggte_weight <- matrix(NA, nrow = num_zeval, ncol = num_gte)

      # Estimated \xi_g_t
      xi_g_t <- array(NA, dim = c(n, num_zeval, num_gte))

      for (id_gte in 1:num_gte) {

        if (type %in% c("dynamic", "calendar")) {

          # Estimated weights \hat w_{g,t}^{es}(e, z) or \hat w_{c}(t)
          aggte_weight[, id_gte] <- mu_G_g[, id_gte] / rowSums(mu_G_g)

          for (i in 1:n) {

            # Estimated Xi
            xi_g_t[i, , id_gte] <- G_g[i, id_gte] / rowSums(mu_G_g) -
              mu_G_g[, id_gte] / (rowSums(mu_G_g))^2 * sum(G_g[i, ])

          }

        } else if (type == "group") {

          # Weights w_{g,t}^{sel}(g', z)
          aggte_weight[, id_gte] <- 1 / sum(gteeval[, 1] == e1 & e1 <= gteeval[, 2] & gteeval[, 2] < gbar)

          # Estimated Xi
          xi_g_t[, , id_gte] <- 0

        } else if (type == "simple") {

          # Estimated weights \hat w_{g,t}^{OW}(z) * \hat \kappa(z)
          aggte_weight[, id_gte] <- mu_G_g[, id_gte] / rowSums(mu_G_g)

        }
      }

      # Matrix of CATT estimates
      catt <- matrix(NA, nrow = num_zeval, ncol = num_gte)

      # Estimated summands in J
      J_temp <- array(NA, dim = c(n, num_zeval, num_gte))

      for (id_gte in 1:num_gte) {

        # Specify (g, t, .)
        g1 <- gteeval[id_gte, 1]
        t1 <- gteeval[id_gte, 2]
        e1 <- gteeval[id_gte, 3]

        # CATT estimates \hat DR_{g,t}(z)
        catt[, id_gte] <- output$Estimate %>%
          filter(g == g1, t == t1) %>%
          select(est) %>%
          as.matrix()

        if (type == "simple") {

          # Estimated \kappa(z)
          aggte_kappa <- rowSums(aggte_weight)

          for (i in 1:n) {

            # Estimated Xi
            xi_g_t[i, , id_gte] <-
              G_g[i, id_gte] / (aggte_kappa * rowSums(mu_G_g)) -
              mu_G_g[, id_gte] / (aggte_kappa * (rowSums(mu_G_g))^2) * sum(G_g[i, ]) -
              mu_G_g[, id_gte] / (aggte_kappa * rowSums(mu_G_g))^2 * rowSums(G_g[i, ] - mu_G_g / rowSums(mu_G_g) * sum(G_g[i, ]))

          }
        }

        for (i in 1:n) {

          # Summands in the estimated J
          J_temp[i, , id_gte] <- aggte_weight[, id_gte] * B_g_t[i, , id_gte] +
            catt[, id_gte] * xi_g_t[i, , id_gte]

        }
      }

      # Compute the estimated J
      J <- matrix(NA, nrow = n, ncol = num_zeval)

      for (r in 1:num_zeval) {
        J[, r] <- rowSums(J_temp[, r, ])
      }

      # Estimated summary parameter
      aggte_est <- rowSums(catt * aggte_weight)

      # Save as lists
      list_aggte_est    <- append(list_aggte_est,    list(aggte_est))
      list_aggte_weight <- append(list_aggte_weight, list(aggte_weight))
      list_B_g_t        <- append(list_B_g_t,        list(B_g_t))
      list_catt         <- append(list_catt,         list(catt))
      list_ci1_lower    <- append(list_ci1_lower,    list(NA))
      list_ci1_upper    <- append(list_ci1_upper,    list(NA))
      list_ci2_lower    <- append(list_ci2_lower,    list(NA))
      list_ci2_upper    <- append(list_ci2_upper,    list(NA))
      list_G_g          <- append(list_G_g,          list(G_g))
      list_J            <- append(list_J,            list(J))
      list_mu_G_g       <- append(list_mu_G_g,       list(mu_G_g))
      list_se           <- append(list_se,           list(NA))
      list_xi_g_t       <- append(list_xi_g_t,       list(xi_g_t))

    }
  }

  #-----------------------------------------------------------------------------
  # Bandwidth selection for uniform inference on the summary parameter
  #-----------------------------------------------------------------------------

  if (bwselect %in% c("IMSE1", "IMSE2", "US1")) {

    # List of bandwidths
    bw <- NULL

    for (id_eval in 1:num_eval) {

      # Last element of (g, t, .)
      e1 <- eval[id_eval]

      if (type %in% c("dynamic", "group", "calendar")) {

        # Number and index of evaluation points (g, t, .) such that last element = this e
        num_gte   <-   sum(all_gteeval[, 3] == e1)
        index_gte <- which(all_gteeval[, 3] == e1)

        # Evaluation points (g, t, .) such that last element = this e
        gteeval <- all_gteeval[index_gte, ]

      } else if (type == "simple") {

        # Same as above but need modifications when "simple" because e = NA
        num_gte   <- nrow(all_gteeval)
        index_gte <- 1:num_gte
        gteeval   <- all_gteeval

      }

      if (num_gte == 1) {

        # Specify (g, t, .)
        g1 <- gteeval[1]
        t1 <- gteeval[2]
        e1 <- gteeval[3]

        # Bandwidth used for CATT inference
        bw_temp <- output$Estimate %>%
          filter(g == g1, t == t1) %>%
          select(bw) %>%
          as.matrix()

        # Save
        bw <- c(bw, bw_temp[1])

      } else if (num_gte > 1) {

        # Initialization
        mathcal_B <- mathcal_V <- NULL

        # Pull J
        J <- list_J[[id_eval]]

        # LLR residuals used for standard errors and bootstrap inference
        U_hat <- matrix(NA, nrow = n, ncol = num_zeval)

        # Estimate \mathcal{B}_{\theta}(z) for the LLR or LQR estimator --------

        mathcal_B <- rep(NA, num_zeval)

        for (r in 1:num_zeval) {
          if (bwselect %in% c("IMSE1", "US1")) {

            # Bandwidth selection for estimating the second order derivative
            mu_J_2_bw <- nprobust::lpbwselect(y = J[, r],
                                              x = Z,
                                              eval = zeval[r],
                                              p = 3,
                                              deriv = 2,
                                              kernel = kernel,
                                              bwselect = "mse-dpi")$bws[, 2]

            # Estimate the second order derivative
            mu_J_2 <- lpr_func(y = J[, r],
                               x = Z,
                               eval = zeval[r],
                               p = 3,
                               deriv = 2,
                               kernel = kernel,
                               h = mu_J_2_bw)

            # Estimated \mathcal{B} for LLR
            mathcal_B[r] <- mu_J_2 * I_2_K1 / 2

          } else if (bwselect == "IMSE2") {

            # Bandwidth selection for estimating the third and fourth order derivatives
            mu_J_3_bw <- nprobust::lpbwselect(y = J[, r],
                                              x = Z,
                                              eval = zeval[r],
                                              p = 4,
                                              deriv = 3,
                                              kernel = kernel,
                                              bwselect = "mse-dpi")$bws[, 2]

            mu_J_4_bw <- nprobust::lpbwselect(y = J[, r],
                                              x = Z,
                                              eval = zeval[r],
                                              p = 5,
                                              deriv = 4,
                                              kernel = kernel,
                                              bwselect = "mse-dpi")$bws[, 2]

            # Estimate the third and fourth order derivatives
            mu_J_3 <- lpr_func(y = J[, r],
                               x = Z,
                               eval = zeval[r],
                               p = 4,
                               deriv = 3,
                               kernel = kernel,
                               h = mu_J_3_bw)

            mu_J_4 <- lpr_func(y = J[, r],
                               x = Z,
                               eval = zeval[r],
                               p = 5,
                               deriv = 4,
                               kernel = kernel,
                               h = mu_J_4_bw)

            # Estimated \mathcal{B} for LQR
            mathcal_B[r] <- (1 / (24 * kd0_Z[r])) *
              (2 * mu_J_3 * kd1_Z[r] + mu_J_4 * kd0_Z[r]) *
              ((I_4_K1^2 - I_2_K1 * I_6_K1) / (I_4_K1 - I_2_K1^2))

          }
        }

        # Estimate \mathcal{V}_{\theta} for the LLR or LQR estimator -----------

        mathcal_V <- rep(NA, num_zeval)

        for (r in 1:num_zeval) {

          # Bandwidth selection for computing LLR residuals
          mu_J_0_bw <- nprobust::lpbwselect(y = J[, r],
                                            x = Z,
                                            eval = Z_supp,
                                            p = 1,
                                            deriv = 0,
                                            kernel = kernel,
                                            bwselect = "imse-dpi")$bws[1, 2]

          # Compute LLR residuals
          mu_J_0 <- lpr_func(y = J[, r],
                             x = Z,
                             eval = Z,
                             p = 1,
                             deriv = 0,
                             kernel = kernel,
                             h = mu_J_0_bw)

          U_hat[, r] <- J[, r] - mu_J_0

          # Bandwidth selection for estimating sigma^2
          sigma2_bw <- nprobust::lpbwselect(y = U_hat[, r]^2,
                                            x = Z,
                                            eval = zeval[r],
                                            p = 1,
                                            deriv = 0,
                                            kernel = kernel,
                                            bwselect = "mse-dpi")$bws[, 2]

          # Estimate sigma^2
          sigma2 <- lpr_func(y = U_hat[, r]^2,
                             x = Z,
                             eval = zeval[r],
                             p = 1,
                             deriv = 0,
                             kernel = kernel,
                             h = sigma2_bw)

          # Estimated \mathcal{V} for LLR or LQR
          if (bwselect %in% c("IMSE1", "US1")) {

            mathcal_V[r] <- const_V1 * sigma2 / kd0_Z[r]

          } else if (bwselect == "IMSE2") {

            mathcal_V[r] <- const_V2 * sigma2 / kd0_Z[r]

          }
        }

        # Compute bandwidth ----------------------------------------------------

        # Compute the integrated \mathcal{B}^2 and \mathcal{V} using trapezoid formula
        int_mathcal_B <- int_mathcal_V <- 0
        for (r in 1:(num_zeval - 1)) {

          int_mathcal_B <- int_mathcal_B +
            (zeval[r + 1] - zeval[r]) * (mathcal_B[r]^2 + mathcal_B[r + 1]^2) / 2

          int_mathcal_V <- int_mathcal_V +
            (zeval[r + 1] - zeval[r]) * (mathcal_V[r] + mathcal_V[r + 1]) / 2

        }

        # Bandwidth
        bw_temp <-
          (bwselect == "IMSE1") * (int_mathcal_V / (4 * int_mathcal_B))^(1/5) * n^(-1/5) +
          (bwselect == "IMSE2") * (int_mathcal_V / (8 * int_mathcal_B))^(1/9) * n^(-1/9) +
          (bwselect == "US1")   * (int_mathcal_V / (4 * int_mathcal_B))^(1/5) * n^(-2/7)

        # Save
        bw <- c(bw, bw_temp)

      }
    }
  }

  # Bandwidth selection for uniform inference over (e, z)
  if (uniformall | (bwselect == "manual" & length(bw) == 1)) {

    bw <- rep(min(bw), num_eval)

  }

  #=============================================================================
  # Estimate CATT with bandwidth tailored for aggregation
  #=============================================================================

  for (id_eval in 1:num_eval) {

    # Last element of (g, t, .)
    e1 <- eval[id_eval]

    if (type %in% c("dynamic", "group", "calendar")) {

      # Number and index of evaluation points (g, t, .) such that last element = this e
      num_gte   <-   sum(all_gteeval[, 3] == e1)
      index_gte <- which(all_gteeval[, 3] == e1)

      # Evaluation points (g, t, .) such that last element = this e
      gteeval <- all_gteeval[index_gte, ]

    } else if (type == "simple") {

      # Same as above but need modifications when "simple" because e = NA
      num_gte   <- nrow(all_gteeval)
      index_gte <- 1:num_gte
      gteeval   <- all_gteeval

    }

    if (num_gte > 1) {

      # CATT estimation
      output2 <- catt_gt_continuous(yname = yname,
                                    tname = tname,
                                    idname = idname,
                                    gname = gname,
                                    zname = zname,
                                    xformla = xformla,
                                    data = data,
                                    zeval = zeval,
                                    gteval = gteeval[, 1:2],
                                    pretrend = FALSE,
                                    control_group = control_group,
                                    anticipation = anticipation,
                                    alp = alp,
                                    bstrap = FALSE,
                                    biters = 1,
                                    porder = porder,
                                    kernel = kernel,
                                    bwselect = "manual",
                                    bw = bw[id_eval],
                                    uniformall = TRUE)

      # Update lists
      list_catt[[id_eval]]   <- matrix(output2$Estimate$est, nrow = num_zeval, ncol = num_gte)

      list_B_g_t[[id_eval]]  <- output2$B_g_t
      list_G_g[[id_eval]]    <- output2$G_g
      list_mu_G_g[[id_eval]] <- output2$mu_G_g

    }
  }

  #-----------------------------------------------------------------------------
  # Uniform inference for the summary parameter with bandwidth tailored for aggregation
  #-----------------------------------------------------------------------------

  if (bstrap) {

    # Kappa
    kappa <- (sqrt(5) + 1) / 2

    # Multiplier bootstrap weights (Mammen)
    mb_weight <- matrix(sample(x = c(2 - kappa, 1 + kappa),
                               size = n * biters,
                               prob = c(kappa/sqrt(5), 1 - kappa/sqrt(5)),
                               replace = TRUE),
                        nrow = biters,
                        ncol = n)

    # Sup-t-statistics
    mb_sup_t <- matrix(NA, nrow = biters, ncol = num_eval)

  }

  for (id_eval in 1:num_eval) {

    # Last element of (g, t, .)
    e1 <- eval[id_eval]

    if (type %in% c("dynamic", "group", "calendar")) {

      # Number and index of evaluation points (g, t, .) such that last element = this e
      num_gte   <-   sum(all_gteeval[, 3] == e1)
      index_gte <- which(all_gteeval[, 3] == e1)

      # Evaluation points (g, t, .) such that last element = this e
      gteeval <- all_gteeval[index_gte, ]

    } else if (type == "simple") {

      # Same as above but need modifications when "simple" because e = NA
      num_gte   <- nrow(all_gteeval)
      index_gte <- 1:num_gte
      gteeval   <- all_gteeval

    }

    if (num_gte > 1) {

      # Subsets used for aggregation
      B_g_t  <- list_B_g_t[[id_eval]]
      G_g    <- list_G_g[[id_eval]]
      mu_G_g <- list_mu_G_g[[id_eval]]

      # Aggregation weights
      aggte_weight <- matrix(NA, nrow = num_zeval, ncol = num_gte)

      # Estimated \xi_g_t
      xi_g_t <- array(NA, dim = c(n, num_zeval, num_gte))

      # LLR residuals used for standard errors and bootstrap inference
      U_hat <- matrix(NA, nrow = n, ncol = num_zeval)

      for (id_gte in 1:num_gte) {

        if (type %in% c("dynamic", "calendar")) {

          # Estimated weights \hat w_{g,t}^{es}(e, z) or \hat w_{c}(t, z)
          aggte_weight[, id_gte] <- mu_G_g[, id_gte] / rowSums(mu_G_g)

          for (i in 1:n) {

            # Estimated Xi
            xi_g_t[i, , id_gte] <- G_g[i, id_gte] / rowSums(mu_G_g) -
              mu_G_g[, id_gte] / (rowSums(mu_G_g))^2 * sum(G_g[i, ])

          }

        } else if (type == "group") {

          # Weights w_{g,t}^{sel}(g', z)
          aggte_weight[, id_gte] <- 1 / sum(gteeval[, 1] == e1 & e1 <= gteeval[, 2] & gteeval[, 2] < gbar)

          # Estimated Xi
          xi_g_t[, , id_gte] <- 0

        } else if (type == "simple") {

          # Estimated weights \hat w_{g,t}^{OW}(z) * \hat \kappa(z)
          aggte_weight[, id_gte] <- mu_G_g[, id_gte] / rowSums(mu_G_g)

        }
      }

      # CATT estimates
      catt <- list_catt[[id_eval]]

      # Estimated summands in J
      J_temp <- array(NA, dim = c(n, num_zeval, num_gte))

      for (id_gte in 1:num_gte) {

        if (type == "simple") {

          # Estimated \kappa(z)
          aggte_kappa <- rowSums(aggte_weight)

          for (i in 1:n) {

            # Estimated Xi
            xi_g_t[i, , id_gte] <-
              G_g[i, id_gte] / (aggte_kappa * rowSums(mu_G_g)) -
              mu_G_g[, id_gte] / (aggte_kappa * (rowSums(mu_G_g))^2) * sum(G_g[i, ]) -
              mu_G_g[, id_gte] / (aggte_kappa * rowSums(mu_G_g))^2 * rowSums(G_g[i, ] - mu_G_g / rowSums(mu_G_g) * sum(G_g[i, ]))

          }
        }

        for (i in 1:n) {

          # Summands in the estimated J
          J_temp[i, , id_gte] <- aggte_weight[, id_gte] * B_g_t[i, , id_gte] +
            catt[, id_gte] * xi_g_t[i, , id_gte]

        }
      }

      # Compute the estimated J
      J <- matrix(NA, nrow = n, ncol = num_zeval)

      for (r in 1:num_zeval) {
        J[, r] <- rowSums(J_temp[, r, ])
      }

      # Estimated summary parameter
      aggte_est <- rowSums(catt * aggte_weight)

      #-------------------------------------------------------------------------
      # Standard Error
      #-------------------------------------------------------------------------

      mathcal_V <- rep(NA, num_zeval)

      for (r in 1:num_zeval) {

        # Bandwidth selection for computing LLR residuals
        mu_J_0_bw <- nprobust::lpbwselect(y = J[, r],
                                          x = Z,
                                          eval = Z_supp,
                                          p = 1,
                                          deriv = 0,
                                          kernel = kernel,
                                          bwselect = "imse-dpi")$bws[1, 2]

        # Compute LLR residuals
        mu_J_0 <- lpr_func(y = J[, r],
                           x = Z,
                           eval = Z,
                           p = 1,
                           deriv = 0,
                           kernel = kernel,
                           h = mu_J_0_bw)

        U_hat[, r] <- J[, r] - mu_J_0

        # Bandwidth selection for estimating sigma^2
        sigma2_bw <- nprobust::lpbwselect(y = U_hat[, r]^2,
                                          x = Z,
                                          eval = zeval[r],
                                          p = 1,
                                          deriv = 0,
                                          kernel = kernel,
                                          bwselect = "mse-dpi")$bws[, 2]

        # Estimate sigma^2
        sigma2 <- lpr_func(y = U_hat[, r]^2,
                           x = Z,
                           eval = zeval[r],
                           p = 1,
                           deriv = 0,
                           kernel = kernel,
                           h = sigma2_bw)

        # Estimated \mathcal{V}_{\theta}
        mathcal_V[r] <- const_V * sigma2 / kd0_Z[r]

      }

      # Standard error
      se <- sqrt(mathcal_V / (n * bw[id_eval]))

      #-------------------------------------------------------------------------
      # Analytical uniform confidence bands
      #-------------------------------------------------------------------------

      # Analytical critical value
      a <- zeval[1]
      b <- zeval[num_zeval]

      a_hat_sq <- 2 * log((b - a) / bw[id_eval]) + 2 * log(sqrt(lambda) / (2 * pi))

      c_hat <- sqrt(a_hat_sq - 2 * log(log(1 / sqrt(1 - alp))))

      # Analytical uniform confidence bands
      ci1_lower <- aggte_est - c_hat * se
      ci1_upper <- aggte_est + c_hat * se

      #-------------------------------------------------------------------------
      # Bootstrap uniform confidence bands
      #-------------------------------------------------------------------------

      if (bstrap) {

        mb_result <- NULL

        for (r in 1:num_zeval) {

          # u_{ih}
          u_value <- (Z - zeval[r]) / bw[id_eval]

          # Kernel
          kernel_value <- kernel_func(u = u_value, kernel = kernel)

          # Psi
          if (porder == 1) {
            Psi <- 1
          } else if (porder == 2) {
            Psi <- (I_4_K1 - u_value^2 * I_2_K1) / (I_4_K1 - I_2_K1^2)
          }

          # Variables
          mb_est_temp <- mb_t_temp <- rep(NA, length = biters)

          for (mb in 1:biters) {

            # Bootstrap estimates
            mb_est_temp[mb] <- aggte_est[r] + (1 / (kd0_Z[r] * n * bw[id_eval])) *
              sum((mb_weight[mb, ] - 1) * Psi * U_hat[, r] * kernel_value)

            # Bootstrap t statistics
            mb_t_temp[mb] <- abs(mb_est_temp[mb] - aggte_est[r]) / se[r]

          }

          # Store
          mb_result <- rbind(mb_result, c(mb_est_temp, mb_t_temp))

        }

        # Bootstrap estimates
        mb_est <- mb_result[, 1:biters]

        # Bootstrap t statistics
        mb_t <- mb_result[, -(1:biters)]

        # Sup-t-statistic for uniform inference over z
        mb_sup_t[, id_eval] <- apply(X = mb_t, MARGIN = 2, FUN = max, na.rm = TRUE)

      }

      #-------------------------------------------------------------------------
      # Save as lists
      #-------------------------------------------------------------------------

      list_aggte_est[[id_eval]]    <- aggte_est
      list_aggte_weight[[id_eval]] <- aggte_weight
      list_catt[[id_eval]]         <- catt
      list_ci1_lower[[id_eval]]    <- ci1_lower
      list_ci1_upper[[id_eval]]    <- ci1_upper
      list_ci2_lower[[id_eval]]    <- NA
      list_ci2_upper[[id_eval]]    <- NA
      list_J[[id_eval]]            <- J
      list_se[[id_eval]]           <- se
      list_xi_g_t[[id_eval]]       <- xi_g_t

    }
  }

  # Construct bootstrap uniform confidence bands
  if (bstrap) {
    if (uniformall) {

      # Sup-t-statistics for uniform inference over (e, z)
      mb_sup_t <- apply(X = mb_sup_t,
                        MARGIN = 1,
                        FUN = max,
                        na.rm = TRUE)

      # Bootstrap critical value for uniform inference over (e, z)
      c_check <- stats::quantile(x = mb_sup_t, probs = 1 - alp)

      c_check <- rep(c_check, num_eval * num_zeval)

    } else {

      # Bootstrap critical values for uniform inference over z
      c_check <- apply(X = mb_sup_t,
                       MARGIN = 2,
                       FUN = stats::quantile,
                       probs = 1 - alp)

      c_check <- rep(c_check, each = num_zeval)

    }

    for (id_eval in 1:num_eval) {

      # Estimates and standard errors
      aggte_est <- list_aggte_est[[id_eval]]
      se        <- list_se[[id_eval]]

      # Bootstrap uniform confidence bands
      ci2_lower <- aggte_est - c_check[((id_eval - 1) * num_zeval + 1):(id_eval * num_zeval)] * se
      ci2_upper <- aggte_est + c_check[((id_eval - 1) * num_zeval + 1):(id_eval * num_zeval)] * se

      # Save as lists
      list_ci2_lower[[id_eval]] <- ci2_lower
      list_ci2_upper[[id_eval]] <- ci2_upper

    }
  }

  #=============================================================================
  # Save the uniform inference results
  #=============================================================================

  Estimate <- NULL

  for (id_eval in 1:num_eval) {

    # Evaluation point specific to the chosen summary parameter
    e1 <- eval[id_eval]

    # Save the uniform inference results
    Estimate <- rbind(Estimate,
                      cbind(e1,
                            zeval,
                            list_aggte_est[[id_eval]],
                            list_se[[id_eval]],
                            list_ci1_lower[[id_eval]],
                            list_ci1_upper[[id_eval]],
                            list_ci2_lower[[id_eval]],
                            list_ci2_upper[[id_eval]],
                            bw[id_eval]))

  }

  # Names
  colnames(Estimate) <- c("eval",
                          "z",
                          "est",
                          "se",
                          "ci1_lower",
                          "ci1_upper",
                          "ci2_lower",
                          "ci2_upper",
                          "bw")

  rownames(Estimate) <- NULL

  # as.data.frame
  Estimate <- as.data.frame(Estimate)

  #-----------------------------------------------------------------------------
  # Figures
  #-----------------------------------------------------------------------------

  Figure1 <- Figure2 <- list()

  for (id_eval in 1:num_eval) {

    # Specify the evaluation point specific to the chosen summary parameter
    e1 <- eval[id_eval]

    # Figure for analytical UCB
    Figure1[[paste0("eval", e1)]] <- graph_func_aggte_continuous(Estimate = Estimate,
                                                                 type     = type,
                                                                 e1       = e1,
                                                                 bstrap   = FALSE)

    if (bstrap) {

      # Figure for bootstrap UCB
      Figure2[[paste0("eval", e1)]] <- graph_func_aggte_continuous(Estimate = Estimate,
                                                                   type     = type,
                                                                   e1       = e1,
                                                                   bstrap   = TRUE)

    }
  }

  #=============================================================================
  # Return
  #=============================================================================

  return(list(Estimate  = Estimate,
              Figure1   = Figure1,
              Figure2   = Figure2))

}
