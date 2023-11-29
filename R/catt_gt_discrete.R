#' Group-Time Conditional Average Treatment Effects Given a Discrete Covariate
#'
#' `catt_gt_discrete` computes doubly robust uniform confidence bands for
#' the group-time conditional average treatment (CATT) function given a discrete
#' pre-treatment covariate in the staggered difference-in-differences setup of
#' Callaway and Sant'Anna (2021).
#' See Imai, Qin, and Yanagi (2023) for details.
#'
#' @param yname The name of the outcome
#' @param tname The name of the time periods
#' @param idname The name of the cross-sectional IDs
#' @param gname The name of the groups.
#' "G = 0" indicates the never treated group.
#' @param zname The name of the scalar discrete covariate
#' for which the group-time conditional average treatment effects are estimated
#' @param xformla
#' A formula for the covariates to include in the model.
#' It should be of the form `~ X1 + X2`.
#' `xformla` should include `zname` as a covariate.
#' @param data The name of data.frame that contains the balanced panel data
#' @param zeval The vector of the evaluation points z
#' @param gteval The vector or matrix of the evaluation points g and t.
#' If it is a vector, the first and second elements indicate g and t, respectively.
#' If it is a matrix, the first and second columns indicate g's and t's, respectively.
#' Default is `NULL`, and `gteval` is automatically constructed.
#' @param control_group Which units to use the control group.
#' Options are "nevertreated" and "notyettreated".
#' Default is "nevertreated".
#' @param anticipation The number of time periods before participating in the
#' treatment where units can anticipate participating in the treatment and
#' therefore it can affect their untreated potential outcomes.
#' Default is 0.
#' @param alp The significance level. Default is 0.05.
#' @param biters
#' The number of bootstrap iterations to use.
#' Default is 1000.
#' @param uniformall
#' Boolean for whether or not to perform the uniform inference over (g, t, z).
#' Default is `FALSE`, and the uniform inference only over z is performed.
#' @param cores The number of cores to use for parallel processing.
#' The number of available cores can be checked with parallel::detectCores().
#' Default is 1.
#'
#' @return A list that contains the following elements:
#' \item{Estimate}{A data.frame that contains the following elements: \cr
#' g: The group. \cr
#' t: The period. \cr
#' z: The covariate value. \cr
#' est: The doubly robust estimate of CATT. \cr
#' se: The standard error. \cr
#' ci_lower: The lower bound of the UCB via multiplier bootstrapping. \cr
#' ci_upper: The upper bound of the UCB via multiplier bootstrapping.}
#' \item{Figure}{A list that contains the ggplot elements for the UCB}
#'
#' @importFrom dplyr arrange group_by filter mutate pull summarize sym
#' @importFrom foreach %dopar%
#' @importFrom purrr %>%
#' @importFrom rlang UQ
#' @importFrom utils globalVariables
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' data <- datageneration(n = 1000, tau = 4, continuous = FALSE)
#' est <- catt_gt_discrete(yname = "Y",
#'                         tname = "period",
#'                         idname = "id",
#'                         gname = "G",
#'                         zname = "Z",
#'                         xformla = ~ Z,
#'                         data = data,
#'                         zeval = c(-1, 0, 1),
#'                         gteval = c(2, 2),
#'                         control_group = "nevertreated",
#'                         anticipation = 0,
#'                         alp = 0.05,
#'                         biters = 1000,
#'                         uniformall = FALSE,
#'                         cores = 1)
#' }
#'
#' @references
#' Callaway, B., & Sant’Anna, P. H. (2021).
#' Difference-in-differences with multiple time periods.
#' Journal of Econometrics, 225(2), 200-230.
#'
#' Imai, S., Qin, L., & Yanagi, T. (2023).
#' Doubly Robust Uniform Confidence Bands
#' for Group-Time Conditional Average Treatment Effects
#' in Difference-in-Differences.
#' arXiv preprint arXiv:2305.02185.
#'
catt_gt_discrete <- function(yname,
                             tname,
                             idname,
                             gname,
                             zname,
                             xformla,
                             data,
                             zeval,
                             gteval = NULL,
                             control_group = "nevertreated",
                             anticipation = 0,
                             alp = 0.05,
                             biters = 1000,
                             uniformall = FALSE,
                             cores = 1) {

  #-----------------------------------------------------------------------------
  # Error handling
  #-----------------------------------------------------------------------------

  error <- errorhandling_discrete(yname = yname,
                                  tname = tname,
                                  idname = idname,
                                  gname = gname,
                                  zname = zname,
                                  xformla = xformla,
                                  data = data,
                                  zeval = zeval,
                                  gteval = gteval,
                                  control_group = control_group,
                                  anticipation = anticipation,
                                  alp = alp,
                                  biters = biters,
                                  uniformall = uniformall,
                                  cores = cores)

  #-----------------------------------------------------------------------------
  # Basic variable definitions
  #-----------------------------------------------------------------------------

  # Basic variables
  data <- data %>%
    mutate(Y = UQ(sym(yname)),
           period = UQ(sym(tname)),
           id = UQ(sym(idname)),
           G = UQ(sym(gname)),
           Z = UQ(sym(zname))) %>%
    arrange(id, period)

  # Evaluation points (g, t)
  supp_g <- data$G %>%
    unique() %>%
    sort()

  supp_t <- data$period %>%
    unique() %>%
    sort()

  period1 <- supp_t[1]

  geval <- intersect(supp_g, supp_g + anticipation) %>%
    setdiff(0)

  teval <- intersect(supp_t, supp_t - anticipation) %>%
    setdiff(period1)

  gbar <- max(supp_g)

  if (is.null(gteval)) {
    for (g1 in geval) {
      for (t1 in teval) {
        if (control_group == "nevertreated") {
          if (t1 >= g1 - anticipation) {

            gteval <- rbind(gteval, c(g1, t1))

          }
        } else if (control_group == "notyettreated") {
          if (t1 >= g1 - anticipation & t1 < gbar - anticipation) {

            gteval <- rbind(gteval, c(g1, t1))

          }
        }
      }
    }
  }

  # Number of evaluation points (g, t)
  if (is.vector(gteval)) {

    num_gteval <- 1

    uniformall <- FALSE

  } else if (is.matrix(gteval)) {

    num_gteval <- nrow(gteval)

  }

  # Number of individuals
  n <- length(unique(data$id))

  # Number of zeval
  zeval <- sort(zeval)
  num_zeval <- length(zeval)

  # Variables used for parallel computing
  cluster <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cluster, cores = cores)
  packages <- c("dplyr")

  # Matrix of the estimates and uniform confidence bands
  Estimate <- NULL

  # A_g_t
  A_g_t <- array(NA, dim = c(n, num_zeval, num_gteval))

  #-----------------------------------------------------------------------------
  # First-stage estimation
  #-----------------------------------------------------------------------------

  # Parametric estimation
  stage1 <- parametric_func(xformla = xformla,
                            data = data,
                            gteval = gteval,
                            control_group = control_group,
                            anticipation = anticipation,
                            period1 = period1)

  # Generalized propensity score (GPS)
  GPS <- stage1$GPS

  # Outcome regression (OR)
  OR  <- stage1$OR

  # Influence functions
  GPS_inffunc <- stage1$GPS_inffunc
  OR_inffunc  <- stage1$OR_inffunc

  #-----------------------------------------------------------------------------
  # Doubly robust estimation of CATT
  #-----------------------------------------------------------------------------

  for (id_gt in 1:num_gteval) {

    # Specify the values of (g, t)
    if (is.vector(gteval)) {

      g1 <- gteval[1]
      t1 <- gteval[2]

    } else if (is.matrix(gteval)) {

      g1 <- gteval[id_gt, 1]
      t1 <- gteval[id_gt, 2]

    }

    # Make the indicator 1{G = g}, never-/not-yet-treated dummy, GPS, and R_g
    if (control_group == "nevertreated") {

      # GPS
      GPScore <- GPS %>%
        filter(g == g1) %>%
        pull(est)

      # Cross-sectional data at the first period
      data1 <- data %>%
        filter(period == period1) %>%
        mutate(G_g = ifelse(G == g1, 1, 0)) %>%
        mutate(nev = ifelse(G == 0, 1, 0)) %>%
        mutate(GPScore = GPScore) %>%
        mutate(R_g = (GPScore * nev) / (1 - GPScore))

    } else if (control_group == "notyettreated") {

      # GPS
      GPScore <- GPS %>%
        filter(g == g1, t == t1) %>%
        pull(est)

      # Cross-sectional data at the first period
      data1 <- data %>%
        filter(period == period1) %>%
        mutate(G_g = ifelse(G == g1, 1, 0)) %>%
        mutate(notyet = (G == 0) + (G > t1 + anticipation)) %>%
        mutate(GPScore = GPScore) %>%
        mutate(R_g = (GPScore * notyet) / (1 - GPScore))

    }

    # Other variables
    G_g <- data1$G_g

    R_g <- data1$R_g

    Z <- data1$Z

    Y_t <- data %>%
      filter(period == t1) %>%
      pull(Y)

    Y_g <- data %>%
      filter(period == g1 - anticipation - 1) %>%
      pull(Y)

    OREG <- OR %>%
      filter(g == g1, t == t1) %>%
      pull(est)

    Y_diff <- Y_t - Y_g - OREG

    # Doubly robust estimation of CATT -----------------------------------------

    for (id_z in 1:num_zeval) {

      # Specify the value of z
      z <- zeval[id_z]

      # Boolean for Z = z1
      boo_z <- (Z == z)

      # Estimate mu_G_z = E[G_g | Z = z]
      mu_G_z <- mean(G_g[boo_z])

      # Estimate mu_R_z = E[R_g | Z = z]
      mu_R_z <- mean(R_g[boo_z])

      # Compute A_g_t
      A_g_t[, id_z, id_gt] <- ((G_g/mu_G_z) - (R_g/mu_R_z)) * Y_diff

      # DR estimate
      est <- sum(A_g_t[boo_z, id_z, id_gt]) / sum(boo_z)

      #-------------------------------------------------------------------------
      # Save the estimates
      #-------------------------------------------------------------------------

      Estimate <- rbind(Estimate,
                        c(g1, t1, z, est))

    }
  }

  # Names
  colnames(Estimate) <- c("g", "t", "z", "est")
  rownames(Estimate) <- NULL

  # as.data.frame
  Estimate <- as.data.frame(Estimate)

  #-----------------------------------------------------------------------------
  # UCB via multiplier bootstrapping
  #-----------------------------------------------------------------------------

  # Bootstrap results
  mb_result <- NULL

  # Multiplier bootstrap point-wise standard error
  mb_se <- matrix(NA, nrow = num_zeval, ncol = num_gteval)

  # Multiplier bootstrap weights (Mammen)
  kappa <- (sqrt(5) + 1) / 2
  mb_weight <- sample(x = c(2 - kappa, 1 + kappa),
                      size = n * biters,
                      prob = c(kappa/sqrt(5), 1 - kappa/sqrt(5)),
                      replace = TRUE) %>%
    matrix(nrow = n, ncol = biters)

  # Bootstrap
  for (id_gt in 1:num_gteval) {

    # Specify the values of (g, t)
    if (is.vector(gteval)) {

      g1 <- gteval[1]
      t1 <- gteval[2]

    } else if (is.matrix(gteval)) {

      g1 <- gteval[id_gt, 1]
      t1 <- gteval[id_gt, 2]

    }

    for (id_z in 1:num_zeval) {

      # Specify the value of z
      z1 <- zeval[id_z]

      # Boolean for Z = z1
      boo_z1 <- (Z == z1)

      # Original estimates
      est_temp <- Estimate %>%
        filter(g == g1, t == t1, z == z1) %>%
        pull(est)

      # Bootstrap
      mb_result_temp <- NULL
      mb_result_temp <- foreach::foreach(mb = 1:biters,
                                         .combine = 'rbind',
                                         .inorder = TRUE,
                                         .packages = packages) %dopar% {

                                           # Bootstrap DR estimates
                                           mb_est_temp <-
                                             sum(mb_weight[boo_z1, mb] * A_g_t[boo_z1, id_z, id_gt]) / sum(mb_weight[boo_z1, mb])

                                           # Bootstrap residuals
                                           mb_res_temp <- mb_est_temp - est_temp

                                           # Return
                                           return(c(g1, t1, z1, mb, mb_est_temp, mb_res_temp))

                                         }

      # Names
      colnames(mb_result_temp) <- c("g", "t", "z", "mb", "mb_est", "mb_res")
      rownames(mb_result_temp) <- NULL

      # as.data.frame
      mb_result_temp <- as.data.frame(mb_result_temp)

      # Multiplier bootstrap empirical quantiles
      q025 <- stats::quantile(x = mb_result_temp$mb_res, probs = 0.25)
      q075 <- stats::quantile(x = mb_result_temp$mb_res, probs = 0.75)

      # Multiplier bootstrap point-wise standard error
      mb_se[id_z, id_gt] <- (q075 - q025) / (stats::qnorm(0.75) - stats::qnorm(0.25))

      # Multiplier bootstrap t-test statistics
      mb_result_temp <- mb_result_temp %>%
        mutate(mb_t_stat = mb_res / mb_se[id_z, id_gt]) %>%
        select(g, t, z, mb, mb_est, mb_res, mb_t_stat)

      # Save the bootstrap results
      mb_result <- rbind(mb_result, mb_result_temp)

    }
  }

  # Names
  colnames(mb_result) <- c("g", "t", "z", "mb", "mb_est", "mb_res", "mb_t_stat")
  rownames(mb_result) <- NULL

  # as.data.frame
  mb_result <- as.data.frame(mb_result)

  # Max-t-test statistics over z
  max_t_stat1 <- mb_result %>%
    group_by(g, t, mb) %>%
    mutate(max_t_stat = max(abs(mb_t_stat)),
           .groups = "drop")

  # Max-t-test statistics over (g, t, z)
  max_t_stat2 <- mb_result %>%
    group_by(mb) %>%
    mutate(max_t_stat = max(abs(mb_t_stat)),
           .groups = "drop")

  # Uniform critical values over z
  mb_cv1 <- max_t_stat1 %>%
    group_by(g, t) %>%
    summarize(mb_cv = stats::quantile(max_t_stat, probs = 1 - alp),
              .groups = "drop") %>%
    pull(mb_cv)

  # Uniform critical values over (g, t, z)
  mb_cv2 <- max_t_stat2 %>%
    pull(max_t_stat) %>%
    stats::quantile(probs = 1 - alp)

  # UCB over (g, t, z) or over z
  if (uniformall) {

    ci_lower <- Estimate$est - mb_cv2 * c(mb_se)
    ci_upper <- Estimate$est + mb_cv2 * c(mb_se)

  } else {

    ci_lower <- Estimate$est - rep(mb_cv1, each = num_zeval) * c(mb_se)
    ci_upper <- Estimate$est + rep(mb_cv1, each = num_zeval) * c(mb_se)

  }

  #-----------------------------------------------------------------------------
  # Results
  #-----------------------------------------------------------------------------

  Estimate <- Estimate %>%
    mutate(se = c(mb_se),
           ci_lower = ci_lower,
           ci_upper = ci_upper) %>%
    select(g, t, z, est, se, ci_lower, ci_upper)

  #-----------------------------------------------------------------------------
  # Figures
  #-----------------------------------------------------------------------------

  Figure <- list()

  for (id_gt in 1:num_gteval) {

    # Specify the values of (g, t)
    if (is.vector(gteval)) {

      g1 <- gteval[1]
      t1 <- gteval[2]

    } else if (is.matrix(gteval)) {

      g1 <- gteval[id_gt, 1]
      t1 <- gteval[id_gt, 2]

    }

    # Figure for multiplier bootstrap UCB
    Figure[[paste0("g", g1, "_t", t1)]] <- graph_func_discrete(Estimate = Estimate,
                                                               g1 = g1,
                                                               t1 = t1)

  }

  # Stop parallel computing ----------------------------------------------------

  parallel::stopCluster(cluster)

  #-----------------------------------------------------------------------------
  # Return
  #-----------------------------------------------------------------------------

  return(list(Estimate = Estimate,
              Figure   = Figure))

}
