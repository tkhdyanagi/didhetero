#' Group-Time Conditional Average Treatment Effects Given a Continuous Covariate
#'
#' `catt_gt_continuous` computes doubly robust uniform confidence bands for
#' the group-time conditional average treatment (CATT) function given a continuous
#' pre-treatment covariate in the staggered difference-in-differences setup of
#' Callaway and Sant'Anna (2021).
#' See Imai, Qin, and Yanagi (2023) for details.
#'
#' @param yname The name of the outcome
#' @param tname The name of the time periods
#' @param idname The name of the cross-sectional IDs
#' @param gname The name of the groups.
#' "G = 0" indicates the never treated group.
#' @param zname The name of the scalar continuous covariate
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
#' @param bstrap
#' Boolean for whether or not to perform the multiplier bootstrap inference.
#' Default is `TRUE`.
#' If bstrap is `FALSE`, only the analytical critical value is used.
#' @param biters
#' The number of bootstrap iterations to use.
#' Default is 1000, which is only applicable if bstrap is `TRUE`.
#' @param porder
#' The polynomial order used for the second- and third-stage estimation.
#' Options are 1 and 2,
#' which correspond to the local linear and quadratic regressions, respectively.
#' Default is 2.
#' @param kernel
#' The kernel function used for the local polynomial regressions.
#' Options are "gau" for the Gaussian kernel and
#' "epa" for the Epanechnikov kernel.
#' Default is "gau".
#' @param bw
#' The scalar bandwidth used for the second- and third-stage estimation.
#' Default is `NULL`, and the bandwidth is automatically selected.
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
#' ci1_lower: The lower bound of the analytical UCB. \cr
#' ci1_upper: The upper bound of the analytical UCB. \cr
#' ci2_lower: The lower bound of the UCB via multiplier bootstrapping. \cr
#' ci2_upper: The upper bound of the UCB via multiplier bootstrapping. \cr
#' bw: The bandwidth.}
#' \item{Figure1}{A list that contains the ggplot elements for the analytical UCB}
#' \item{Figure2}{A list that contains the ggplot elements for the UCB via multiplier bootstrapping}
#'
#' @importFrom dplyr arrange filter mutate pull sym
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
#' data <- datageneration(n = 1000, tau = 4, continuous = TRUE)
#' est <- catt_gt_continuous(yname = "Y",
#'                           tname = "period",
#'                           idname = "id",
#'                           gname = "G",
#'                           zname = "Z",
#'                           xformla = ~ Z,
#'                           data = data,
#'                           zeval = seq(-1, 1, by = 0.1),
#'                           gteval = c(2, 2),
#'                           control_group = "nevertreated",
#'                           anticipation = 0,
#'                           alp = 0.05,
#'                           bstrap = TRUE,
#'                           biters = 1000,
#'                           porder = 2,
#'                           kernel = "gau",
#'                           bw = NULL,
#'                           uniformall = FALSE,
#'                           cores = 1)
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
catt_gt_continuous <- function(yname,
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
                               bstrap = TRUE,
                               biters = 1000,
                               porder = 2,
                               kernel = "gau",
                               bw = NULL,
                               uniformall = FALSE,
                               cores = 1) {

  #-----------------------------------------------------------------------------
  # Error handling
  #-----------------------------------------------------------------------------

  error <- errorhandling_continuous(yname = yname,
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
                                    bstrap = bstrap,
                                    biters = biters,
                                    porder = porder,
                                    kernel = kernel,
                                    bw = bw,
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
  packages <- c("nprobust")

  # Matrix of the estimates and uniform confidence bands
  Estimate <- NULL

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

  # Constant used for estimating \mathcal{V}_{p,g,t,\delta}
  if (porder == 1) {

    const_V <- I_0_K2

  } else if (porder == 2) {

    const_V <-
      (I_4_K1^2 * I_0_K2 - 2 * I_2_K1 * I_4_K1 * I_2_K2 + I_2_K1^2 * I_4_K2) /
      (I_4_K1 - I_2_K1^2)^2

  }

  # A_g_t and B_g_t
  A_g_t <- B_g_t <- array(NA, dim = c(n, num_zeval, num_gteval))

  # Matrices for mu_E_g_t and mu_F_g_t
  mu_E_g_t <- mu_F_g_t <- matrix(NA, nrow = num_zeval, ncol = num_gteval)

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

  #-----------------------------------------------------------------------------
  # Kernel density estimation
  #-----------------------------------------------------------------------------

  # Z
  Z <- data %>%
    filter(period == period1) %>%
    pull(Z)

  # Support
  Z_supp <- seq(min(Z), max(Z), length = 100)

  # Bandwidth selection
  kd_Z_bw <- nprobust::kdbwselect(x = Z,
                                  eval = zeval,
                                  kernel = "epa",
                                  bwselect = "mse-dpi")$bws[, 2]
  # Kernel density estimate
  kd_Z <- kde_func(x = Z,
                   eval = zeval,
                   kernel = "epa",
                   h = kd_Z_bw)

  #-----------------------------------------------------------------------------
  # Bandwidth selection for second and third stage estimation
  #-----------------------------------------------------------------------------

  bwselect <- bw

  if (is.null(bwselect)) {

    for (id_gt in 1:num_gteval) {

      # Initialization
      mathcal_B <- mathcal_V <- NULL

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

      E_g_t <- R_g * Y_diff
      F_g_t <- G_g * Y_diff

      # LLR estimation for mu_G_g, mu_R_g, mu_E_g_t, and mu_F_g_t --------------

      # Bandwidth selection
      mu_G_g_bw <- nprobust::lpbwselect(y = G_g,
                                        x = Z,
                                        eval = zeval,
                                        p = 1,
                                        deriv = 0,
                                        kernel = kernel,
                                        bwselect = "mse-dpi")$bws[, 2]

      mu_R_g_bw <- nprobust::lpbwselect(y = R_g,
                                        x = Z,
                                        eval = zeval,
                                        p = 1,
                                        deriv = 0,
                                        kernel = kernel,
                                        bwselect = "mse-dpi")$bws[, 2]

      mu_E_g_t_bw <- nprobust::lpbwselect(y = E_g_t,
                                          x = Z,
                                          eval = zeval,
                                          p = 1,
                                          deriv = 0,
                                          kernel = kernel,
                                          bwselect = "mse-dpi")$bws[, 2]

      mu_F_g_t_bw <- nprobust::lpbwselect(y = F_g_t,
                                          x = Z,
                                          eval = zeval,
                                          p = 1,
                                          deriv = 0,
                                          kernel = kernel,
                                          bwselect = "mse-dpi")$bws[, 2]

      # LLR estimation
      mu_G_g <- lpr_func(y = G_g,
                         x = Z,
                         eval = zeval,
                         p = 1,
                         deriv = 0,
                         kernel = kernel,
                         h = mu_G_g_bw)

      mu_R_g <- lpr_func(y = R_g,
                         x = Z,
                         eval = zeval,
                         p = 1,
                         deriv = 0,
                         kernel = kernel,
                         h = mu_R_g_bw)

      mu_E_g_t[, id_gt] <- lpr_func(y = E_g_t,
                                    x = Z,
                                    eval = zeval,
                                    p = 1,
                                    deriv = 0,
                                    kernel = kernel,
                                    h = mu_E_g_t_bw)

      mu_F_g_t[, id_gt] <- lpr_func(y = F_g_t,
                                    x = Z,
                                    eval = zeval,
                                    p = 1,
                                    deriv = 0,
                                    kernel = kernel,
                                    h = mu_F_g_t_bw)

      # Variable definitions ---------------------------------------------------

      for (r in 1:num_zeval) {

        p_diff <- (G_g / mu_G_g[r]) - (R_g / mu_R_g[r])

        A_g_t[, r, id_gt] <- p_diff * Y_diff

        B_g_t[, r, id_gt] <- p_diff * Y_diff +
          ((mu_E_g_t[r, id_gt] / mu_R_g[r]^2) * R_g) -
          ((mu_F_g_t[r, id_gt] / mu_G_g[r]^2) * G_g)

      }

      # Estimate \mathcal{B}_{1,g,t,\delta} for the LLR estimator --------------

      mathcal_B <- foreach::foreach(r = 1:num_zeval,
                                    .combine = 'c',
                                    .export = c('kernel_func', 'lpr_func', 'polynomial_func'),
                                    .inorder = TRUE,
                                    .packages = packages) %dopar% {

                                      # Bandwidth selection for estimating the second order derivative
                                      mu_B_g_t_2_bw <- nprobust::lpbwselect(y = B_g_t[, r, id_gt],
                                                                            x = Z,
                                                                            eval = zeval[r],
                                                                            p = 3,
                                                                            deriv = 2,
                                                                            kernel = kernel,
                                                                            bwselect = "mse-dpi")$bws[, 2]

                                      # Estimate the second order derivative
                                      mu_B_g_t_2 <- lpr_func(y = B_g_t[, r, id_gt],
                                                             x = Z,
                                                             eval = zeval[r],
                                                             p = 3,
                                                             deriv = 2,
                                                             kernel = kernel,
                                                             h = mu_B_g_t_2_bw)

                                      # Estimated \mathcal{B} for LLR
                                      mu_B_g_t_2 * I_2_K1 / 2

                                    }


      # Estimate \mathcal{V}_{1,g,t,\delta} for the LLR estimator --------------

      mathcal_V <- foreach::foreach(r = 1:num_zeval,
                                    .combine = 'c',
                                    .export = c('kernel_func', 'lpr_func', 'polynomial_func'),
                                    .inorder = TRUE,
                                    .packages = packages) %dopar% {

                                      # Bandwidth selection for computing LLR residuals
                                      mu_B_g_t_0_bw <- nprobust::lpbwselect(y = B_g_t[, r, id_gt],
                                                                            x = Z,
                                                                            eval = Z_supp,
                                                                            p = 1,
                                                                            deriv = 0,
                                                                            kernel = kernel,
                                                                            bwselect = "imse-dpi")$bws[1, 2]

                                      # Compute LLR residuals
                                      mu_B_g_t_0 <- lpr_func(y = B_g_t[, r, id_gt],
                                                             x = Z,
                                                             eval = Z,
                                                             p = 1,
                                                             deriv = 0,
                                                             kernel = kernel,
                                                             h = mu_B_g_t_0_bw)

                                      U_hat <- B_g_t[, r, id_gt] - mu_B_g_t_0

                                      # Bandwidth selection for estimating sigma^2
                                      sigma2_bw <- nprobust::lpbwselect(y = U_hat^2,
                                                                        x = Z,
                                                                        eval = zeval[r],
                                                                        p = 1,
                                                                        deriv = 0,
                                                                        kernel = kernel,
                                                                        bwselect = "mse-dpi")$bws[, 2]

                                      # Estimate sigma^2
                                      sigma2 <- lpr_func(y = U_hat^2,
                                                         x = Z,
                                                         eval = zeval[r],
                                                         p = 1,
                                                         deriv = 0,
                                                         kernel = kernel,
                                                         h = sigma2_bw)

                                      # Estimated \mathcal{V} for LLR
                                      I_0_K2 * sigma2 / kd_Z[r]

                                    }

      # LLR-IMSE-optimal bandwidth and under-smoothing bandwidth ---------------

      # Compute the integrated \mathcal{B}^2 and \mathcal{V} using trapezoid formula
      int_mathcal_B <- int_mathcal_V <- 0
      for (r in 1:(num_zeval - 1)) {

        int_mathcal_B <- int_mathcal_B +
          (zeval[r + 1] - zeval[r]) * (mathcal_B[r]^2 + mathcal_B[r + 1]^2) / 2

        int_mathcal_V <- int_mathcal_V +
          (zeval[r + 1] - zeval[r]) * (mathcal_V[r] + mathcal_V[r + 1]) / 2

      }

      # LLR-IMSE-optimal bandwidth
      bw_imse <- (int_mathcal_V / (4 * int_mathcal_B))^(1/5) * n^(-1/5)

      # Under-smoothing bandwidth
      bw_us <- bw_imse * n^(1/5) * n^(-2/7)

      # bandwidth selection for LLR and LQR
      bw <- c(bw,
              bw_us * (porder == 1) + bw_imse * (porder == 2))

    }
  }

  # Bandwidth selection for uniform inference over (g, t, z)
  if (uniformall) {

    bw <- rep(min(bw), num_gteval)

  }

  #-----------------------------------------------------------------------------
  # Second and third stage estimation
  #-----------------------------------------------------------------------------

  for (id_gt in 1:num_gteval) {

    # Variable definitions
    est <- mathcal_V <- NULL

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

    # Second stage estimation: p-th order LPR estimation -----------------------

    # Estimation of mu_E_g_t and mu_F_g_t (if needed)
    if (!is.null(bwselect)) {

      E_g_t <- R_g * Y_diff
      F_g_t <- G_g * Y_diff

      mu_E_g_t_bw <- nprobust::lpbwselect(y = E_g_t,
                                          x = Z,
                                          eval = zeval,
                                          p = 1,
                                          deriv = 0,
                                          kernel = kernel,
                                          bwselect = "mse-dpi")$bws[, 2]

      mu_F_g_t_bw <- nprobust::lpbwselect(y = F_g_t,
                                          x = Z,
                                          eval = zeval,
                                          p = 1,
                                          deriv = 0,
                                          kernel = kernel,
                                          bwselect = "mse-dpi")$bws[, 2]

      mu_E_g_t[, id_gt] <- lpr_func(y = E_g_t,
                                    x = Z,
                                    eval = zeval,
                                    p = 1,
                                    deriv = 0,
                                    kernel = kernel,
                                    h = mu_E_g_t_bw)

      mu_F_g_t[, id_gt] <- lpr_func(y = F_g_t,
                                    x = Z,
                                    eval = zeval,
                                    p = 1,
                                    deriv = 0,
                                    kernel = kernel,
                                    h = mu_F_g_t_bw)

    }


    # LPR estimation
    mu_G_g <- lpr_func(y = G_g,
                       x = Z,
                       eval = zeval,
                       p = porder,
                       deriv = 0,
                       kernel = kernel,
                       h = bw[id_gt])

    mu_R_g <- lpr_func(y = R_g,
                       x = Z,
                       eval = zeval,
                       p = porder,
                       deriv = 0,
                       kernel = kernel,
                       h = bw[id_gt])

    # Variable definitions
    for (r in 1:num_zeval) {

      p_diff <- (G_g / mu_G_g[r]) - (R_g / mu_R_g[r])

      A_g_t[, r, id_gt] <- p_diff * Y_diff

      B_g_t[, r, id_gt] <- p_diff * Y_diff +
        ((mu_E_g_t[r, id_gt] / mu_R_g[r]^2) * R_g) -
        ((mu_F_g_t[r, id_gt] / mu_G_g[r]^2) * G_g)

    }

    # Third stage estimation: p-th order LPR estimation
    est <- foreach::foreach(r = 1:num_zeval,
                            .combine = 'c',
                            .export = c('kernel_func', 'lpr_func', 'polynomial_func'),
                            .inorder = TRUE,
                            .packages = packages) %dopar% {

                              lpr_func(y = A_g_t[, r, id_gt],
                                       x = Z,
                                       eval = zeval[r],
                                       p = porder,
                                       deriv = 0,
                                       kernel = kernel,
                                       h = bw[id_gt])

                            }

    #---------------------------------------------------------------------------
    # Standard error
    #---------------------------------------------------------------------------

    # Estimate \mathcal{V}_{p,g,t,\delta}
    mathcal_V <- foreach::foreach(r = 1:num_zeval,
                                  .combine = 'c',
                                  .export = c('kernel_func', 'lpr_func', 'polynomial_func'),
                                  .inorder = TRUE,
                                  .packages = packages) %dopar% {

                                    # Bandwidth selection for computing LPR residuals
                                    mu_B_g_t_0_bw <- nprobust::lpbwselect(y = B_g_t[, r, id_gt],
                                                                          x = Z,
                                                                          eval = Z_supp,
                                                                          p = porder,
                                                                          deriv = 0,
                                                                          kernel = kernel,
                                                                          bwselect = "imse-dpi")$bws[1, 2]

                                    # Compute LPR residuals
                                    mu_B_g_t_0 <- lpr_func(y = B_g_t[, r, id_gt],
                                                           x = Z,
                                                           eval = Z,
                                                           p = porder,
                                                           deriv = 0,
                                                           kernel = kernel,
                                                           h = mu_B_g_t_0_bw)

                                    U_hat <- B_g_t[, r, id_gt] - mu_B_g_t_0

                                    # Bandwidth selection for estimating sigma^2
                                    sigma2_bw <- nprobust::lpbwselect(y = U_hat^2,
                                                                      x = Z,
                                                                      eval = zeval[r],
                                                                      p = 1,
                                                                      deriv = 0,
                                                                      kernel = kernel,
                                                                      bwselect = "mse-dpi")$bws[, 2]

                                    # Estimate sigma^2
                                    sigma2 <- lpr_func(y = U_hat^2,
                                                       x = Z,
                                                       eval = zeval[r],
                                                       p = 1,
                                                       deriv = 0,
                                                       kernel = kernel,
                                                       h = sigma2_bw)

                                    # Estimated \mathcal{V}_{p,g,t,\delta}
                                    const_V * sigma2 / kd_Z[r]

                                  }

    # Standard error
    se <- sqrt(mathcal_V / (n * bw[id_gt]))

    #---------------------------------------------------------------------------
    # Analytical uniform confidence bands
    #---------------------------------------------------------------------------

    # Analytical critical value
    a <- zeval[1]
    b <- zeval[num_zeval]

    a_hat_sq <- 2 * log((b - a) / bw[id_gt]) + 2 * log(sqrt(lambda) / (2 * pi))

    c_hat <- sqrt(a_hat_sq - 2 * log(log(1 / sqrt(1 - alp))))

    # Analytical uniform confidence bands
    ci1_lower <- est - c_hat * se
    ci1_upper <- est + c_hat * se

    #---------------------------------------------------------------------------
    # Save the estimates
    #---------------------------------------------------------------------------

    Estimate <- rbind(Estimate,
                      cbind(g1,
                            t1,
                            zeval,
                            est,
                            se,
                            ci1_lower,
                            ci1_upper,
                            bw[id_gt]))

  }

  # Names
  colnames(Estimate) <- c("g",
                          "t",
                          "z",
                          "est",
                          "se",
                          "ci1_lower",
                          "ci1_upper",
                          "bw")

  rownames(Estimate) <- NULL

  # as.data.frame
  Estimate <- as.data.frame(Estimate)

  #-----------------------------------------------------------------------------
  # Uniform confidence bands via multiplier bootstrapping
  #-----------------------------------------------------------------------------

  if (bstrap) {

    # Multiplier bootstrap LPR estimates and t statistics
    mb_est <- mb_t <- array(NA, dim = c(biters, num_zeval, num_gteval))

    # Sup-t-statistics
    mb_sup_t <- matrix(NA, nrow = biters, ncol = num_gteval)

    # Kappa
    kappa <- (sqrt(5) + 1) / 2

    # Multiplier bootstrap
    for (mb in 1:biters) {

      # Multiplier bootstrap weights (Mammen)
      mb_weight <- sample(x = c(2 - kappa, 1 + kappa),
                          size = n,
                          prob = c(kappa/sqrt(5), 1 - kappa/sqrt(5)),
                          replace = TRUE)

      for (id_gt in 1:num_gteval) {

        # Specify the values of (g, t)
        if (is.vector(gteval)) {

          g1 <- gteval[1]
          t1 <- gteval[2]

        } else if (is.matrix(gteval)) {

          g1 <- gteval[id_gt, 1]
          t1 <- gteval[id_gt, 2]

        }

        # Original estimates
        est_temp <- Estimate %>%
          filter(g == g1, t == t1) %>%
          pull(est)

        # Original standard errors
        se_temp <- Estimate %>%
          filter(g == g1, t == t1) %>%
          pull(se)

        # Bootstrapped LPR estimate
        mb_est[mb, , id_gt] <- foreach::foreach(r = 1:num_zeval,
                                                .combine = 'c',
                                                .export = c('kernel_func', 'lpr_func', 'polynomial_func'),
                                                .inorder = TRUE,
                                                .packages = packages) %dopar% {

                                                  lpr_func(y = A_g_t[, r, id_gt],
                                                           x = Z,
                                                           eval = zeval[r],
                                                           p = porder,
                                                           deriv = 0,
                                                           kernel = kernel,
                                                           h = bw[id_gt],
                                                           weight = mb_weight)

                                                }

        # Multiplier bootstrap t statistics
        mb_t[mb, , id_gt] <- abs(mb_est[mb, , id_gt] - est_temp) / se_temp

        # Sup-t-statistic for uniform inference over z
        mb_sup_t[mb, id_gt] <- max(mb_t[mb, , id_gt], na.rm = TRUE)

      }
    }

    if (uniformall) {

      # Sup-t-statistic for uniform inference over (g, t, z)
      mb_sup_t <- apply(X = mb_sup_t,
                        MARGIN = 1,
                        FUN = max,
                        na.rm = TRUE)

      # Multiplier bootstrap critical value for uniform inference over (g, t, z)
      c_check <- stats::quantile(x = mb_sup_t, probs = 1 - alp)

    } else {

      # Multiplier bootstrap critical values for uniform inference over z
      c_check <- apply(X = mb_sup_t,
                       MARGIN = 2,
                       FUN = stats::quantile,
                       probs = 1 - alp)

      c_check <- rep(c_check, each = num_zeval)

    }

    # Multiplier bootstrap uniform confidence bands
    ci2_lower <- Estimate$est - c_check * Estimate$se
    ci2_upper <- Estimate$est + c_check * Estimate$se

  } else if (!bstrap) {

    ci2_lower <- ci2_upper <- NA

  }

  # Mutate and arrange
  Estimate <- Estimate %>%
    mutate(ci2_lower = ci2_lower,
           ci2_upper = ci2_upper) %>%
    select(g, t, z, est, se, ci1_lower, ci1_upper, ci2_lower, ci2_upper, bw)

  #-----------------------------------------------------------------------------
  # Figures
  #-----------------------------------------------------------------------------

  Figure1 <- Figure2 <- list()

  for (id_gt in 1:num_gteval) {

    # Specify the values of (g, t)
    if (is.vector(gteval)) {

      g1 <- gteval[1]
      t1 <- gteval[2]

    } else if (is.matrix(gteval)) {

      g1 <- gteval[id_gt, 1]
      t1 <- gteval[id_gt, 2]

    }

    # Figure for analytical UCB
    Figure1[[paste0("g", g1, "_t", t1)]] <- graph_func_continuous(Estimate = Estimate,
                                                                  g1 = g1,
                                                                  t1 = t1,
                                                                  bstrap = FALSE)

    if (bstrap) {

      # Figure for multiplier bootstrap UCB
      Figure2[[paste0("g", g1, "_t", t1)]] <- graph_func_continuous(Estimate = Estimate,
                                                                    g1 = g1,
                                                                    t1 = t1,
                                                                    bstrap = TRUE)

    }
  }

  # Stop parallel computing ----------------------------------------------------

  parallel::stopCluster(cluster)

  # Return ---------------------------------------------------------------------

  return(list(Estimate = Estimate,
              Figure1  = Figure1,
              Figure2  = Figure2))

}
