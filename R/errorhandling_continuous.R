#' Error handling for a continuous covariate
#'
#' @param yname The name of the outcome
#' @param tname The name of the time periods
#' @param idname The name of the cross-sectional IDs
#' @param gname The name of the groups
#' @param zname The name of the scalar continuous covariate
#' for which the group-time conditional average treatment effects are estimated
#' @param xformla
#' A formula for the covariates to include in the model
#' @param data The name of data.frame that contains the data
#' @param zeval The vector of the evaluation points z
#' @param gteval The vector or matrix of the evaluation points g and t
#' @param control_group Which units to use the control group
#' @param anticipation The number of time periods before participating in the
#' treatment where units can anticipate participating in the treatment and
#' therefore it can affect their untreated potential outcomes
#' @param alp The significance level. Default is 0.05
#' @param bstrap
#' Boolean for whether or not to perform the multiplier bootstrap inference
#' @param biters
#' The number of bootstrap iterations to use
#' @param porder
#' The polynomial order used for the second- and third-stage estimation
#' @param kernel
#' The kernel function used for the local polynomial regressions
#' @param bw
#' The scalar bandwidth used for the second- and third-stage estimation
#' @param uniformall
#' Boolean for whether or not to perform the uniform inference over (g, t, z)
#' @param cores The number of cores to use for parallel processing
#'
#' @returns NULL if there are no errors
#'
#' @noRd
#'
errorhandling_continuous <- function(yname,
                                     tname,
                                     idname,
                                     gname,
                                     zname,
                                     xformla,
                                     data,
                                     zeval,
                                     gteval,
                                     control_group,
                                     anticipation,
                                     alp,
                                     bstrap,
                                     biters,
                                     porder,
                                     kernel,
                                     bw,
                                     uniformall,
                                     cores) {

  # yname, tname, idname, gname, zname -----------------------------------------

  if (!is.character(yname)  | length(yname)  != 1 |
      !is.character(tname)  | length(tname)  != 1 |
      !is.character(idname) | length(idname) != 1 |
      !is.character(gname)  | length(gname)  != 1 |
      !is.character(zname)  | length(zname)  != 1) {

    # Character
    stop(paste("yname, tname, idname, gname, and zname must be characters."))

  }

  # xformla --------------------------------------------------------------------

  if (!rlang::is_formula(xformla)) {

    # Formula
    stop(paste("xformla must be a formula of the form `xformla = ~ X1 + X2`."))

  }

  # data -----------------------------------------------------------------------

  if (!is.data.frame(data)) {

    # data.frame
    stop(paste("data must be a data.frame."))

  }

  # zeval ----------------------------------------------------------------------

  if (!is.numeric(zeval) | !is.vector(zeval)) {

    # Numeric
    stop(paste("zeval must be a numeric scalar or vector."))

  }

  # gteval ---------------------------------------------------------------------

  if (!is.null(gteval)) {
    if (!is.numeric(gteval)) {
      if (!is.vector(gteval) & !is.matrix(gteval)) {

        # Numeric
        stop(paste("gteval must be a numeric vector or matrix."))

      }
    }
  }

  # control_group --------------------------------------------------------------

  if (control_group != "nevertreated" & control_group != "notyettreated") {

    # Options
    stop(paste("control_group must be 'nevertreated' or 'notyettreated'."))

  }

  # anticipation ---------------------------------------------------------------

  if (!is.numeric(anticipation) | length(anticipation) != 1) {

    # Numeric
    stop(paste("anticipation must be an integer."))

  }

  # alp ------------------------------------------------------------------------

  if (!is.numeric(alp)) {

    # Numeric
    stop(paste("alp must be a positive number between 0 and 1."))

  } else {

    # Scalar
    if (length(alp) != 1) {
      stop(paste("alp must be a positive number between 0 and 1."))
    }

    # Positive value between 0 and 1
    if (alp <= 0 | alp >= 1) {
      stop(paste("alp must be a positive number between 0 and 1."))
    }

  }

  # bstrap ---------------------------------------------------------------------

  if (!is.logical(bstrap)) {

    # Boolean
    stop(paste("bstrap must be Boolean."))

  }

  # biters ---------------------------------------------------------------------

  # Positive scalar
  if (length(biters) != 1 | biters <= 0) {

    stop(paste("biters must be a positive number."))

  }

  # porder ---------------------------------------------------------------------

  # Options
  if (porder != 1 & porder != 2) {

    stop(paste("porder must be 1 or 2."))

  }

  # kernel ---------------------------------------------------------------------

  # Options
  if (kernel != "epa" & kernel != "gau") {

    stop(paste("kernel must be 'epa' or 'gau'."))

  }

  # bw -------------------------------------------------------------------------

  # Positive number
  if (!is.null(bw)) {
    if (!is.numeric(bw) | length(bw) != 1 | bw <= 0) {

      stop(paste("bw must be a positive number"))

    }
  }

  # uniformall -----------------------------------------------------------------

  if (!is.logical(uniformall)) {

    # Boolean
    stop(paste("uniformall must be Boolean."))

  }

  # cores ----------------------------------------------------------------------

  # Positive number
  if (!is.numeric(cores) | length(cores) != 1 | cores <= 0) {

    stop(paste("cores must be a positive integer"))

  }

  # Return ---------------------------------------------------------------------

  return(NULL)

}
