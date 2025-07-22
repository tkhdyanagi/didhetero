#' Error handling for aggregation with a continuous covariate
#'
#' @param output The output of the `catt_gt_continuous` function. 
#' @param type Which type of the summary parameter is of interest. 
#' @param eval The vector of the evaluation point specific to the chosen summary parameter.
#' @param bstrap Boolean for whether or not to perform weighted bootstrapping.
#' @param biters The number of bootstrap iterations.
#' @param bwselect The bandwidth selection method used for the aggregation.
#' @param bw The bandwidth used for the aggregation.
#' @param uniformall Boolean for whether or not to perform the uniform inference over both eval and z.
#'
#' @returns NULL if there are no errors
#'
#' @noRd
#'
errorhandling_aggte_continuous <- function(output,
                                           type,
                                           eval,
                                           bstrap,
                                           biters,
                                           bwselect,
                                           bw,
                                           uniformall) {
  
  # output ---------------------------------------------------------------------
  
  # type -----------------------------------------------------------------------
  
  if (type != "simple" & type != "dynamic" & type != "group" & type != "calendar") {
    
    stop(paste("type must be 'simple', 'dynamic', 'group', or 'calendar'."))
    
  }
  
  # eval ---------------------------------------------------------------------
  
  if (!is.null(eval)) {
    if (!is.numeric(eval)) {
      if (!is.vector(eval)) {
        
        # Numeric
        stop(paste("eval must be NULL or a numeric vector."))
        
      }
    }
  }
  
  # bstrap ---------------------------------------------------------------------
  
  if (!is.logical(bstrap)) {
    
    # Boolean
    stop(paste("bstrap must be Boolean."))
    
  }
  
  # biters ---------------------------------------------------------------------
  
  # Positive scalar
  if (bstrap) {
    if (length(biters) != 1 | biters <= 0) {
      
      stop(paste("When bstrap = TRUE, biters must be a positive number."))
      
    }
  }
  
  # bwselect -------------------------------------------------------------------
  
  # Options
  if (bwselect != "IMSE1" & bwselect != "IMSE2" & bwselect != "US1" & bwselect != "manual") {
    
    stop(paste("bwselect must be 'IMSE1', 'IMSE2', 'US1', or 'manual'."))
    
  }
  
  # For "manual"
  if (bwselect == "manual" & is.null(bw)) {
    
    stop(paste("bw must be specified manually when bwselect = 'manual'."))
    
  }
  
  # For other options
  if (bwselect != "manual" & !is.null(bw)) {
    
    stop(paste("bw can be specified manually only when bwselect = 'manual'."))
    
  }
  
  # bw -------------------------------------------------------------------------
  
  if (!is.null(bw)) {
    
    # Positive number
    if (!is.numeric(bw)) {
      
      stop(paste("bw must be a positive scalar or vector whose length equals to the number of eval."))
      
    }
    
    # Length of bw
    if (length(bw) != length(eval)) {
      
      stop(paste("bw must be a positive scalar or vector whose length equals to the number of eval."))
      
    }
  }
  
  # uniformall -----------------------------------------------------------------
  
  if (!is.logical(uniformall)) {
    
    # Boolean
    stop(paste("uniformall must be Boolean."))
    
  }
  
  # Return ---------------------------------------------------------------------
  
  return(NULL)
  
}
