#' Make a graph for summary parameter estimates by ggplot2
#'
#' @param Estimate The name of the data.frame that contains summary parameter estimates
#' @param type Which type of the summary parameter is of interest
#' @param e1 The evaluation point e
#' @param bstrap
#' Boolean for whether or not to perform the multiplier bootstrap inference.
#'
#' @returns A ggplot2 figure for the summary parameter estimates and the corresponding uniform confidence bands
#'
#' @importFrom ggplot2 aes element_line element_rect element_text geom_line geom_ribbon ggplot labs theme
#' @importFrom purrr %>%
#'
#' @noRd
#'
graph_func_aggte_continuous <- function(Estimate, type, e1, bstrap) {
  
  mytitle <- paste0("Evaluation point ", e1, ".")
  
  # mycolor1 <- mycolor2 <- "#00BFC4"
  
  mycolor1 <- "black"
  mycolor2 <- "grey"
  
  if (!bstrap) {
    
    Estimate <- Estimate %>%
      mutate(ci_lower = ci1_lower,
             ci_upper = ci1_upper)
    
  } else {
    
    Estimate <- Estimate %>%
      mutate(ci_lower = ci2_lower,
             ci_upper = ci2_upper)
    
  }
  
  if (type %in% c("dynamic", "group", "calendar")) {
    
    Estimate <- Estimate %>%
      filter(eval == e1)
    
  }
  
  graph <- Estimate %>%
    ggplot(aes(x = z, group = 1)) +
    geom_line(aes(y = est), color = mycolor1, linewidth = 2) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = mycolor2, alpha = 0.6) +
    theme(
      panel.grid.major = element_line(color = "grey"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0)
    ) +
    labs(title = mytitle, x = expression(z), y = "", parse = TRUE, size = 1)
  
  return(graph)
  
}
