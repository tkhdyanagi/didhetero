#' Make a graph of continuous CATT estimates by ggplot2
#'
#' @param Estimate The name of the data.frame that contains CATT estimates
#' @param g1 The group being treated at time g
#' @param t1 The evaluated time period
#' @param bstrap
#' Boolean for whether or not to perform the multiplier bootstrap inference.
#'
#' @returns A ggplot2 figure for the CATT estimates and the corresponding uniform confidence bands
#'
#' @importFrom ggplot2 aes element_line element_rect element_text geom_line geom_ribbon ggplot labs theme
#' @importFrom purrr %>%
#'
#' @noRd
#'
graph_func_continuous <- function(Estimate, g1, t1, bstrap) {

  mytitle <- paste0("Group ", g1, ". Time ", t1, ".")

  mycolor <- "#00BFC4"

  if (!bstrap) {

    Estimate <- Estimate %>%
      mutate(ci_lower = ci1_lower,
             ci_upper = ci1_upper)

  } else {

    Estimate <- Estimate %>%
      mutate(ci_lower = ci2_lower,
             ci_upper = ci2_upper)

  }

  graph <- Estimate %>%
    filter(g == g1 & t == t1) %>%
    ggplot(aes(x = z, group = 1)) +
    geom_line(aes(y = est), color = mycolor, linewidth = 1) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = mycolor, alpha = 0.2) +
    theme(
      panel.grid.major = element_line(color = "grey"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0)
    ) +
    labs(title = mytitle, x = expression(z), y = "CATT", parse = TRUE, size = 1)

  return(graph)

}
