#' Make a graph of discrete CATT estimates by ggplot2
#'
#' @param Estimate The name of the data.frame that contains CATT estimates
#' @param g1 The group being treated at time g
#' @param t1 The evaluated time period
#'
#' @returns A ggplot2 figure for the CATT estimates and the corresponding uniform confidence bands
#'
#' @importFrom ggplot2 aes element_line element_rect element_text geom_errorbar geom_hline geom_point ggplot labs position_dodge theme
#' @importFrom purrr %>%
#'
#' @noRd
#'
graph_func_discrete <- function(Estimate, g1, t1) {

  mytitle <- paste0("Group ", g1, ". Time ", t1, ".")

  mycolor <- "#00BFC4"

  graph <- Estimate %>%
    filter(g == g1 & t == t1) %>%
    ggplot(aes(x = z, y = est, group = 1)) +
    geom_errorbar(
      aes(
        ymin = est - est + ci_lower,
        ymax = est - est + ci_upper
      ),
      position = position_dodge(0.5),
      width = 0.1,
      linewidth = 2,
      color = mycolor
    ) +
    geom_point(
      position = position_dodge(0.5),
      size = 6,
      color = mycolor
    ) +
    theme(
      panel.grid.major = element_line(color = "grey"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(color = "black", fill = NA),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0),
      legend.position = "none"
    ) +
    labs(title = mytitle, x = expression(z), y = "CATT", parse = TRUE, size = 1) +
    geom_hline(
      yintercept = 0,
      color = "black",
      linetype = "dashed",
      linewidth = 1.5
    )

  return(graph)

}
