#' Statistics on cross-scale correlations
#' @param graphics_df a `data.frame` output from `crossscale_correlations()`
#' @return a list with two elements, the plots and the graphics_df
#' @export
#' @importFrom rlang .data
summarise_crossscale_correlations <- function(graphics_df, stat = "mean", threshold = 0, filter_xy = FALSE){
  if(class(graphics_df) != "data.frame") stop("`graphics_df` is not of class `data.frame`")
  .list = list(min = min, mean = mean, max = max, sd = sd)
  if (filter_xy){
    graphics_df <- graphics_df %>% dplyr::filter(!(grepl("\\.x", var1) | grepl("\\.x", var2) | grepl("\\.y", var1) | grepl("\\.y", var2) | var1 == "zeta1" | var2 == "zeta1" | var1 == "zeta2" | var2 == "zeta2"))
  }
  new_graphics_df <- graphics_df %>% 
    dplyr::mutate(vars = paste0(.data$var1, "-", .data$var2)) %>%
    dplyr::group_by(vars) %>%
    dplyr::summarise(dplyr::across(dplyr::contains("value"), .list)) %>%
    reshape2::melt(id.vars = "vars") %>%
    dplyr::mutate(variable = gsub("value_", "", .data$variable))
    p <- lapply(names(.list), function(n){
      ggpubr::ggdotchart(
      new_graphics_df %>% dplyr::filter(variable == n), 
      x = "vars",
      y = "value",
      add = "segments",
      sorting = "descending",
      rotate = TRUE,
      ggtheme = ggplot2::theme_bw(),
      title = n,
      dot.size = 1)
    })
    names(p) <- names(.list)
    distance_matrix <- new_graphics_df %>% 
      dplyr::mutate(value = abs(.data$value)) %>%
      dplyr::filter(variable == stat & value >= threshold) %>%
      dplyr::mutate(var1 = sapply(vars, function(x) unlist(strsplit(x, "-"))[[1]])) %>%
      dplyr::mutate(var2 = sapply(vars, function(x) unlist(strsplit(x, "-"))[[2]])) %>% 
      dplyr::select(.data$var1, .data$var2, .data$value) %>% 
      dplyr::arrange(.data$var1, .data$var2)
    lvls <- unique(c(distance_matrix$var1, distance_matrix$var2))
    distance_matrix <- distance_matrix %>%
      dplyr::mutate(var1 = factor(.data$var1, levels = lvls)) %>%
      dplyr::mutate(var2 = factor(.data$var2, levels = lvls))
    distance_matrix <- rbind(distance_matrix, distance_matrix %>% dplyr::mutate(var1 = distance_matrix$var2, var2 = distance_matrix$var1)) %>% 
    dplyr::arrange(var2) %>%
    tidyr::spread(.data$var2, .data$value, drop = FALSE)
  distance_matrix$var1 <- NULL
  distance_matrix <- as.matrix(distance_matrix)
  distance_matrix <- distance_matrix[order(colnames(distance_matrix)),order(colnames(distance_matrix))]
  groupColors <- viridisLite::viridis(nrow(distance_matrix), begin = .1, end = .9)
  groupColors <- sapply(groupColors, function(x) substring(x, first = 1, last = 7)) %>% unname()
  cdiag <- chorddiag::chorddiag(distance_matrix, groupColors = groupColors, groupnamePadding = 20, precision = 3)
  return(list(dotplots = p, cdiag = cdiag))
}