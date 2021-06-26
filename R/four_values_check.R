#' Check that at least four values of the raster stack are not `NA` for each element of `raster_list`
#' @param raster_list a `list` of `stars` objects
#' @param selected `character`, the selected attributes
#' @export
#' @return `logical`
#' @keywords postprocessing
four_values_check <- function(raster_list, selected) {
  sapply(raster_list, function(s) {
    # notNA <- sapply(seq(dim(s)["band"]), function(i) sum(!is.na(c(s[,,,i][[1]]))))
    # all(notNA > 4)
    df <- s %>%
      stars::st_set_dimensions("band", values = seq_along(selected)) %>%
      as.data.frame() %>%
      dplyr::select(-c("x", "y")) %>%
      dplyr::group_by(.data$band) %>%
      dplyr::group_map(~.)
    df <- do.call(cbind, df)
    min(apply(df, MARGIN = 2, FUN = function(x) sum(!is.na(x)))) > 4
  })
}
