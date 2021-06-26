#' Visualizes values of a matrix
#' @param m a `matrix`
#' @param ply `logical`
#' @return if `ply = TRUE` returns a `plotly::ggploty()` object, if `FALSE` returns a `ggplot` object
#' @importFrom rlang .data
#' @export
#' @keywords fft2d_viz
view_matrix <- function(m, ply = TRUE) {
  m <- m %>% t()
  m <- reshape2::melt(m)
  p <- ggplot2::ggplot(m, ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(na.value = "transparent") +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = "columns", y = "rows") +
    ggplot2::theme(legend.key.width = ggplot2::unit(1, "in")) +
    ggpubr::theme_pubr() +
    ggplot2::coord_fixed(ratio = 1)
  if (ply) {
    return(plotly::ggplotly(p))
  } else {
    return(p)
  }
}

#' Transforms matrix for be used by `rayshader`.
#' @param m a `matrix`
#' @return a `matrix`
#' @export
#' @keywords fft2d_viz
matrix_for_rayshader <- function(m) {
  m <- t(m)
  return(m)
}

#' Plot a radial Fourier spectrum
#' @param binned_power_spectrum `matrix`, a binned spectrum from `bin()`
#' @param FT2D a `list` from `fft2d()`
#' @param xdecades `numeric`, how decades should be plotted in the x-direction
#' @param ydecades `numeric`, how decades should be plotted in the y-direction
#' @return a `ggplot` object
#' @importFrom rlang .data
#' @export
#' @keywords fft2d_viz
spectrum_plot <- function(binned_power_spectrum, FT2D, xdecades = 3, ydecades = 3) {
  .x <- NULL
  df_bin <- data.frame(frequency = 10^binned_power_spectrum[, 1], normalized_power = 10^binned_power_spectrum[, 2])
  df_complete <- data.frame(frequency = FT2D$radial_frequency_vector, normalized_power = FT2D$spectral_power_vector)
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df_complete, ggplot2::aes(x = .data$frequency, y = .data$normalized_power, alpha = 0.1)) +
    ggplot2::geom_point(data = df_bin, ggplot2::aes(x = .data$frequency, y = .data$normalized_power, color = "red")) +
    ggplot2::labs(title = "Topography power spectrum") +
    ggplot2::labs(y = bquote("DFT mean squared amplitude " ~ (m^2)), x = bquote("Wavenumber " ~ (m^{
      -1
    }))) +
    ggplot2::scale_x_log10(
      breaks = scales::trans_breaks(n = xdecades, "log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      sec.axis = ggplot2::sec_axis(~ 1 / .,
        breaks = scales::trans_breaks(n = xdecades, "log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        name = "Wavelength (m)"
      )
    ) +
    ggplot2::scale_y_log10(
      breaks = scales::trans_breaks(n = ydecades, "log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggpubr::theme_pubr() +
    ggplot2::guides(color = FALSE, alpha = FALSE)
  return(p1)
}
