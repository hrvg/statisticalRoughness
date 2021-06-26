#' Fast auto-covariance with fftw
#' @param data_vector a vector of centered data
#' @importFrom fftwtools fftw
#' @return a vector of auto-covariance
acv_fft <- function(data_vector) {
  # Inspired by this https://gist.github.com/FHedin/05d4d6d74e67922dfad88038b04f621c

  #  need to pad with zeroes first ; pad to a power of 2 will give faster FFT
  len_dat <- length(data_vector)
  len_opt <- 2^(ceiling(log2(len_dat))) - len_dat
  fft_dat <- c(data_vector, rep.int(0, len_opt))

  # fft using fast fftw library as backend
  fft_dat <- fftw(fft_dat)

  # take the inverse transform of the power spectral density
  fft_dat <- (abs(fft_dat))^2
  fft_dat <- fftw(fft_dat, inverse = 1)

  # We repeat the same process (except for centering) on a ‘mask_dat’ signal,
  # in order to estimate the error made on the previous computation.
  mask_dat <- rep.int(1, len_dat)
  mask_dat <- c(mask_dat, rep.int(0, len_opt))
  mask_dat <- fftw(mask_dat)
  mask_dat <- fftw((abs(mask_dat))^2, inverse = 1)

  # The “error” made can now be easily corrected by an element-wise division
  fft_dat <- fft_dat / mask_dat

  # keep real parts only (although there should be not imaginary part or really small ones) and remove padding data
  fft_dat <- Re(fft_dat[1:len_dat])

  return(fft_dat)
}
