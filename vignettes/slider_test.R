library(slider)

foo <- function(x){
	(x - mean(x, na.rm = TRUE))^2
}

N = 1e3
H = 0.7
x <- seq(1:N)
df <- data.frame(x = x, fbm = as.numeric(somebm::fbm(hurst = H, n = length(x) - 1)))
fit <- lm(fbm ~ x, data = df)
df$original_fbm <- df$fbm
df$fbm <- df$fbm - fit$fitted.values

row <- df$fbm
tictoc::tic()
ACV <- stats::acf(row, plot = FALSE, type = "covariance", demean = FALSE, lag.max = length(row) - 1, na.action = na.pass)
tictoc::toc()

tictoc::tic()
w = sapply(seq(1, N - 1), function(L){
	sqrt(mean(unlist(slider::slide(df$fbm, ~stats::var(.x, na.rm = TRUE), .after = L)), na.rm = TRUE))
})
tictoc::toc()

plot(seq(1, N-1), w, log = "xy")


get_alpha_(w, 1, 200)