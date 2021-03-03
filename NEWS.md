# statisticalRoughness 0.4


New features

- The following functions were added to compute the height-height correlation function, roughness exponents and anisotropy exponents, and power spectrum decay slope:
	+ `get_hhcf()_`: based on the auto-covariance function
	+ `get_alpha()_`: using the auto-correlation length as break point (instead of a bootstrapped segmented regression)
	+ `get_all_alpha()_`: wrapper around `get_alpha_()`
	+ `get_beta()`: still uses a bootstrapped segmented regression but is only performed once per tile
- The following wrapper functions were added:
	+ `get_zeta_df()`
	+ `get_zeta_raster()`
- The following functions were added for loading and visualizing the results:
	+ `clamp_raster_sigmas()`
	+ `crossscale_correlations()`
	+ `four_values_check()`
	+ `get_kruskal_flag()`
	+ `make_all_plots()`
	+ `make_angular()`
	+ `make_leaflet_map()`
	+ `make_stacked_density_plot()`
	+ `modes_from_stacked_density()`
	+ `raster_select()`
	+ `read_zeta_raster()`
	+ `reduce_spatial_noise()`
	+ `slice_clamp()`


Deprecation

The following functions are not used anymore:
	+ `get_hhcf()`
	+ `get_alpha()`
	+ `get_all_alpha()`


- The following functions were modified

# statisticalRoughness 0.3

New features

- The following functions were added to compute the height-height correlation function:
	+ `get_hhcf()`
- The following functions were added to compute roughness exponents:
	+ `get_alpha()`
	+ `alpha_plot()`
	+ `filter_alpha()`
	+ `summarise_alpha()`
	+ `get_all_alpha()`
- The following functions were added to compute anisotropy exponents:
	+ `get_zeta()`
- Some convenience functions:
	+ `signifNA()` to better handle piping into `signif()` in vignettes

Enhancements

- added one vignette describing the derivation of roughness and anisotropy exponents
- added one vignette comparing derivations of roughness and anisotropy exponents in three constrasting landscapes

# statisticalRoughness 0.2

New features

- The following functions were added to allow a 2D Fourier analysis of topography:
	+ `detrend_dem()`
	+ `pad_mat()`
	+ `Hann2D()`
	+ `fft2D()`
- The following functions were added to rotate raster according to a 2D Fourier analysis:
	+ `bin()`
	+ `get_normalized_spectral_power_matrix()`
	+ `filter_spectral_power_matrix()`
	+ `get_fourier_angle()`
	+ `rotate_raster()`
- Some visualization functions:
	+ `view_matrix()`
	+ `spectrum_plot()`	

Enhancements

- re-organized the Reference and the Articles
- added one vignette describing the basic use of Fourier analysis functions
- added one vignette describing how to rotate rasters along 2DFT main components

# statisticalRoughness 0.1

Initial release.

- Contains the following functions
	+ `all_rasters_to_polygons()`
	+ `compute_Hurst_rasters()`
	+ `compute_Hurst_rasters_internal()`
	+ `get_mean_res()`
	+ `get_mean_sd_rasters()`
	+ `merge_Hurst_rasters()`
	+ `get_all_R_L()`
	+ `get_factors()`
	