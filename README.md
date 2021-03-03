# `statisticalRoughness`

<!-- [![Build Status](https://travis-ci.com/hrvg/statisticalRoughness.svg?token=Dx1gYTrTiuxgW9Sq3s3q&branch=master)](https://travis-ci.com/hrvg/statisticalRoughness) -->
[![CircleCI](https://circleci.com/gh/hrvg/statisticalRoughness.svg?style=svg)](https://circleci.com/gh/hrvg/statisticalRoughness)
[![codecov](https://codecov.io/gh/hrvg/statisticalRoughness/branch/master/graph/badge.svg?token=8VGWGTXUWG)](https://codecov.io/gh/hrvg/statisticalRoughness)

## Purpose

`statisticalRoughness` is an `R` package designed to efficiently compute roughness exponents for large areas and over a range of specified scales.

## How to install

```
# Install development version from GitHub
devtools::install_github("hrvg/statisticalRoughness")
```

It is also likely you'll need to install the latest `terra` version which has dependencies on GDAL (>= 3.0.0), GEOS (>= 3.3.0) and Proj.4 (>= 6.0.0).

On Linux, you can do:

```
sudo add-apt-repository -y ppa:ubuntugis/ubuntugis-unstable
sudo apt-get -q update
sudo apt-get -y install libgdal-dev libgeos-dev libproj-dev 
```

The `terra` [github page](https://github.com/rspatial/terra) has some insights for other operating systems.

## How to use

Check out the [package's vignette to get started](articles/getting-started.html).