# hubble-fit

`hubble-fit` is a collection of functions to fit the Hubble diagram with a Cosmological model.

The initial version was written in [Python](https://python.org).
This is a reimplementation in [Go](https://golang.org).

## Introduction

The `data/jla_lcparams.txt` file contains light-curve parameters as described by the [Joint Light-curve Analysis](http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html) package:

```
name: name of the SN
zcmb: CMB frame redshift (including peculiar velocity corrections for
      nearby supernova based on the models of M.J. Hudson)
zhel: Heliocentric redshift (note both zcmb and zhel are needed
      to compute the luminosity distance)
dz: redshift error (no longer used by the plugin)
mb: B band peak magnitude
dmb: Error in mb (includes contributions from intrinsic dispersion,
     lensing, and redshift uncertainty)
x1: SALT2 shape parameter
dx1: Error in shape parameter
colour: Colour parameter
dcolour: Error in colour
3rdvar: In these files, the log_10 host stellar mass
d3rdvar: Error in 3rdvar
tmax: Date of peak brightness (mjd)
dtmax: Error in tmax
cov_m_s: The covariance between mb and x1
cov_m_c: The covariance between mb and colour
cov_s_c: The covariance between x1 and colour
set: A number indicating which sample this SN belongs to, with
   1 - SNLS, 2 - SDSS, 3 - low-z, 4 - Riess HST
ra: Right Ascension in degree (J2000)
dec: Declination in degree (J2000)
biascor: The correction for analysis bias applied to measured magnitudes
	 (this correction is already applied to mb, original measurements
	  can be obtained by subtracting this term to mb)
```

Note that you will need the covariance matrices from JLA to run the `hubble-fit` code.

## Example

Once you have installed the [Go](https://golang.org) toolchain:

```sh
$> go get github.com/lsst-lpc/hubble-fit
$> cd $GOPATH/github.com/lsst-lpc/hubble-fit
$> hubble-fit
```
