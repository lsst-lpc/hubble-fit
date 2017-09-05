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

$> time hubble-fit

Number of supernovae :  740
mean of the residuals  0.014391538550558799
mean of the absolute residuals  0.12726543895212666
hubble: initial parameters: [0.295 0.141 3.101 -19.05 -0.07]
res=&{Location:{X:[0.295 0.141 3.101 -19.05 -0.07] F:652.7606127419591 Gradient:[] Hessian:<nil>} Stats:{MajorIterations:19 FuncEvaluations:41 GradEvaluations:0 HessEvaluations:0 Runtime:24.851189993s} Status:FunctionConvergence}
hubble: status = FunctionConvergence
hubble: func   = 652.760613
==== results ====
Omega M = +0.2950
Alpha   = +0.1410
Beta    = +3.1010
Mb      = -19.0500
Delta M = -0.0700

real 0m26.481s
user 0m27.824s
sys  0m0.496s
```

## Python reference

```sh
$> cd $GOPATH/github.com/lsst-lpc/hubble-fit
$> time python2 ./py/hubblefit.py

./py/hubblefit.py:387: InitialParamWarning: errordef is not given. Default to 1.
  limit_Mb=(-20., -18.), limit_delta_M=(-0.1, -0.0), fix_omgM=False, fix_alpha=False, fix_beta=False, fix_Mb=False, fix_delta_M=False, print_level=1)
**************************************************
*                     MIGRAD                     *
**************************************************

**********************************************************************
---------------------------------------------------------------------------------------
fval = 682.8910357476947 | total call = 90 | ncalls = 90
edm = 5.6711138504257765e-05 (Goal: 1e-05) | up = 1.0
---------------------------------------------------------------------------------------
|          Valid |    Valid Param | Accurate Covar |         Posdef |    Made Posdef |
---------------------------------------------------------------------------------------
|           True |           True |           True |           True |          False |
---------------------------------------------------------------------------------------
|     Hesse Fail |        Has Cov |      Above EDM |                |  Reach calllim |
---------------------------------------------------------------------------------------
|          False |           True |          False |            u'' |          False |
---------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------
|      |  Name   |  Value   | Para Err |   Err-   |   Err+   |  Limit-  |  Limit+  |          |
------------------------------------------------------------------------------------------------
|    0 |    omgM =  0.2951  |  0.03327 |          |          |  0.2     |  0.4     |          |
|    1 |   alpha =  0.1412  |  0.00657 |          |          |  0.1     |  0.2     |          |
|    2 |    beta =  3.102   |  0.08064 |          |          |  2       |  4       |          |
|    3 |      Mb = -19.05   |  0.02317 |          |          | -20      | -18      |          |
|    4 | delta_M = -0.07008 |  0.02216 |          |          | -0.1     | -0       |          |
------------------------------------------------------------------------------------------------

**********************************************************************

real 1m17.719s
user 1m12.372s
sys  0m4.808s
```
