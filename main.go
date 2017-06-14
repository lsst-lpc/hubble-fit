// Copyright 2017 The hubble-fit Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// hubble-fit is a program to fit the Hubble diagram with a Cosmological model.
package main

import (
	"fmt"
	"image/color"
	"io/ioutil"
	"log"
	"math"
	"os"
	"strconv"
	"strings"

	"github.com/astrogo/fitsio"
	"github.com/gonum/diff/fd"
	"github.com/gonum/integrate"
	"github.com/gonum/matrix/mat64"
	"github.com/gonum/optimize"
	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/vg"
	"github.com/gonum/plot/vg/draw"
	"github.com/gonum/stat"
//	"go-hep.org/x/hep/fit"
)

/*
#####
Constants
#####
*/

const (
	c       = 299792458	// speed of light
	H       = 0.070		// Hubble constant
	omgM    = 0.295		// Mass density parameter, = (8*pi*G*rho_mass)/(3*H²)
	alpha   = 0.141		// free parameter for the Hubble fit (factor for stretch)
	beta    = 3.101		// free parameter for the Hubble fit (factor for color)
	Mb      = -19.05	// absolute blue magnitude of a 1A SN
	delta_M = -0.070	// uncertainty on the absolute blue magnitude
)

/*
#####
Basic analysis functions
#####
*/

func Sx(s []float64) float64 {

	// compute the standard deviation of a given distribution

	var sum, sigma float64
	for _, i := range s {
		sum += i * i
	}

	sigma = math.Sqrt(sum / (float64(len(s)) - 1))

	return sigma
}

func RMS(s []float64) float64 {

	// compute the root mean square of a given distribution

	var sum, rms float64
	for _, i := range s {
		sum = +i * i
	}

	rms = math.Sqrt(sum / float64(len(s)))

	return rms
}

func RMSerr(s []float64) float64 {

	// compute the error on the root mean square of a given distribution

	var rmserr float64
	rmserr = Sx(s) / math.Sqrt(2*float64(len(s)))

	return rmserr
}

func MEANerr(s []float64) float64 {

	// compute the error on the mean of a given distribution

	var meanerr float64
	meanerr = Sx(s) / math.Sqrt(float64(len(s)))

	return meanerr
}

func difference(s1, s2 []float64) ([]float64, []float64) {

	// compute the difference between two slices, value by value

	diff := make([]float64, len(s1))
	abs := make([]float64, len(s1))

	for i := range(s1) {
		val := s2[i] - s1[i]
		diff[i] = val
		abs[i] = math.Abs(val)
	}

	return diff, abs
}

func AtoF (s string) float64 {

	// convert a string to float64

	a := strings.Trim(s, " ")
	v, err := strconv.ParseFloat(a, 64)
	if err != nil {
		log.Panic(err)
	}
	return v
}

/*
#####
Theoretical mu
#####
*/

func integral(z, n float64) float64 {

	// compute the integral of the D_l function, based on the trapeze method
	// integrate.Trapezoidal takes two slices (abscissa and ordinate) from which to compute the integral
	// inputs : z = redshift, n = number of steps for the trapeze method

	abs := make([]float64, int(n+1))
	ord := make([]float64, int(n+1))

	for i := 0; i <= int(n); i++ {
		abs[i] = z * float64(i) / n
	}

	for j := 0; j <= int(n); j++ {
		ord[j] = 1 / math.Sqrt((1+abs[j]*omgM)*(1+abs[j])*(1+abs[j])-(1-omgM)*(2+abs[j])*abs[j])
	}

	return integrate.Trapezoidal(abs, ord)

}

func mu_th(z, n float64) float64 {

	// compute the theoretical mu for each SN
	// inputs : z = redshift, n = number of steps for the integral, to be passed to integral function

	return 5 * math.Log10(((1+z)*c/(10*H))*integral(z, n))
}

func mu_th_slice(s []float64, n float64) []float64 {

	// builds a slice containing the theoretical distances
	// inputs : s = slice of the zcmb, n = number of steps for the integral, to be passed to the integral function through mu_th

	diag_ord := make([]float64, len(s))

	for i := 0; i < len(s); i++ {
		diag_ord[i] = mu_th(s[i], n)
	}

	return diag_ord
}

/*
######
Experimental mu
######
*/

func mu_exp_slice(m_stell, mb, stretch, colour []float64) []float64 {

	// compute experimental mus and add them in an array
	// inputs : mB = B band peak magnitude, stretch and colour = SALT2 parameters

	slice := make([]float64, len(m_stell))
	for i := 0; i < len(mb); i++ {
		if m_stell[i] < 10 {
			slice[i] = mb[i] - Mb + alpha*stretch[i] - beta*colour[i]
		} else {
			slice[i] = mb[i] - Mb - delta_M + alpha*stretch[i] - beta*colour[i]
		}
	}

	return slice
}

func mu_exp_err(dmb, dstretch, dcolour []float64) []float64 {

	// compute errors on experimental mus and add them in an array
	// inputs : dmb, dstretch, dcolour = errors on mb, stretch and colour from JLA data

	slice := make([]float64, len(dmb))
	for i := range(dmb) {
		ai := dmb[i]
		bi := alpha*dstretch[i]
		ci := beta*dcolour[i]
		slice[i] = math.Sqrt(ai*ai + bi*bi + ci*ci)
	}

	return slice
}


/*
######
Functions for diagrams : plots, fit, etc
######
*/

func newPoints(xs, ys []float64) plotter.XYs {

	// creates a plotter.XYs (cloud of points) from given 2D coordinates
	// inputs : xs, ys = slices containing the x and y values respectively

	pts := make(plotter.XYs, len(xs))
	for i := range pts {
		pts[i].X = xs[i]
		pts[i].Y = ys[i]
	}

	return pts
}

func newHubblePlot(points plotter.XYs) *plot.Plot {

	// plot the given points
	// input : points = plotter.XYs created by the points function

	p, err := plot.New()
	if err != nil {
		log.Panic(err)
	}

	p.Title.Text = "Hubble diagram"
	p.X.Label.Text = "Redshift z"
	p.Y.Label.Text = "Distance mu"
	p.X.Scale = plot.LogScale{}
	p.X.Tick.Marker = plot.LogTicks{}
	p.Add(plotter.NewGrid())

	s, err := plotter.NewScatter(points)
	if err != nil {
		log.Panic(err)
	}

	s.GlyphStyle.Color = color.RGBA{R: 255, B: 0, A: 255}
	s.GlyphStyle.Radius = vg.Points(3)
	s.GlyphStyle.Shape = draw.PlusGlyph{}

	p.Add(s)

	return p
}

func modified_integral(z, omega float64) float64 {

	// modified integral function where omgM is a variablef that can be optimized
	// inputs : z = redshift, omega = omgM as a free parameter

	const n = 1000
	xs := make([]float64, n+1)
	ys := make([]float64, n+1)

	for i := range(xs) {
		x := z * float64(i) / float64(n)
		xs[i] = x
		ys[i] = 1 / math.Sqrt((1+x*omega)*(1+x)*(1+x)-(1-omega)*(2+x)*x)
	}

	return integrate.Trapezoidal(xs, ys)
}

func newHistogram(s []float64) *plot.Plot {

	// builds the histogram
	// input : s = slice containing the values to plot

	p, err := plot.New()
	if err != nil {
		log.Panic(err)
	}

	val := make(plotter.Values, len(s))
	for i := range s {
		val[i] = s[i]
	}

	p.Title.Text = "Histogram"
	h, err := plotter.NewHist(val, len(val))
	if err != nil {
		log.Panic(err)
	}
	p.Add(h)

	return p
}

/*
######
Matrices and chi2
######
*/

func Chi2(ps []float64) float64 {

	// reads the raw datas and format it into a slice of strings, each containing the information for a SN
	datas, err := ioutil.ReadFile("./data/jla_lcparams.txt")
	if err != nil {
		fmt.Println(err)
	}
	str := strings.Trim(string(datas), "\n")		// gets rid, if needed, of a new line symbol at the end of the file
	supernovae := strings.Split(str, "\n")			// builds a slice with data for a SN in each item
	supernovae = append(supernovae[:0], supernovae[1:]...)	// suppress the first item (labels of the columns)

	// creates slices with the needed values from JLA datas

	N := len(supernovae)
	zcmb := make([]float64, N)	// redshift in the cmb frame
	mb := make([]float64, N)	// b band peak magnitude
	stretch := make([]float64, N)	// stretch factor
	colour := make([]float64, N)	// colour factor
	m_stell := make([]float64, N)	// log10 of the host galaxy stellar mass

	for i, v := range supernovae {
		split_str := strings.Split(v, " ")
		if len(split_str) != 21 {
			os.Exit(1)
		}

		zcmb[i]	= AtoF(split_str[1])
		mb[i] = AtoF(split_str[4])
		stretch[i] = AtoF(split_str[6])
		colour[i] = AtoF(split_str[8])
		m_stell[i] = AtoF(split_str[10])
	}

	mu_exp, mu_th, mu_diff := make([]float64, N), make([]float64, N), make([]float64, N)

	for i := 0; i < N; i++ {
		if m_stell[i] < 10 {
			mu_exp[i] = mb[i] - ps[3] + ps[1]*stretch[i] - ps[2]*colour[i]
		} else {
			mu_exp[i] = mb[i] - ps[3] - ps[4] + ps[1]*stretch[i] - ps[2]*colour[i]
		}
	}

	var n int = 1000
	for i := 0; i < N; i++ {
		abs := make([]float64, n+1)
		ord := make([]float64, n+1)

		for j := 0; j <= n; j++ {
			abs[j] = zcmb[i] * float64(j) / float64(n)
		}

		for k := 0; k <= n; k++ {
			ord[k] = 1 / math.Sqrt((1+abs[k]*ps[0])*(1+abs[k])*(1+abs[k])-(1-ps[0])*(2+abs[k])*abs[k])
		}

		mu_th[i] = 5 * math.Log10(((1+zcmb[i])*c/(10*H))*integrate.Trapezoidal(abs,ord))
	}

	for i := 0; i < N; i++ {
		mu_diff[i] = mu_exp[i] - mu_th[i]
	}

	c_eta := mat64.NewDense(2220, 2220, nil)
	c_eta_temp := mat64.NewDense(2220, 2220, nil)

	c_eta_temp = matrice("./covmat/C_bias.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_cal.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_dust.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_host.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_model.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_nonia.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_pecvel.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_stat.fits")
	c_eta.Add(c_eta, c_eta_temp)

	a_vector := mat64.NewVector(3, []float64{1, ps[1], -ps[2]})

	c_mat_elements := make([]float64, N*N)
	i, j := 0, 0
	for k := 0; k < N*N; k++ {
		submat := c_eta.Slice(i, i+3, j, j+3)
		element := mat64.Inner(a_vector, submat, a_vector)
		c_mat_elements[k] = element
		if i+3 < N*3-1 {
			i = i+3
		} else {
			j = j+3
			i = 0
		}
	}

	c_mat := mat64.NewDense(N, N, c_mat_elements)

	data, err := ioutil.ReadFile("./covmat/sigma_mu.txt")
	if err != nil {
		fmt.Println(err)
	}
	str2 := strings.Trim(string(data), " \n")
	elements := strings.Split(str2, "\n")
	elements = append(elements[:0], elements[4:]...)

	sigma_mat := mat64.NewDense(N, N, nil)
	for i :=0; i<N; i++ {
		words := strings.Split(elements[i], " ")
		sz, slen, scoh := AtoF(words[4]), AtoF(words[2]), AtoF(words[0])
		val := ((5*150000)/(sz*c))*((5*150000)/(sz*c)) + slen*slen + scoh*scoh
		sigma_mat.Set(i, i, val)
	}

	covmat := mat64.NewDense(N, N, nil)
	covmat.Add(c_mat, sigma_mat)

	mu_vector := mat64.NewVector(len(mu_diff), mu_diff)
	rows, cols := covmat.Dims()
	if rows != cols {
		log.Fatalf("cov-matrix not square")
	}
	inv := mat64.NewSymDense(rows, nil)
	for i := 0; i < rows; i++ {
		for j := i; j< rows; j++ {
			inv.SetSym(i, j, covmat.At(i, j))
		}
	}
	var chol mat64.Cholesky
	if ok := chol.Factorize(inv); !ok {
		log.Fatalf("cov-matrix not positive semi-definite")
	}

	err2 := inv.InverseCholesky(&chol)
	if err2 != nil {
		log.Fatal(err2)
	}

	return mat64.Inner(mu_vector, inv, mu_vector)
}

/*

// reads a FITS fil, extracts the data and convert them into the corresponding matrix
func matrice (file string) *mat64.Dense {

	r, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer r.Close()

	fits, err := fitsio.Open(r)
	if err != nil {
		panic(err)
	}
	defer fits.Close()

	hdu := fits.HDU(0)
	img := hdu.(fitsio.Image)
	hdr := img.Header()
	rows := hdr.Axes()[0]
	cols := hdr.Axes()[1]
	raw := make([]float64, rows*cols)
	err = img.Read(&raw)

	Mat := mat64.NewDense(2220, 2220, raw)

	return Mat
}

// builds the covariance matrix (equation 13)
func covar_matrix (file string, alpha, beta float64, N int, c_eta *mat64.Dense) *mat64.Dense {

	a_vector := mat64.NewVector(3, []float64{1, alpha, -beta})

	c_mat_elements := make([]float64, N*N)
	i, j := 0, 0
	for k := 0; k < N*N; k++ {
		submat := c_eta.Slice(i, i+3, j, j+3)
		element := mat64.Inner(a_vector, submat, a_vector)
		c_mat_elements[k] = element
		if i+3 < N*3-1 {
			i = i+3
		} else {
			j = j+3
			i = 0
		}
	}

	c_mat := mat64.NewDense(N, N, c_mat_elements)

	data, err := ioutil.ReadFile(file)
	if err != nil {
		fmt.Println(err)
	}
	str := strings.Trim(string(data), " \n")
	elements := strings.Split(str, "\n")
	elements = append(elements[:0], elements[4:]...)

	sigma_mat := mat64.NewDense(N, N, nil)
	for i :=0; i<N; i++ {
		words := strings.Split(elements[i], " ")
		sz, slen, scoh := AtoF(words[4]), AtoF(words[2]), AtoF(words[0])
		val := ((5*150000)/(sz*c))*((5*150000)/(sz*c)) + slen*slen + scoh*scoh
		sigma_mat.Set(i, i, val)
	}

	covar_mat := mat64.NewDense(N, N, nil)
	covar_mat.Add(c_mat, sigma_mat)

	return covar_mat
}

func chi2 (diff_data []float64, covmat *mat64.Dense) float64 {

	mu_vector := mat64.NewVector(len(diff_data), diff_data)
	rows, cols := covmat.Dims()
	if rows != cols {
		log.Fatalf("cov-matrix not square")
	}
	inv := mat64.NewSymDense(rows, nil)
	for i := 0; i < rows; i++ {
		for j := i; j< rows; j++ {
			inv.SetSym(i, j, covmat.At(i, j))
		}
	}
	var chol mat64.Cholesky
	if ok := chol.Factorize(inv); !ok {
		log.Fatalf("cov-matrix not positive semi-definite")
	}

	err := inv.InverseCholesky(&chol)
	if err != nil {
		log.Fatal(err)
	}

	return mat64.Inner(mu_vector, inv, mu_vector)
}
*/
// Curve1D returns the result of a non-linear least squares to fit
// a function f to the underlying data with method m.
func Curve1D(f Func1D, settings *optimize.Settings, m optimize.Method) (*optimize.Result, error) {
	f.init()

	p := optimize.Problem{
		Func: f.fct,
		Grad: f.grad,
	}

	if m == nil {
		m = &optimize.NelderMead{}
	}

	p0 := make([]float64, len(f.Ps))
	copy(p0, f.Ps)
	return optimize.Local(p, p0, settings, m)
}

//Func1D describes a 1D function to fit some data
// this is a modified version of the fit function from the go-hep package
type Func1D struct {
	// F is the function to minimize.
	// ps is the slice of parameters to optimize during the fit.
	F func(x float64, ps []float64) float64

	// N is the number of parameters to optimize during the fit.
	// If N is 0, Ps must not be nil.
	N int

	// Ps is the initial values for the parameters.
	// If Ps is nil, the set of initial parameters values is a slice of
	// length N filled with zeros.
	Ps []float64

	X   []float64
	Y   []float64
	Err []float64

	sig2 []float64 // inverse of squares of measurement errors along Y.

	fct  func(ps []float64) float64 // cost function (objective function)
	grad func(grad, ps []float64)
}

func (f *Func1D) init() {

	f.sig2 = make([]float64, len(f.Y))
	switch {
	default:
		for i := range f.Y {
			f.sig2[i] = 1
		}
	case f.Err != nil:
		for i, v := range f.Err {
			f.sig2[i] = 1 / (v * v)
		}
	}

	if f.Ps == nil {
		f.Ps = make([]float64, f.N)
	}

	if len(f.Ps) == 0 {
		panic("fit: invalid number of initial parameters")
	}

	if len(f.X) != len(f.Y) {
		panic("fit: mismatch length")
	}

	if len(f.sig2) != len(f.Y) {
		panic("fit: mismatch length")
	}

	f.fct = func(ps []float64) float64 {
		var chi2 float64
		for i := range f.X {
			res := f.F(f.X[i], ps) - f.Y[i]
			chi2 += res * res * f.sig2[i]
		}
		return 0.5 * chi2
	}

	f.grad = func(grad, ps []float64) {
		fd.Gradient(grad, f.fct, ps, nil)
	}
}

/*
######
Main function
######
*/

func main() {

	// reads the raw datas and format it into a slice of strings, each containing the information for a SN
	datas, err := ioutil.ReadFile("./data/jla_lcparams.txt")
	if err != nil {
		fmt.Println(err)
	}
	str := strings.Trim(string(datas), "\n")		// gets rid, if needed, of a new line symbol at the end of the file
	supernovae := strings.Split(str, "\n")			// builds a slice with data for a SN in each item
	supernovae = append(supernovae[:0], supernovae[1:]...)	// suppress the first item (labels of the columns)

	// creates slices with the needed values from JLA datas

	N := len(supernovae)
	sn_names := make([]string, N)	// names
	zcmb := make([]float64, N)	// redshift in the cmb frame
	zhel := make([]float64, N)	// redshift in the heliocentric frame
	mb := make([]float64, N)	// b band peak magnitude
	dmb := make([]float64, N)	// b band peak magnitude error
	stretch := make([]float64, N)	// stretch factor
	dstretch := make([]float64, N)	// stretch factor error
	colour := make([]float64, N)	// colour factor
	dcolour := make([]float64, N)	// colour factor error
	m_stell := make([]float64, N)	// log10 of the host galaxy stellar mass
	exp := make([]float64, N)	// supernova sample identifier : 1 = SNLS, 2 = SDSS, 3 = lowz, 4 = Riess HST
	ra_jla := make([]float64, N)	// right ascension (in degrees)
	de_jla := make([]float64, N)	// declination (in degrees)

	for i, v := range supernovae {
		split_str := strings.Split(v, " ")
		if len(split_str) != 21 {
			os.Exit(1)
		}

		sn_names[i] = split_str[0]
		zcmb[i]	= AtoF(split_str[1])
		zhel[i] = AtoF(split_str[2])
		mb[i] = AtoF(split_str[4])
		dmb[i] = AtoF(split_str[5])
		stretch[i] = AtoF(split_str[6])
		dstretch[i] = AtoF(split_str[7])
		colour[i] = AtoF(split_str[8])
		dcolour[i] = AtoF(split_str[9])
		m_stell[i] = AtoF(split_str[10])
		exp[i] = AtoF(split_str[17])
		ra_jla[i] = AtoF(split_str[18])
		de_jla[i] = AtoF(split_str[19])
	}

	fmt.Println("Number of supernovae : ", len(sn_names))

	// compute and store slices needed for the plots and analysis

	muexp_data := mu_exp_slice(m_stell, mb, stretch, colour)	// experimental mus from SALT2 model
	muth_data := mu_th_slice(zcmb, 10000)				// theoretical mus from lambdaCDM model
	diff_data, abs_diff_data := difference(muth_data, muexp_data)	// difference and absolute difference between exp and theoretical mus
	mean_residuals := stat.Mean(diff_data, nil)			// mean of the differences
	abs_mean_residuals := stat.Mean(abs_diff_data, nil)		// mean of the absolute differences
//	muexp_error := mu_exp_err(dmb, dstretch, dcolour)		// errors on mu_exp, from error propagation in SALT2 model

	fmt.Println("mean of the residuals ", mean_residuals)
	fmt.Println("mean of the absolute residuals ", abs_mean_residuals)

	// compute, draw and store the various diagrams

	hubblePlot := newHubblePlot(newPoints(zcmb, muexp_data))
	err = hubblePlot.Save(1080, 720, "hubble_diagram.png")
	if err != nil {
		log.Fatalf("Error saving 'hubble_diagram.png' : %v", err)
	}

	hubblePlot2 := newHubblePlot(newPoints(zcmb, diff_data))
	err = hubblePlot2.Save(720, 200, "hubble_diagram2.png")
	if err != nil {
		log.Fatalf("Error saving 'hubble_diagram2.png' : %v", err)
	}

	histogram := newHistogram(diff_data)
	err = histogram.Save(500, 500, "histogram.png")
	if err != nil {
		log.Fatalf("Error saving 'histogram.png' : %v", err)
	}

/*	// reads the FITS files and build the matrix C_eta

	c_eta := mat64.NewDense(2220, 2220, nil)
	c_eta_temp := mat64.NewDense(2220, 2220, nil)

	c_eta_temp = matrice("./covmat/C_bias.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_cal.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_dust.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_host.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_model.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_nonia.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_pecvel.fits")
	c_eta.Add(c_eta, c_eta_temp)

	c_eta_temp = matrice("./covmat/C_stat.fits")
	c_eta.Add(c_eta, c_eta_temp)


	covariance_matrix := covar_matrix("./covmat/sigma_mu.txt", alpha, beta, N, c_eta)

	chi2 := chi2(diff_data, covariance_matrix)
	fmt.Println("Chi2 : ", chi2)

	empty := make([]float64, N)
	params := []float64{0.295, 0.141, 3.101, -19.05, -0.70}
	res, err := Curve1D(
		Func1D{
			F: func(z float64, ps []float64) float64 {
				var i int = -1
				for p, v := range zcmb {
					if v == z {
						i = p
					}
				}
				if i < 0 {
					panic("error in retrieving index value in the fitting function")
				}

				out := mb[i]-ps[3]+ps[1]*stretch[i]-ps[2]*colour[i] - 5*math.Log10(((1+z)*c/(10*H))*modified_integral(z, ps[0]))
				if m_stell[i] >= 10 {
					out -= ps[4]
				}
				return out
			},
			X: zcmb,
			Y: empty,
			Err: muexp_error,
			Ps: params,
		},
		nil, &optimize.NelderMead{},
	)
	if err != nil {
		log.Fatalf("Error calling Curve1D : %v", err)
	}

	fmt.Println("Omega M : ", res.X[0])
	fmt.Println("Alpha : ", res.X[1])
	fmt.Println("Beta : ", res.X[2])
	fmt.Println("Mb : ", res.X[3])
	fmt.Println("Delta M : ", res.X[4])
*/
	var params = []float64{0.295, 0.141, 3.101, -19.05, -0.070}
	temp := Chi2(params)
	fmt.Println(temp)
}
