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

	//	"go-hep.org/x/hep/fit"

	"gonum.org/v1/gonum/diff/fd"
	"gonum.org/v1/gonum/integrate"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/optimize"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

const (
	c       = 299792458 // speed of light
	H       = 0.070     // Hubble constant
	omgM    = 0.295     // Mass density parameter, = (8*pi*G*rho_mass)/(3*H²)
	alpha   = 0.141     // free parameter for the Hubble fit (factor for stretch)
	beta    = 3.101     // free parameter for the Hubble fit (factor for color)
	Mb      = -19.05    // absolute blue magnitude of a 1A SN
	delta_M = -0.070    // uncertainty on the absolute blue magnitude
)

func main() {

	// reads the raw datas and format it into a slice of strings, each containing the information for a SN
	datas, err := ioutil.ReadFile("./data/jla_lcparams.txt")
	if err != nil {
		fmt.Println(err)
	}
	str := strings.Trim(string(datas), "\n")               // gets rid, if needed, of a new line symbol at the end of the file
	supernovae := strings.Split(str, "\n")                 // builds a slice with data for a SN in each item
	supernovae = append(supernovae[:0], supernovae[1:]...) // suppress the first item (labels of the columns)

	// creates slices with the needed values from JLA datas

	var (
		N        = len(supernovae)
		sn_names = make([]string, N)  // names
		zcmb     = make([]float64, N) // redshift in the cmb frame
		zhel     = make([]float64, N) // redshift in the heliocentric frame
		mb       = make([]float64, N) // b band peak magnitude
		dmb      = make([]float64, N) // b band peak magnitude error
		stretch  = make([]float64, N) // stretch factor
		dstretch = make([]float64, N) // stretch factor error
		colour   = make([]float64, N) // colour factor
		dcolour  = make([]float64, N) // colour factor error
		m_stell  = make([]float64, N) // log10 of the host galaxy stellar mass
		exp      = make([]float64, N) // supernova sample identifier : 1 = SNLS, 2 = SDSS, 3 = lowz, 4 = Riess HST
		ra_jla   = make([]float64, N) // right ascension (in degrees)
		de_jla   = make([]float64, N) // declination (in degrees)
	)

	for i, v := range supernovae {
		split_str := strings.Split(v, " ")
		if len(split_str) != 21 {
			log.Fatalf("invalid number of jla_lcparams fields. got=%d, want=%d", len(split_str), 21)
		}

		sn_names[i] = split_str[0]
		zcmb[i] = atof(split_str[1])
		zhel[i] = atof(split_str[2])
		mb[i] = atof(split_str[4])
		dmb[i] = atof(split_str[5])
		stretch[i] = atof(split_str[6])
		dstretch[i] = atof(split_str[7])
		colour[i] = atof(split_str[8])
		dcolour[i] = atof(split_str[9])
		m_stell[i] = atof(split_str[10])
		exp[i] = atof(split_str[17])
		ra_jla[i] = atof(split_str[18])
		de_jla[i] = atof(split_str[19])
	}

	fmt.Println("Number of supernovae : ", len(sn_names))

	// compute and store slices needed for the plots and analysis

	muexp_data := mu_exp_slice(m_stell, mb, stretch, colour)      // experimental mus from SALT2 model
	muth_data := mu_th_slice(zcmb, 10000)                         // theoretical mus from lambdaCDM model
	diff_data, abs_diff_data := difference(muth_data, muexp_data) // difference and absolute difference between exp and theoretical mus
	mean_residuals := stat.Mean(diff_data, nil)                   // mean of the differences
	abs_mean_residuals := stat.Mean(abs_diff_data, nil)           // mean of the absolute differences
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

	// computes the best values for the five parameters by minimizing chi2
	var params = []float64{0.295, 0.141, 3.101, -19.05, -0.070}
	res, err := FitChi2(Chi2, params, nil, nil)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("res=%+v\n", res)
	fmt.Println("Omega M : ", res.X[0])
	fmt.Println("Alpha : ", res.X[1])
	fmt.Println("Beta : ", res.X[2])
	fmt.Println("Mb : ", res.X[3])
	fmt.Println("Delta M : ", res.X[4])
}

// FitChi2 is the test function for the computation of chi2
func FitChi2(f func(ps []float64) float64, ps []float64, settings *optimize.Settings, m optimize.Method) (*optimize.Result, error) {
	grad := func(grad, ps []float64) {
		fd.Gradient(grad, f, ps, nil)
	}

	if m == nil {
		m = &optimize.NelderMead{}
	}

	p0 := make([]float64, len(ps))
	copy(p0, ps)

	p := optimize.Problem{
		Func: f,
		Grad: grad,
	}

	return optimize.Local(p, p0, settings, m)
}

// Chi2 is the test function for a chi2 computation that takes into account the covariance matrix
func Chi2(ps []float64) float64 {

	// reads the raw datas and format it into a slice of strings, each containing the information for a SN

	datas, err := ioutil.ReadFile("./data/jla_lcparams.txt")
	if err != nil {
		log.Fatal(err)
	}
	str := strings.Trim(string(datas), "\n")               // gets rid, if needed, of a new line symbol at the end of the file
	supernovae := strings.Split(str, "\n")                 // builds a slice with data for a SN in each item
	supernovae = append(supernovae[:0], supernovae[1:]...) // suppress the first item (labels of the columns)

	// creates slices with the needed values from JLA datas

	var (
		N       = len(supernovae)
		zcmb    = make([]float64, N) // redshift in the cmb frame
		mb      = make([]float64, N) // b band peak magnitude
		stretch = make([]float64, N) // stretch factor
		colour  = make([]float64, N) // colour factor
		m_stell = make([]float64, N) // log10 of the host galaxy stellar mass
	)

	for i, v := range supernovae {
		split_str := strings.Split(v, " ")
		if len(split_str) != 21 {
			os.Exit(1)
		}

		zcmb[i] = atof(split_str[1])
		mb[i] = atof(split_str[4])
		stretch[i] = atof(split_str[6])
		colour[i] = atof(split_str[8])
		m_stell[i] = atof(split_str[10])
	}

	// creates slices for theoretical and experimental mus and the associated errors

	var (
		mu_exp  = make([]float64, N)
		mu_th   = make([]float64, N)
		mu_diff = make([]float64, N)
	)

	for i := 0; i < N; i++ {
		if m_stell[i] < 10 {
			mu_exp[i] = mb[i] - ps[3] + ps[1]*stretch[i] - ps[2]*colour[i]
		} else {
			mu_exp[i] = mb[i] - ps[3] - ps[4] + ps[1]*stretch[i] - ps[2]*colour[i]
		}
	}

	for i := 0; i < N; i++ {
		const n = 1000
		xs := make([]float64, n+1)
		ys := make([]float64, n+1)

		for j := range xs {
			xs[j] = zcmb[i] * float64(j) / float64(n)
			ys[j] = 1 / math.Sqrt((1+xs[j]*ps[0])*(1+xs[j])*(1+xs[j])-(1-ps[0])*(2+xs[j])*xs[j])
		}

		mu_th[i] = 5 * math.Log10(((1+zcmb[i])*c/(10*H))*integrate.Trapezoidal(xs, ys))
	}

	for i := 0; i < N; i++ {
		mu_diff[i] = mu_exp[i] - mu_th[i]
	}

	// reads the FITS files and builds C_eta

	c_eta := mat.NewDense(2220, 2220, nil)
	c_eta.Add(c_eta, newDenseFrom("./covmat/C_bias.fits"))
	c_eta.Add(c_eta, newDenseFrom("./covmat/C_cal.fits"))
	c_eta.Add(c_eta, newDenseFrom("./covmat/C_dust.fits"))
	c_eta.Add(c_eta, newDenseFrom("./covmat/C_host.fits"))
	c_eta.Add(c_eta, newDenseFrom("./covmat/C_model.fits"))
	c_eta.Add(c_eta, newDenseFrom("./covmat/C_nonia.fits"))
	c_eta.Add(c_eta, newDenseFrom("./covmat/C_pecvel.fits"))
	c_eta.Add(c_eta, newDenseFrom("./covmat/C_stat.fits"))

	// computes the A^t C_eta A part of the covariance matrix
	a_vector := mat.NewVecDense(3, []float64{1, ps[1], -ps[2]})

	c_mat_elements := make([]float64, N*N)
	i, j := 0, 0
	for k := 0; k < N*N; k++ {
		submat := c_eta.Slice(i, i+3, j, j+3)
		element := mat.Inner(a_vector, submat, a_vector)
		c_mat_elements[k] = element
		if i+3 < N*3-1 {
			i = i + 3
		} else {
			j = j + 3
			i = 0
		}
	}

	c_mat := mat.NewDense(N, N, c_mat_elements)

	// builds the covariance matrix by addind the diagonal matrices

	data, err := ioutil.ReadFile("./covmat/sigma_mu.txt")
	if err != nil {
		log.Fatal(err)
	}
	str2 := strings.Trim(string(data), " \n")
	elements := strings.Split(str2, "\n")
	elements = append(elements[:0], elements[4:]...)

	sigma_mat := mat.NewDense(N, N, nil)
	for i := 0; i < N; i++ {
		var (
			words = strings.Split(elements[i], " ")
			sz    = atof(words[4])
			slen  = atof(words[2])
			scoh  = atof(words[0])
			val   = ((5*150000)/(sz*c))*((5*150000)/(sz*c)) + slen*slen + scoh*scoh
		)
		sigma_mat.Set(i, i, val)
	}

	covmat := mat.NewDense(N, N, nil)
	covmat.Add(c_mat, sigma_mat)

	// builds the vector mu_hat - mu_lambdaCDM

	mu_vector := mat.NewVecDense(len(mu_diff), mu_diff)

	// inverses the matrix C using a Cholesky decomposition

	rows, cols := covmat.Dims()
	if rows != cols {
		log.Fatalf("cov-matrix not square")
	}
	inv := mat.NewSymDense(rows, nil)
	for i := 0; i < rows; i++ {
		for j := i; j < rows; j++ {
			inv.SetSym(i, j, covmat.At(i, j))
		}
	}
	var chol mat.Cholesky
	if ok := chol.Factorize(inv); !ok {
		log.Fatalf("cov-matrix not positive semi-definite")
	}

	err2 := chol.InverseTo(inv)
	if err2 != nil {
		log.Fatal(err2)
	}

	return mat.Inner(mu_vector, inv, mu_vector)
}

// difference computes the differences between two slices, value by value
func difference(s1, s2 []float64) ([]float64, []float64) {

	// input : s1, s2 : two slices
	// output : diff : slice of the term by term differences, abs : slice of the absolute values of the differences

	diff := make([]float64, len(s1))
	abs := make([]float64, len(s1))

	for i := range s1 {
		val := s2[i] - s1[i]
		diff[i] = val
		abs[i] = math.Abs(val)
	}

	return diff, abs
}

// atof converts a string to float64
func atof(s string) float64 {

	a := strings.Trim(s, " ")
	v, err := strconv.ParseFloat(a, 64)
	if err != nil {
		log.Panic(err)
	}
	return v
}

// integral computes the integral of the D_l function, based on the trapeze method
func integral(z, n float64) float64 {

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

// mu_th computes the theoretical mu for each supernova.
// z is the redshift.
// n is the number of steps for the integral.
func mu_th(z, n float64) float64 {
	return 5 * math.Log10(((1+z)*c/(10*H))*integral(z, n))
}

// mu_th_slice builds a slice containing the theoretical distances.
// s is a slice of the zcmb.
// n is the number of steps for the integral.
func mu_th_slice(s []float64, n float64) []float64 {
	diag_ord := make([]float64, len(s))
	for i := range s {
		diag_ord[i] = mu_th(s[i], n)
	}
	return diag_ord
}

// mu_exp_slice computes the experimental mus and add them in an array.
func mu_exp_slice(m_stell, mb, stretch, colour []float64) []float64 {

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

// mu_exp_err computes the errors on experimental mus and add them in an array
func mu_exp_err(dmb, dstretch, dcolour []float64) []float64 {

	// inputs : dmb, dstretch, dcolour = errors on mb, stretch and colour from JLA data

	slice := make([]float64, len(dmb))
	for i := range dmb {
		ai := dmb[i]
		bi := alpha * dstretch[i]
		ci := beta * dcolour[i]
		slice[i] = math.Sqrt(ai*ai + bi*bi + ci*ci)
	}

	return slice
}

// newPoints creates a plotter.XYs (cloud of points) from given 2D coordinates
func newPoints(xs, ys []float64) plotter.XYs {

	// inputs : xs, ys = slices containing the x and y values respectively

	pts := make(plotter.XYs, len(xs))
	for i := range pts {
		pts[i].X = xs[i]
		pts[i].Y = ys[i]
	}

	return pts
}

// newHubblePlot plots the given points
func newHubblePlot(points plotter.XYs) *plot.Plot {

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

// modified integral function where omgM is a variable that can be optimized
func modified_integral(z, omega float64) float64 {

	// inputs : z = redshift, omega = omgM as a free parameter

	const n = 1000
	xs := make([]float64, n+1)
	ys := make([]float64, n+1)

	for i := range xs {
		x := z * float64(i) / float64(n)
		xs[i] = x
		ys[i] = 1 / math.Sqrt((1+x*omega)*(1+x)*(1+x)-(1-omega)*(2+x)*x)
	}

	return integrate.Trapezoidal(xs, ys)
}

// builds a histogram from a given slice of values
func newHistogram(s []float64) *plot.Plot {

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

// reads a FITS fil, extracts the data and convert them into the corresponding matrix
func newDenseFrom(file string) *mat.Dense {

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
	if err != nil {
		log.Fatal(err)
	}

	return mat.NewDense(2220, 2220, raw)
}
