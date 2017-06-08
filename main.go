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

	"github.com/gonum/integrate"
	"github.com/gonum/optimize"
	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/vg"
	"github.com/gonum/plot/vg/draw"
	"github.com/gonum/stat"
	"go-hep.org/x/hep/fit"
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

	v, err := strconv.ParseFloat(s, 64)
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


//	for _,v := range supernovae {
//		split_str := strings.Split(v, " ")


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
	muexp_error := mu_exp_err(dmb, dstretch, dcolour)		// errors on mu_exp, from error propagation in SALT2 model

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

	// fitting function without covariance matrix

	empty := make([]float64, len(zcmb))
	params := []float64{0.295, 0.141, 3.101, -19.05, -0.070}
	res, err := fit.Curve1D(
		fit.Func1D{
			F: func(z float64, ps []float64) float64 {
				var i int = -1
				for p, v := range zcmb {
					if v == z {
						i = p
					}
				}
				if i < 0 {
					panic("error in retrieving the index of zcmb in the fitting function")
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
	fmt.Println("Delta_M : ", res.X[4])
}
