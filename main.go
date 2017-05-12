// Copyright 2017 The hubble-fit authors. All rights reserved.
// Use of this source code is governed by a BSD-style 
// license that can be found in the LICENSE file.

// hubble-fit is a program to fit the Hubble diagram with a Cosmological model

package main

import (
	"fmt"
	"github.com/gonum/integrate"
	"github.com/gonum/optimize"
	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/vg"
	"github.com/gonum/plot/vg/draw"
	"github.com/gonum/stat"
	"go-hep.org/x/hep/fit"
	"image/color"
	"io/ioutil"
	"log"
	"math"
	"strconv"
	"strings"
)

/*
#####
Constants
#####
*/

const (
	c       = 299792458 // speed of light
	H       = 0.070     // Hubble constant
	omgM    = 0.295     // Mass density parameter, = (8*pi*G*rho_mass)/(3*H²)
	alpha   = 0.141     // free parameter for the Hubble fit (factor for stretch)
	beta    = 3.101     // free parameter for the Hubble fit (factor for color)
	Mb      = -19.05
	delta_M = -0.070
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

	var diff_slice []float64
	var abs_diff_slice []float64

	for i,_ := range(s1) {
		val := s2[i] - s1[i]
		diff_slice = append(diff_slice, val)
		abs_diff_slice = append(abs_diff_slice, math.Abs(val))
	}

	return diff_slice, abs_diff_slice
}

/*
#####
Theoretical mu
#####
*/

func integral(z, n float64) float64 {

	// compute the integral of the D_l function, based on the trapeze method
	// integrate.Trapezoidal takes two slices (abscissa and ordinate) from which to compute the integral
	// input z : redshift, n : number of steps for the trapeze method

	var inte float64
	abs := make([]float64, 0)
	ord := make([]float64, 0)

	var i float64
	for i = 0; i <= n; i++ {
		val := z * i / n
		abs = append(abs, val)
	}

	for j := 0; j <= int(n); j++ {
		val := 1 / math.Sqrt((1+abs[j]*omgM)*(1+abs[j])*(1+abs[j])-(1-omgM)*(2+abs[j])*abs[j])
		ord = append(ord, val)
	}

	inte = integrate.Trapezoidal(abs, ord)

	return inte

}

func mu_th(z, n float64) float64 {

	// compute the theoretical mu for each SN
	// input z : redshift, n : number of steps for the integral, to be passed to integral function

	var mu_th float64
	mu_th = 5 * math.Log10(((1+z)*c/(10*H))*integral(z, n))

	return mu_th
}

func mu_th_slice(s []float64, n float64) []float64 {

	// compute the theoretical distances
	// input s : slice of the zcmb, n : number of steps for the integral, to be passed to the integral function through mu_th

	diag_ord := make([]float64, 0)

	for i := 0; i < len(s); i++ {
		diag_ord = append(diag_ord, mu_th(s[i], n))
	}

	return diag_ord
}

/*
######
Experimental mu
######
*/

func mu_exp_slice(m_stell, mb, stretch, colour []float64) []float64 {

	// compute experimental mu and add them in an array
	// input mB : B band peak magnitude, stretch and colour : SALT2 parameters

	slice := make([]float64, 0)
	for i := 0; i < len(mb); i++ {
		if m_stell[i] < 10 {
			val := mb[i] - Mb + alpha*stretch[i] - beta*colour[i]
			slice = append(slice, val)
		} else {
			val := mb[i] - Mb - delta_M + alpha*stretch[i] - beta*colour[i]
			slice = append(slice, val)
		}
	}

	return slice
}

func mu_exp_err(dmb, dstretch, dcolour []float64) []float64 {

	// compute errors on experimental mus and add them in an array

	slice := make([]float64, 0)
	for i := 0; i < len(dmb); i++ {
		val := math.Sqrt(math.Pow(dmb[i], 2) + math.Pow(alpha*dstretch[i], 2) + math.Pow(beta*dcolour[i], 2))
		slice = append(slice, val)
	}

	return slice
}

/*
######
Chi squared function
######
*/

func chi_squared(mu_th_slice, mu_exp_slice, mu_exp_err []float64) float64 {

	//

	var chi2 float64 = 0

	for i := 0; i < len(mu_th_slice); i++ {
		val := (mu_exp_slice[i] - mu_th_slice[i]) / mu_exp_err[i]
		chi2 = chi2 + math.Pow(val, 2)
	}

	return chi2
}

/*
######
Hubble diagram
######
*/

func points(x_values, y_values []float64) plotter.XYs {

	// creates the points from the datas
	// input n : number of supernovae

	pts := make(plotter.XYs, len(x_values))
	for i := range pts {
		pts[i].X = x_values[i]
		pts[i].Y = y_values[i]
	}

	return pts
}

func diag(points plotter.XYs) *plot.Plot {

	// plot the given points

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

func blabla(z, omega float64) float64 {

	var inte float64
	abs := make([]float64, 0)
	ord := make([]float64, 0)

	var i float64
	for i = 0; i <= 1000; i++ {
		val := z * i / 1000
		abs = append(abs, val)
	}

	for j := 0; j <= 1000; j++ {
		val := 1 / math.Sqrt((1+abs[j]*omega)*(1+abs[j])*(1+abs[j])-(1-omega)*(2+abs[j])*abs[j])
		ord = append(ord, val)
	}

	inte = integrate.Trapezoidal(abs, ord)

	return inte

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
	str := strings.Trim(string(datas), "# \n")		//gets rid of first and last symbols
	supernovae := strings.Split(str, "\n")
	supernovae = append(supernovae[:0], supernovae[1:]...)	//suppress the first item (labels of the columns)

	// creates slices with the needed values from JLA datas

	sn_names := make([]string, 0)
	zcmb := make([]float64, 0)
	zhel := make([]float64, 0)
	mb := make([]float64, 0)
	dmb := make([]float64, 0)
	stretch := make([]float64, 0)
	dstretch := make([]float64, 0)
	colour := make([]float64, 0)
	dcolour := make([]float64, 0)
	m_stell := make([]float64, 0)
	exp := make([]float64, 0)
	ra_jla := make([]float64, 0)
	de_jla := make([]float64, 0)

	for _, v := range supernovae {
		split_str := strings.Split(v, " ")
		for x, y := range split_str {
			switch x {
			case 0:
				sn_names = append(sn_names, y)
			case 1:
				z, _ := strconv.ParseFloat(y, 64)
				zcmb = append(zcmb, z)
			case 2:
				z, _ := strconv.ParseFloat(y, 64)
				zhel = append(zhel, z)
			case 4:
				z, _ := strconv.ParseFloat(y, 64)
				mb = append(mb, z)
			case 5:
				z, _ := strconv.ParseFloat(y, 64)
				dmb = append(dmb, z)
			case 6:
				z, _ := strconv.ParseFloat(y, 64)
				stretch = append(stretch, z)
			case 7:
				z, _ := strconv.ParseFloat(y, 64)
				dstretch = append(dstretch, z)
			case 8:
				z, _ := strconv.ParseFloat(y, 64)
				colour = append(colour, z)
			case 9:
				z, _ := strconv.ParseFloat(y, 64)
				dcolour = append(dcolour, z)
			case 10:
				z, _ := strconv.ParseFloat(y, 64)
				m_stell = append(m_stell, z)
			case 17:
				z, _ := strconv.ParseFloat(y, 64)
				exp = append(exp, z)
			case 18:
				z, _ := strconv.ParseFloat(y, 64)
				ra_jla = append(ra_jla, z)
			case 19:
				z, _ := strconv.ParseFloat(y, 64)
				de_jla = append(de_jla, z)
			}
		}
	}

	// compute and store slices needed for the plots and analysis

	var muexp_data []float64		// experimental mus from SALT2 model
	muexp_data = mu_exp_slice(m_stell, mb, stretch, colour)

	var muth_data []float64			// theoretical mus from lambdaCDM model
	muth_data = mu_th_slice(zcmb, 10000)

	var diff_data []float64			// difference between experimental and theoretical mus
	var abs_diff_data []float64		// absolute values of the differences
	diff_data, abs_diff_data = difference(muth_data, muexp_data)

	var mean_residuals float64		// mean of the differences
	mean_residuals = stat.Mean(diff_data, nil)

	var abs_mean_residuals float64		// mean of the absolute differences
	abs_mean_residuals = stat.Mean(abs_diff_data, nil)

	var muexp_error []float64		// slice containing the error on mu_exp, from error propagation in SALT2 model
	muexp_error = mu_exp_err(dmb, dstretch, dcolour)

	// compute needed parameters

	var khi_squared float64
	khi_squared = chi_squared(muth_data, muexp_data, muexp_error)

	fmt.Println("khi squared ", khi_squared)
	fmt.Println("mean of the residuals ", mean_residuals)
	fmt.Println("mean of the absolute residuals ", abs_mean_residuals)

	ps := make([]float64, 0)
	ps = append(ps, 0.295)
	res,_ := fit.Curve1D(
		fit.Func1D{
			F: func(z float64, ps []float64) float64 {
				return 5*math.Log10( ((1+z)*299792458/(10*0.070)) * blabla(z, ps[0]) )
			},
			X: zcmb,
			Y: muexp_data,
			Err: muexp_error,
			N: 1,
		},
		nil, &optimize.NelderMead{},
	)

	fmt.Println("Omega M : ", res.X[0])

	// compute, draw and store the various diagrams

	scatter_data := points(zcmb, muexp_data)
	diagram := diag(scatter_data)
	diagram.Save(1080, 720, "hubble_diagram.png")

	scatter2_data := points(zcmb, diff_data)
	diagram2 := diag(scatter2_data)
	diagram2.Save(720, 200, "hubble_diagram2.png")

	o, err := plot.New()
	if err != nil {
		log.Panic(err)
	}
	plop := make(plotter.Values, len(diff_data))
	for i,_ := range diff_data {
		plop[i] = diff_data[i]
	}
	o.Title.Text = "Histogram"
	h, err := plotter.NewHist(plop, len(plop))
	if err != nil {
		log.Panic(err)
	}
	o.Add(h)

	err = o.Save(200, 200, "histogram.png")
	if err!= nil {
		log.Panic(err)
	}

}
