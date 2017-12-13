// Copyright 2017 The hubble-fit Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/astrogo/fitsio"
	"gonum.org/v1/gonum/diff/fd"
	"gonum.org/v1/gonum/integrate"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/optimize"
)

type context struct {
	p     optimize.Problem
	jla   jla
	ceta  *mat.Dense
	sigma *mat.Dense

	pools struct {
		mu_exp  sync.Pool
		mu_th   sync.Pool
		mu_diff sync.Pool

		covbuf sync.Pool
		elmts  sync.Pool
		xs     sync.Pool
		ys     sync.Pool
	}
}

func newContext() (*context, error) {
	var (
		err error
		ctx context
	)
	ctx.jla, err = newJLA("./data/jla_lcparams.txt")
	if err != nil {
		return nil, err
	}

	n := ctx.jla.n
	nn := n * n

	pool := func(f func() interface{}) sync.Pool { return sync.Pool{New: f} }

	ctx.pools.mu_exp = pool(func() interface{} { return make([]float64, n) })
	ctx.pools.mu_th = pool(func() interface{} { return make([]float64, n) })
	ctx.pools.mu_diff = pool(func() interface{} { return make([]float64, n) })

	ctx.pools.xs = pool(func() interface{} {
		const N = 1000
		o := make([][]float64, n)
		for i := range o {
			o[i] = make([]float64, N+1)
		}
		return o
	})
	ctx.pools.ys = pool(func() interface{} {
		o := make([][]float64, n)
		const N = 1000
		for i := range o {
			o[i] = make([]float64, N+1)
		}
		return o
	})

	ctx.pools.elmts = pool(func() interface{} { return make([]float64, nn) })
	ctx.pools.covbuf = pool(func() interface{} { return make([]float64, nn) })

	ctx.ceta, err = readCeta("./data/covmat")
	if err != nil {
		return nil, err
	}

	ctx.sigma, err = readSigmaMu("./data/covmat/sigma_mu.txt", ctx.jla.n)
	if err != nil {
		return nil, err
	}

	return &ctx, nil
}

func (ctx *context) fit(ps []float64, settings *optimize.Settings, method optimize.Method) (*optimize.Result, error) {
	res, err := FitChi2(ctx, ctx.chi2, ps, settings, method)
	if err != nil {
		return res, err
	}

	if res.Hessian == nil {
		var fds fd.Settings
		fds.Concurrent = settings != nil && settings.Concurrent > 0
		res.Hessian = mat.NewSymDense(len(ps), nil)
		fd.Hessian(res.Hessian, ctx.chi2, res.X, &fds)
	}

	return res, err
}

// chi2 is the test function for a chi2 computation that takes into account the covariance matrix
func (ctx *context) chi2(ps []float64) float64 {
	N := ctx.jla.n
	mu_exp := ctx.pools.mu_exp.Get().([]float64)
	mu_th := ctx.pools.mu_th.Get().([]float64)
	mu_diff := ctx.pools.mu_diff.Get().([]float64)
	defer ctx.pools.mu_exp.Put(mu_exp)
	defer ctx.pools.mu_th.Put(mu_th)
	defer ctx.pools.mu_diff.Put(mu_diff)

	xs := ctx.pools.xs.Get().([][]float64)
	ys := ctx.pools.ys.Get().([][]float64)
	defer ctx.pools.xs.Put(xs)
	defer ctx.pools.ys.Put(ys)

	for i := 0; i < N; i++ {
		if ctx.jla.m_stell[i] < 10 {
			mu_exp[i] = ctx.jla.mb[i] - ps[3] + ps[1]*ctx.jla.stretch[i] - ps[2]*ctx.jla.colour[i]
		} else {
			mu_exp[i] = ctx.jla.mb[i] - ps[3] - ps[4] + ps[1]*ctx.jla.stretch[i] - ps[2]*ctx.jla.colour[i]
		}

		const n = 1000
		xsi := xs[i]
		ysi := ys[i]

		for j := range xsi {
			xsi[j] = ctx.jla.zcmb[i] * float64(j) / float64(n)
			ysi[j] = 1 / math.Sqrt((1+xsi[j]*ps[0])*(1+xsi[j])*(1+xsi[j])-(1-ps[0])*(2+xsi[j])*xsi[j])
		}

		mu_th[i] = 5 * math.Log10(((1+ctx.jla.zcmb[i])*c/(10*H))*integrate.Trapezoidal(xsi, ysi))

		mu_diff[i] = mu_exp[i] - mu_th[i]
	}

	c_eta := ctx.ceta

	// computes the A^t C_eta A part of the covariance matrix
	a_vector := mat.NewVecDense(3, []float64{1, ps[1], -ps[2]})

	c_mat_elements := ctx.pools.elmts.Get().([]float64)
	defer ctx.pools.elmts.Put(c_mat_elements)
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
	covbuf := ctx.pools.covbuf.Get().([]float64)
	defer ctx.pools.covbuf.Put(covbuf)
	covmat := mat.NewDense(N, N, covbuf)
	covmat.Add(c_mat, ctx.sigma)

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

	err := chol.InverseTo(inv)
	if err != nil {
		log.Fatal(err)
	}

	return mat.Inner(mu_vector, inv, mu_vector)
}

type jla struct {
	n        int       // number of supernovae
	names    []string  // names of the supernovae
	zcmb     []float64 // redshift in the cmb frame
	zhel     []float64 // redshift in the heliocentric frame
	mb       []float64 // b band peak magnitude
	dmb      []float64 // b band peak magnitude error
	stretch  []float64 // stretch factor
	dstretch []float64 // stretch factor error
	colour   []float64 // colour factor
	dcolour  []float64 // colour factor error
	m_stell  []float64 // log10 of the host galaxy stellar mass
	exp      []float64 // supernova sample identifier : 1 = SNLS, 2 = SDSS, 3 = lowz, 4 = Riess HST
	ra_jla   []float64 // right ascension (in degrees)
	de_jla   []float64 // declination (in degrees)
}

func newJLA(name string) (jla, error) {
	var params jla
	data, err := ioutil.ReadFile(name)
	if err != nil {
		return params, err
	}

	str := strings.Trim(string(data), "\n")                // gets rid, if needed, of a new line symbol at the end of the file
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
			return params, fmt.Errorf("hubble: invalid number of jla_lcparams fields. got=%d, want=%d", len(split_str), 21)
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

	params.n = N
	params.names = sn_names
	params.zcmb = zcmb
	params.zhel = zhel
	params.mb = mb
	params.dmb = dmb
	params.stretch = stretch
	params.dstretch = dstretch
	params.colour = colour
	params.dcolour = dcolour
	params.m_stell = m_stell
	params.exp = exp
	params.ra_jla = ra_jla
	params.de_jla = de_jla
	return params, nil
}

func readCeta(dir string) (*mat.Dense, error) {
	var err error
	ceta := mat.NewDense(2220, 2220, nil)
	for _, name := range []string{
		"C_bias.fits",
		"C_cal.fits",
		"C_dust.fits",
		"C_host.fits",
		"C_model.fits",
		"C_nonia.fits",
		"C_pecvel.fits",
		"C_stat.fits",
	} {
		name = filepath.Join(dir, name)
		m, err := newDenseFrom2(name)
		if err != nil {
			return nil, err
		}
		ceta.Add(ceta, m)
	}

	return ceta, err
}

func readSigmaMu(name string, N int) (*mat.Dense, error) {
	data, err := ioutil.ReadFile(name)
	if err != nil {
		return nil, err
	}
	str := strings.Trim(string(data), " \n")
	elements := strings.Split(str, "\n")
	elements = append(elements[:0], elements[4:]...)

	sigma := mat.NewDense(N, N, nil)
	for i := 0; i < N; i++ {
		var (
			words = strings.Split(elements[i], " ")
			sz    = atof(words[4])
			slen  = atof(words[2])
			scoh  = atof(words[0])
			val   = ((5*150000)/(sz*c))*((5*150000)/(sz*c)) + slen*slen + scoh*scoh
		)
		sigma.Set(i, i, val)
	}
	return sigma, nil
}

// reads a FITS fil, extracts the data and convert them into the corresponding matrix
func newDenseFrom2(file string) (*mat.Dense, error) {

	r, err := os.Open(file)
	if err != nil {
		return nil, err
	}
	defer r.Close()

	fits, err := fitsio.Open(r)
	if err != nil {
		return nil, err
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
		return nil, err
	}

	return mat.NewDense(2220, 2220, raw), nil
}
