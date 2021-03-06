#include <Rcpp.h>
#include <RcppGSL.h>
//
#include <cstdlib>
#include <map>
#include <algorithm>
//
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
//
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <boost/algorithm/cxx11/iota.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

namespace { // anonymous

  using namespace Rcpp;

  // declare types
  typedef boost::tuple<int,int,double,double> P_ij_key;
  typedef std::map<P_ij_key,gsl_vector*> dP_ij_t;
  //forward declaration(s)
  int funcReclassified (double t, const double y[], double f[],
			void *model);
  int funcReclassifiedB (double t, const double y[], double f[],
			   void *model);
  int funcReclassifiedExtended (double t, const double y[], double f[], void *model);
  // Base class for splines
  class splineBasis {
  public:
    virtual double calc(double t, gsl_vector * beta, bool centering = true) = 0;
    virtual void calcBasis(double x, gsl_vector * basis, bool centering = true) = 0;
    virtual ~splineBasis() {}
  };
  // enum's cause problems with RcppGSL::vector, so instead we use int's 
  int Never=0, Current=1, Former=2, Reclassified=3, Death=4; // Issue: Death=3 or Death=4?
  int CurrentStatus=0,Recall=1,FormerWithCessation=2;
  // utility functions
  void Rprint(gsl_matrix * mat) {
    size_t i, j;
    for (i=0; i < mat->size1; ++i) {
      for (j=0; j < mat->size2; ++j) 
	Rprintf("%f ", gsl_matrix_get(mat, i, j));
      Rprintf("\n");
    }
  }
  void Rprint(gsl_vector * vec) {
    size_t i;
    for (i=0; i < vec->size; ++i) {
      Rprintf("%f ", gsl_vector_get(vec, i));
    }
  }
  // template<class T> T bounds(T x, T lo, T hi) { return x<lo ? lo : (x>hi ? hi : x); } 
  inline
  double bounds(double x, double lo, double hi) { return x<lo ? lo : (x>hi ? hi : x); } 
  gsl_vector * as_vector(SEXP sexp) {
    Rcpp::NumericVector nv = as<Rcpp::NumericVector>(sexp);
    gsl_vector *out;
    out = gsl_vector_alloc(nv.size());
    for (int i=0; i<nv.size(); ++i)
      gsl_vector_set(out,i,nv[i]);
    return out;
  }
  gsl_matrix * as_matrix(SEXP sexp) {
    Rcpp::NumericMatrix nm = as<Rcpp::NumericMatrix>(sexp);
    gsl_matrix *out;
    out = gsl_matrix_alloc(nm.nrow(), nm.ncol());
    for (int i=0; i<nm.nrow(); ++i)
      for (int j=0; j<nm.ncol(); ++j)
	gsl_matrix_set(out,i,j,nm(i,j));
    return out;
  }
  gsl_vector_int * as_vector_int(SEXP sexp) {
    Rcpp::IntegerVector iv = as<Rcpp::IntegerVector>(sexp);
    gsl_vector_int *out;
    out = gsl_vector_int_alloc(iv.size());
    for (int i=0; i<iv.size(); ++i)
      gsl_vector_int_set(out,i,iv[i]);
    return out;
  }
  SEXP wrap_vector(gsl_vector * vec) {
    Rcpp::NumericVector nv(vec->size);
    for (size_t i=0; i<vec->size; ++i)
      nv[i] = gsl_vector_get(vec,i);
    return wrap(nv);
  }
  
  class ns : public splineBasis {
  public:
    ns(gsl_vector *knots, bool intercept = false, bool centred = false, double centre = 0.0, size_t k = 4) : 
      intercept(intercept), centred(centred) {
      size_t i, j;
      gsl_bspline_deriv_workspace *dw; 
      gsl_matrix *dB;
      gsl_matrix *const_matrix;
      gsl_vector *tau;
      gsl_matrix *Q_full;
      gsl_matrix *R;
      offset = intercept ? 0 : 1;
      nbreak = knots->size;
      ncoeffs_base = nbreak + k - 2;
      ncoeffs_int = intercept ? ncoeffs_base : ncoeffs_base - 1;
      ncoeffs = ncoeffs_int - 2;
      min_knot = gsl_vector_get(knots,0);
      max_knot = gsl_vector_get(knots, knots->size - 1);
      /* allocate workspaces */
      bw = gsl_bspline_alloc(k, nbreak);
      dw = gsl_bspline_deriv_alloc(k);
      dB = gsl_matrix_alloc(ncoeffs_base,3);
      const_matrix = gsl_matrix_alloc(ncoeffs_int,2);
      ttl = gsl_matrix_alloc(2,ncoeffs_base);
      ttu = gsl_matrix_alloc(2,ncoeffs_base);
      tau = gsl_vector_alloc(2);
      Q_full = gsl_matrix_alloc(ncoeffs_int,ncoeffs_int);
      Q_sub = gsl_matrix_alloc(ncoeffs_int,ncoeffs);
      R = gsl_matrix_alloc(ncoeffs_int,2);
      B = gsl_vector_alloc(ncoeffs_base);
      X = gsl_vector_alloc(ncoeffs_int);
      xv = gsl_vector_alloc(2);
      nsC = gsl_vector_alloc(ncoeffs);
      _basis = gsl_vector_alloc(ncoeffs);
      // set the knots (which adds the repeated knots at the boundaries)
      gsl_bspline_knots(knots, bw);
      gsl_bspline_deriv_eval (min_knot, 2, dB, bw, dw);
      for (i = 0; i < ncoeffs_int; ++i) 
	{
	  gsl_matrix_set(const_matrix, i, 0, gsl_matrix_get(dB, i+offset, 2));
	}
      for (i = 0; i < ncoeffs_base; ++i) 
	{
	  gsl_matrix_set(ttl, 0, i, gsl_matrix_get(dB, i, 0));
	  gsl_matrix_set(ttl, 1, i, gsl_matrix_get(dB, i, 1));
	}
      gsl_bspline_deriv_eval (max_knot, 2, dB, bw, dw);
      for (i = 0; i < ncoeffs_int; ++i) 
	{
	  gsl_matrix_set(const_matrix, i, 1, gsl_matrix_get(dB, i+offset, 2));
	}
      for (i = 0; i < ncoeffs_base; ++i) 
	{
	  gsl_matrix_set(ttu, 0, i, gsl_matrix_get(dB, i, 0));
	  gsl_matrix_set(ttu, 1, i, gsl_matrix_get(dB, i, 1));
	}
      gsl_linalg_QR_decomp (const_matrix, tau);
      gsl_linalg_QR_unpack (const_matrix, tau, Q_full, R);
      // drop the first two columns of Q_full
      for (i = 0; i<Q_full->size1; ++i) {
	for (j = 2; j<Q_full->size2; ++j) {
	  gsl_matrix_set(Q_sub, i, j - 2, gsl_matrix_get(Q_full, i, j));
	}
      }
      if (centred) {
	calcBasis(centre, nsC, false);
      }
      gsl_bspline_deriv_free(dw);
      gsl_matrix_free(dB);
      gsl_matrix_free(const_matrix);
      gsl_matrix_free(R);
      gsl_matrix_free(Q_full);
      gsl_vector_free(tau);
    }
    virtual ~ns() {
      gsl_bspline_free(bw);
      gsl_matrix_free(Q_sub);
      gsl_matrix_free(ttl);
      gsl_matrix_free(ttu);
      gsl_vector_free(X);
      gsl_vector_free(B);
      gsl_vector_free(xv);
      gsl_vector_free(nsC);
      gsl_vector_free(_basis);
    }
    void calcBasis(double x, gsl_vector * basis, bool centering = true) {
      size_t j;
      if (x < min_knot) {
	gsl_vector_set(xv, 0, 1.0);
	gsl_vector_set(xv, 1, x - min_knot);
	gsl_blas_dgemv(CblasTrans, 1.0, ttl, xv, 0.0, B);
      } else if (x > max_knot) {
	gsl_vector_set(xv, 0, 1.0);
	gsl_vector_set(xv, 1, x - max_knot);
	gsl_blas_dgemv(CblasTrans, 1.0, ttu, xv, 0.0, B);
      } else {
	/* compute B_j(xi) for all j */
	gsl_bspline_eval(x, B, bw);
      }
      /* fill row i of X */
      for (j = 0; j < ncoeffs_int; ++j)
	{
	  gsl_vector_set(X, j, gsl_vector_get(B, j+offset));
	}
      // calculate basis=x^T*Q_sub
      gsl_blas_dgemv (CblasTrans,
		      1.0, Q_sub, X,
		      0.0, basis);
      if (centred && centering)
	gsl_blas_daxpy (-1.0, nsC, basis);
    }
    double calc(double x, gsl_vector * beta, bool centering = true) {
      double value;
      calcBasis(x, _basis, centering);
      int status = gsl_blas_ddot(beta, _basis, &value);
      if (status != GSL_SUCCESS)
	REprintf ("error, return value=%d\n", status);
      return value;
    }
    bool intercept, centred;
    size_t nbreak, ncoeffs_base, ncoeffs_int, ncoeffs;
    int offset;
    gsl_bspline_workspace *bw;
    gsl_matrix *Q_sub, *ttl, *ttu;
    gsl_vector *B, *X, *xv, *nsC, *_basis;
    double min_knot, max_knot;
  };

  class ps : public splineBasis {
  public:
    ps(double boundaryL, double boundaryU, bool intercept = false, bool centred = false, double centre = 0.0, size_t nterm = 10, size_t k = 4) : 
      intercept(intercept), centred(centred), nterm(nterm), boundaryL(boundaryL), boundaryU(boundaryU) {
      size_t i, degree;
      int m;
      gsl_bspline_deriv_workspace *dw; 
      gsl_matrix *dB;
      double dx;
      degree = k - 1;
      nbreak = nterm + 2*degree + 1;
      knots = gsl_vector_alloc(nbreak);
      dx = (boundaryU - boundaryL) / double(nterm);
      for (i = 0, m=-int(degree); m<=(int(nterm) - 1); ++i, ++m) 
	gsl_vector_set(knots,i,boundaryL + dx*m);
      for (i = 0, m=nterm+degree; i<=degree; ++i, ++m) 
	gsl_vector_set(knots,m,boundaryU + dx*int(i));
      offset = intercept ? 0 : 1;
      ncoeffs_base = nbreak + k - 2;
      ncoeffs_int = intercept ? ncoeffs_base : ncoeffs_base - 1;
      ncoeffs = ncoeffs_int - 6;
      // Rprintf("ncoeffs=%i\n",ncoeffs);
      /* allocate workspaces */
      bw = gsl_bspline_alloc(k, nbreak);
      dw = gsl_bspline_deriv_alloc(k);
      dB = gsl_matrix_alloc(ncoeffs_base,2);
      ttl = gsl_matrix_alloc(2,ncoeffs_base);
      ttu = gsl_matrix_alloc(2,ncoeffs_base);
      B = gsl_vector_alloc(ncoeffs_base);
      xv = gsl_vector_alloc(2);
      nsC = gsl_vector_alloc(ncoeffs);
      _basis = gsl_vector_alloc(ncoeffs);
      // set the knots (which adds the repeated knots at the boundaries)
      gsl_bspline_knots(knots, bw);
      gsl_bspline_deriv_eval (boundaryL, 1, dB, bw, dw);
      for (i = 0; i < (size_t) ncoeffs_base; ++i) 
	{
	  gsl_matrix_set(ttl, 0, i, gsl_matrix_get(dB, i, 0));
	  gsl_matrix_set(ttl, 1, i, gsl_matrix_get(dB, i, 1));
	}
      gsl_bspline_deriv_eval (boundaryU, 1, dB, bw, dw);
      for (i = 0; i < (size_t) ncoeffs_base; ++i) 
	{
	  gsl_matrix_set(ttu, 0, i, gsl_matrix_get(dB, i, 0));
	  gsl_matrix_set(ttu, 1, i, gsl_matrix_get(dB, i, 1));
	}
      if (centred) {
	calcBasis(centre, nsC, false);
      }
      gsl_bspline_deriv_free(dw);
      gsl_matrix_free(dB);
    }
    virtual ~ps() {
      gsl_bspline_free(bw);
      gsl_matrix_free(ttl);
      gsl_matrix_free(ttu);
      gsl_vector_free(B);
      gsl_vector_free(xv);
      gsl_vector_free(knots);
      gsl_vector_free(nsC);
      gsl_vector_free(_basis);
    }
    void calcBasis(double x, gsl_vector * basis, bool centering = true) {
      size_t j;
      if (x < boundaryL) {
	gsl_vector_set(xv, 0, 1.0);
	gsl_vector_set(xv, 1, x - boundaryL);
	gsl_blas_dgemv(CblasTrans, 1.0, ttl, xv, 0.0, B);
      } else if (x > boundaryU) {
	gsl_vector_set(xv, 0, 1.0);
	gsl_vector_set(xv, 1, x - boundaryU);
	gsl_blas_dgemv(CblasTrans, 1.0, ttu, xv, 0.0, B);
      } else {
	/* compute B_j(xi) for all j */
	gsl_bspline_eval(x, B, bw);
      }
      /* fill basis */
      for (j = 0; j < (size_t) ncoeffs; ++j)
	{
	  gsl_vector_set(basis, j, gsl_vector_get(B, j+3+offset));
	}
      if (centred && centering)
	gsl_blas_daxpy (-1.0, nsC, basis);
    }
    double calc(double x, gsl_vector * beta, bool centering = true) {
      double value;
      calcBasis(x, _basis, centering);
      int status = gsl_blas_ddot(beta, _basis, &value);
      if (status != GSL_SUCCESS)
	REprintf ("error, return value=%d\n", status);
      return value;
    }
    bool intercept, centred;
    int nbreak, ncoeffs_base, ncoeffs_int, ncoeffs, nterm;
    int offset;
    double boundaryL, boundaryU;
    gsl_bspline_workspace *bw;
    gsl_matrix *ttl, *ttu;
    gsl_vector *knots, *B, *xv, *nsC, *_basis;
  };

  /** @brief Base class for purged Markov chains to model smoking.
      It is specialised for smoking, including the data fields, the smoothers 
   **/
  class Purged {
  public:
    bool debug;
    // Params fields
    double int01, int12, beta20, maxage;
    size_t N, nstate, ncoef, nbeta01, nbeta12;
    gsl_spline *spline0, *spline1, *spline2;
    gsl_interp_accel *acc0, *acc1, *acc2;
    splineBasis *s01;
    splineBasis *s12;
    std::map<P_ij_key,double> P_ijs;
    std::map<P_ij_key,gsl_vector *> dP_ijs; // map of pointers: care required
    gsl_odeiv2_system sys;
    gsl_odeiv2_driver *d;
    // data fields
    gsl_vector_int *recall, *finalState;
    gsl_vector *time1, *time2, *time3, *freq, *beta01, *beta12, *ages0, *mu0, *mu1, *mu2;
    Purged(SEXP sexp) {
      List args = as<List>(sexp);
      debug = args("debug");
      finalState = as_vector_int(args("finalState"));
      recall = as_vector_int(args("recall"));
      time1 = as_vector(args("time1"));
      time2 = as_vector(args("time2"));
      time3 = as_vector(args("time3"));
      freq = as_vector(args("freq"));
      beta01 = as_vector(args("beta01"));
      beta12 = as_vector(args("beta12"));
      ages0 = as_vector(args("ages0"));
      mu0 = as_vector(args("mu0"));
      mu1 = as_vector(args("mu1"));
      mu2 = as_vector(args("mu2"));
      int01 = args("int01"); 
      int12 = args("int12"); 
      beta20 = args("beta20"); 
      N = time1->size;
      nstate = 4;
      nbeta01 = beta01->size;
      nbeta12 = beta12->size;
      ncoef = 3+nbeta01+nbeta12;
      acc0 = gsl_interp_accel_alloc ();
      acc1 = gsl_interp_accel_alloc ();
      acc2 = gsl_interp_accel_alloc ();
      spline0 = gsl_spline_alloc (gsl_interp_cspline, ages0->size);
      spline1 = gsl_spline_alloc (gsl_interp_cspline, ages0->size);
      spline2 = gsl_spline_alloc (gsl_interp_cspline, ages0->size);
      gsl_spline_init (spline0, ages0->data, mu0->data, ages0->size);
      gsl_spline_init (spline1, ages0->data, mu1->data, ages0->size);
      gsl_spline_init (spline2, ages0->data, mu2->data, ages0->size);
      maxage = gsl_vector_get(ages0,ages0->size-1)+0.5;
      if (debug)
	Rprintf("maxage=%f\n",maxage);
      sys.jacobian = 0; 
      sys.params = (void *) this;
      sys.function = funcReclassified;
      sys.dimension = nstate;
      d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
					 1e-10, 1e-10, 0.0);
    }
    virtual ~Purged() {
      gsl_vector_int_free(recall);
      gsl_vector_int_free(finalState);
      gsl_vector_free(time1);
      gsl_vector_free(time2);
      gsl_vector_free(time3);
      gsl_vector_free(freq);
      gsl_vector_free(beta01);
      gsl_vector_free(beta12);
      gsl_vector_free(ages0);
      gsl_vector_free(mu0);
      gsl_vector_free(mu1);
      gsl_vector_free(mu2);
      gsl_spline_free (spline0);
      gsl_spline_free (spline1);
      gsl_spline_free (spline2);
      gsl_interp_accel_free (acc0);
      gsl_interp_accel_free (acc1);
      gsl_interp_accel_free (acc2);
      gsl_odeiv2_driver_free (d);
      P_ijs.clear();
      dP_ijs_clear();
    }
    void dP_ijs_clear() {
      for (dP_ij_t::iterator it = dP_ijs.begin(); it != dP_ijs.end(); it++)
	gsl_vector_free(it->second);
      dP_ijs.clear();
    }
    void default_sys() {
      sys.function = funcReclassified;
      sys.dimension = nstate;
      gsl_odeiv2_driver_free (d);
      d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
					 1e-10, 1e-10, 0.0);
    }
    void extended_sys() {
      sys.function = funcReclassifiedExtended;
      sys.dimension = nstate*(ncoef+1);
      gsl_odeiv2_driver_free (d);
      d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
					 1e-10, 1e-10, 0.0);
    }
    double P_ij(int i, int j, double s, double t, double eps = 1e-10) {
      // do we have the values cached?
      P_ij_key key = P_ij_key(i,j,s,t);
      if (P_ijs.count(key)>0) return P_ijs[key];
      // otherwise...
      default_sys();
      double y[sys.dimension];
      for (size_t k=0; k<nstate; ++k)
	y[k] = 0.0;
      y[i] = 1.0;
      double tm = s;
      int status = gsl_odeiv2_driver_apply (d, &tm, t, y);
      if (status != GSL_SUCCESS)
	REprintf ("error, return value=%d\n", status);
      // cache the results
      for (size_t k=0; k<nstate; ++k)
	P_ijs[P_ij_key(i,k,s,t)] = y[k]<eps && -eps<y[k] ? eps : y[k];
      return P_ijs[key];
    }
    double P_iK(int i, double s, double t) {
      double total = 0.0;
      for (size_t k=0; k<nstate; ++k)
	total += P_ij(i,k,s,t);
      return total;
    }
    gsl_vector *
    dP_ij(int i, int j, double s, double t, double eps = 1e-10) {
      // do we have the values cached?
      P_ij_key key = P_ij_key(i,j,s,t);
      if (dP_ijs.find(key) != dP_ijs.end()) return dP_ijs.find(key)->second;
      // otherwise...
      extended_sys();
      double y[sys.dimension];
      for (size_t k=0; k<sys.dimension; ++k)
	y[k] = 0.0;
      y[i] = 1.0;
      double tm = s;
      int status = gsl_odeiv2_driver_apply (d, &tm, t, y);
      if (status != GSL_SUCCESS)
	REprintf ("error, return value=%d\n", status);
      // cache the results
      for (size_t k=0; k<nstate; ++k) { // state
	gsl_vector * out; 
	out = gsl_vector_alloc(ncoef);
	P_ij_key key_ik = P_ij_key(i,k,s,t);
	P_ijs[key_ik] = y[k]<eps && -eps<y[k] ? eps : y[k];
	for (size_t m=0; m<ncoef; ++m) 
	  gsl_vector_set(out, m, y[nstate*(m+1)+k]); // (y, P1', P2', ..., Pncoef')
	dP_ijs.insert(dP_ij_t::value_type(key_ik, out));
      }
      return dP_ijs.find(key)->second;
    }
    gsl_vector * dP_iK(int i, double s, double t) {
      gsl_vector * total;
      total = gsl_vector_calloc(ncoef); // freed in negll_gradient_component()
      for (size_t k=0; k<nstate; ++k) {
	gsl_vector_add(total,dP_ij(i,k,s,t));
      }
      return total;
    }
    // Excess mortality model
    double
    mu_ij (int i, int j, double t, double eps = 1e-16)
    {
      if (i==Never && j==Current)  return exp(int01 + s01->calc(t, beta01));
      if (i==Current && j==Former) return exp(int12 + s12->calc(t, beta12));
      if (i==Former && j==Reclassified)  return exp(beta20);
      if (j==Death) {
	if (i==Never) return 0.0;
	if (i==Current) 
	  return std::max(eps,
			  gsl_spline_eval(spline1, bounds(t,0.5,maxage), acc1) -
			  gsl_spline_eval(spline0, bounds(t,0.5,maxage), acc0));
	if (i==Former) 
	  return std::max(eps,
			  gsl_spline_eval(spline2, bounds(t,0.5,maxage), acc2) - 
			  gsl_spline_eval(spline0, bounds(t,0.5,maxage), acc0));
	if (i==Reclassified) return 0.0;
      }
      REprintf("Unexpected combinations of states in mu_ij: %i to %j.\n",i,j);
      return 0.0;
    }
    // Model B: observed mortality
    double
    mu_ijB(int i, int j, double t)
    {
      if (i==Never && j==Current)  return exp(int01 + s01->calc(t, beta01));
      if (i==Current && j==Former) return exp(int12 + s12->calc(t, beta12));
      if (i==Former && j==Reclassified)  return exp(beta20);
      if (j==Death) {
	if (i==Never) return gsl_spline_eval(spline0, bounds(t,0.5,maxage), acc0);
	if (i==Current) return gsl_spline_eval(spline1, bounds(t,0.5,maxage), acc1);
	if (i==Former) return gsl_spline_eval(spline2, bounds(t,0.5,maxage), acc2);
	if (i==Reclassified) return gsl_spline_eval(spline0, bounds(t,0.5,maxage), acc0);
      }
      REprintf("Unexpected combinations of states in mu_ijB: %i to %j.\n",i,j);
      return 0.0;
    }
    /** @brief Rate of change wrt theta_u (cf. partial) of mu_ij = (partial mu_ij / partial theta_u) / mu_ij 
     **/
    gsl_vector * rmu_ij(int i, int j, double t)
    {
      gsl_vector * x;
      x = gsl_vector_calloc(ncoef); // freed in negll_gradient_component()
      int k = 0;
      // set theta and x 
      // intercept
      if (i==Never && j==Current) {
	k=0;
	gsl_vector_set(x, k++, 1);
	// spline parameters
	gsl_vector * _beta01; 
	_beta01 = gsl_vector_alloc(nbeta01);
	s01->calcBasis(t, _beta01);
	for (size_t n = 0; n<nbeta01; ++n) {
	  gsl_vector_set(x, k++, gsl_vector_get(_beta01,n));
	}
	gsl_vector_free(_beta01);
      }
      else if (i==Current && j==Former) {
	k=nbeta01+1-1;
	// intercept
	gsl_vector_set(x, k++, 1);
	// spline parameters
	gsl_vector * _beta12;
	_beta12 = gsl_vector_alloc(nbeta12);
	s12->calcBasis(t, _beta12);
	for (size_t n = 0; n<nbeta12; ++n) {
	  gsl_vector_set(x, k++, gsl_vector_get(_beta12,n));
	}
	gsl_vector_free(_beta12);
      }
      else if (i==Former && j==Reclassified) {
	k=nbeta01+nbeta12+2-1;
	gsl_vector_set(x, k++, 1);
      }
      else REprintf("Unexpected combinations of states in rmu_ij: %i to %j.\n",i,j);
      return x;
    }
    double negll_component(int state, double s, double t, double u, double freq = 1.0, int recall = 1) {
      double ll = 0.0;
      if (recall == Recall) { // recall of smoking history available
	if (state == Never) // Never->Never _and_ Never->Reclassified (equiv to current status)
	  ll = log(P_ij(Never,Never,0.0,u)+
		   P_ij(Never,Reclassified,0.0,u)); // ignore s and t
	if (state == Current) {
	  if (debug) 
	    Rprintf("Components: (%f,%f,%f,%f,%f)\n",P_ij(Never,Never,0.0,s),mu_ij(Never,Current,s),P_ij(Current,Current,s,u));
	  ll = log(P_ij(Never,Never,0.0,s))+
	    log(mu_ij(Never,Current,s))+
	    log(P_ij(Current,Current,s,u)); // ignore t
	}
	if (state == Former) {
	  if (debug) 
	    Rprintf("Components: (%f,%f,%f,%f,%f)\n",P_ij(Never,Never,0.0,s),mu_ij(Never,Current,s),P_ij(Current,Current,s,t),mu_ij(Current,Former,t),P_ij(Former,Former,t,u));
	  ll = log(P_ij(Never,Never,0.0,s))+
	    log(mu_ij(Never,Current,s))+
	    log(P_ij(Current,Current,s,t))+
	    log(mu_ij(Current,Former,t))+
	    log(P_ij(Former,Former,t,u));
	}
      }
      if (recall == CurrentStatus) {// current status only
	if (state == Never) // Never->Never _and_ Never->Reclassified (equiv to current status)
	  ll = log(P_ij(Never,Never,0.0,u)+
		   P_ij(Never,Reclassified,0.0,u)); // ignore s and t
	else 
	  ll = log(P_ij(Never,state,0.0,u));
      }
      if (recall == FormerWithCessation && state == Former) {// recall of age quit (not initiation) for former smokers
	ll = log(P_ij(Never,Current,0.0,t)) + 
	  log(mu_ij(Current,Former,t))+
	  log(P_ij(Former,Former,t,u)); // ignores s
      }
      if (ll == 0.0) REprintf("ll==0.0? (state=%i, recall=%i)\n",state,recall);
      if (debug) Rprintf("(state=%i, recall=%i, s=%f,t=%f,u=%f) ll=%f\n",state,recall,s,t,u,ll);
      ll -= log(P_iK(Never,0.0,u));
      return -ll*freq;
    }
    gsl_vector * negll_gradient_component(int state, double s, double t, double u, double freq = 1.0, int recall = 1) {
      gsl_vector * gradient;
      gradient = gsl_vector_calloc(ncoef);
      if (recall == Recall) { // recall of smoking history available
	if (state == Never) // Never->Never _and_ Never->Reclassified (equiv to current status)
	  for (size_t k=0; k<ncoef; ++k)
	    gsl_vector_set(gradient, k,
			   (gsl_vector_get(dP_ij(Never,Never,0.0,u), k)+
			    gsl_vector_get(dP_ij(Never,Reclassified,0.0,u), k)) / 
			   (P_ij(Never,Never,0.0,u)+P_ij(Never,Reclassified,0.0,u))); // ignore s and t
	if (state == Current) {
	  gsl_vector * rmu01 = rmu_ij(Never,Current,s);
	  for (size_t k=0; k<ncoef; ++k) 
	    gsl_vector_set(gradient, k, 
			   gsl_vector_get(dP_ij(Never,Never,0.0,s),k) / P_ij(Never,Never,0.0,s)+
			   gsl_vector_get(rmu01,k)+
			   gsl_vector_get(dP_ij(Current,Current,s,u),k)/P_ij(Current,Current,s,u)); // ignore t
	  gsl_vector_free(rmu01);
	}
	if (state == Former) {
	  gsl_vector *rmu01 = rmu_ij(Never,Current,s), 
	    *rmu12 = rmu_ij(Current,Former,t);
	  for (size_t k=0; k<ncoef; ++k)
	    gsl_vector_set(gradient, k,
			   gsl_vector_get(dP_ij(Never,Never,0.0,s),k)/P_ij(Never,Never,0.0,s)+
			   gsl_vector_get(rmu01,k) +
			   gsl_vector_get(dP_ij(Current,Current,s,t),k)/P_ij(Current,Current,s,t)+
			   gsl_vector_get(rmu12,k) +
			   gsl_vector_get(dP_ij(Former,Former,t,u),k)/P_ij(Former,Former,t,u));
	  gsl_vector_free(rmu01); gsl_vector_free(rmu12);
	}
      }
      if (recall == CurrentStatus) {// current status only
	if (state == Never) // Never->Never _and_ Never->Reclassified (equiv to current status)
	  for (size_t k=0; k<ncoef; ++k)
	    gsl_vector_set(gradient, k,
			   (gsl_vector_get(dP_ij(Never,Never,0.0,u),k)+
			    gsl_vector_get(dP_ij(Never,Reclassified,0.0,u),k)) / 
			   (P_ij(Never,Never,0.0,u)+P_ij(Never,Reclassified,0.0,u))); // ignore s and t
	else 
	  for (size_t k=0; k<ncoef; ++k)
	    gsl_vector_set(gradient, k, 
			   gsl_vector_get(dP_ij(Never,state,0.0,u),k)/P_ij(Never,state,0.0,u));
      }
      if (recall == FormerWithCessation && state == Former) {// recall of age quit (not initiation) for former smokers
	gsl_vector * rmu12 = rmu_ij(Current,Former,t);
	for (size_t k=0; k<ncoef; ++k)
	  gsl_vector_set(gradient, k, 
			 gsl_vector_get(dP_ij(Never,Current,0.0,t),k)/P_ij(Never,Current,0.0,t) + 
			 gsl_vector_get(rmu12,k)+
			 gsl_vector_get(dP_ij(Former,Former,t,u),k)/P_ij(Former,Former,t,u)); // ignores s
	gsl_vector_free(rmu12);
      }
      // if (ll == 0.0) REprintf("ll==0.0? (state=%i, recall=%i)\n",state,recall);
      gsl_vector * dP = dP_iK(Never,0.0,u);
      for (size_t k=0; k<ncoef; ++k) {
	gsl_vector_set(gradient, k, gsl_vector_get(gradient,k) - 
		       gsl_vector_get(dP,k) / P_iK(Never,0.0,u));
      }
      gsl_vector_free(dP);
      for (size_t k=0; k<ncoef; ++k)
	gsl_vector_set(gradient,k, -gsl_vector_get(gradient,k)*freq);
      return gradient;
    }
    gsl_vector *
    negll_gradient(bool clear_values = true) {
      // clear the stored values
      if (clear_values) {
	P_ijs.clear();
	dP_ijs_clear();
      }
      gsl_vector * total;
      total = gsl_vector_calloc(ncoef);
      for (size_t i = 0; i< finalState->size; ++i) {
	gsl_vector * component = 
	  negll_gradient_component(gsl_vector_int_get(finalState,i), 
				   gsl_vector_get(time1,i), 
				   gsl_vector_get(time2,i), 
				   gsl_vector_get(time3,i), 
				   gsl_vector_get(freq,i), 
				   gsl_vector_int_get(recall,i));
	gsl_vector_add(total,component);
	gsl_vector_free(component);
      }
      return total;
    }
    double negll(bool clear_values = true) {
      if (clear_values) P_ijs.clear();
      double total = 0.0;
      for (size_t i = 0; i< finalState->size; ++i) {
	total += negll_component(gsl_vector_int_get(finalState,i), 
				 gsl_vector_get(time1,i), 
				 gsl_vector_get(time2,i), 
				 gsl_vector_get(time3,i), 
				 gsl_vector_get(freq,i), 
				 gsl_vector_int_get(recall,i));
	if (debug) Rprintf("i=%i,negll=%g\n",i,total);
      }
      if (debug) {
	Rprintf("Size of Pij map = %i\n",P_ijs.size());
	Rprintf("ll(Never,0.0,0.0,70.0)=%f\n",-negll_component(Never,0.0,0.0,70.0));
	Rprintf("ll(Current,20.0,0.0,70.0)=%f\n",-negll_component(Current,20.0,0.0,70.0));
	Rprintf("ll(Former,20.0,50.0,70.0)=%f\n",-negll_component(Former,20.0,50.0,70.0));
	Rprintf("P_ij(Never,Never,0,50)=%f\n",P_ij(Never, Never, 0.0, 50.0));
	Rprintf("P_ij(Never,Current,0,50)=%f\n",P_ij(Never, Current, 0.0, 50.0));
	Rprintf("P_ij(Never,Former,0,50)=%f\n",P_ij(Never, Former, 0.0, 50.0));
	Rprintf("P_ij(Never,Reclassified,0,50)=%f\n",P_ij(Never, Reclassified, 0.0, 20.0));
	Rprintf("P_ij(Never,Death,0,50)=%f\n",1.0 - (P_ij(Never, Never, 0.0, 50.0)+P_ij(Never, Current, 0.0, 50.0)+
						     P_ij(Never, Former, 0.0, 50.0)+P_ij(Never, Reclassified, 0.0, 50.0))); // Death not included in (live) state-space
	Rprintf("mu_ij(Never,Current,50)=%f\n",mu_ij(Never, Current, 50.0));
	Rprintf("mu_ij(Never,Death,50)=%f\n",mu_ij(Never, Death, 50.0));
	Rprintf("mu_ij(Current,Former,50)=%f\n",mu_ij(Current, Former, 50.0));
	Rprintf("mu_ij(Current,Death,50)=%f\n",mu_ij(Current, Death, 50.0));
	Rprintf("mu_ij(Former,Reclassified,50)=%f\n",mu_ij(Former, Reclassified, 50.0));
	Rprintf("mu_ij(Former,Death,50)=%f\n",mu_ij(Former, Death, 50.0));
	Rprintf("mu_ij(Reclassified,Death,50)=%f\n",mu_ij(Reclassified, Death, 50.0));
	Rprintf("P_iK(Never,0.0,70.0)=%f\n",P_iK(Never,0.0,70.0));
	Rprintf("negll(Current,20,70)=%f\n",negll_component(Current, 20.0, 0.0, 70.0));
	Rprintf("negll(Former,20,50,70)=%f\n",negll_component(Former, 20.0, 50.0, 70.0));
	Rprintf("P_ij(Current,Current,9,71)=%f\n",P_ij(Current,Current, 9.0, 71.0));
	Rprintf("P_ij(Current,Former,9,71)=%f\n",P_ij(Current,Former, 9.0, 71.0));
	Rprintf("P_ij(Current,Reclassified,9,71)=%f\n",P_ij(Current,Reclassified, 9.0, 71.0));
      }
      return total;
    }
  };

  class PurgedNS : public Purged {
  public:
    PurgedNS(SEXP sexp) : 
      Purged(sexp) 
    {
      List args = as<List>(sexp);
      gsl_vector * knots01 = as_vector(args("knots01"));
      gsl_vector * knots12 = as_vector(args("knots12"));
      s01 = new ns(knots01, 
		   false, // intercept
		   true, // centred
		   20.0 // centre
		   );
      // knots for cessation
      s12 = new ns(knots12,
		   false, // intercept
		   true, // centred
		   40.0 // centre
		   );
      gsl_vector_free(knots01); gsl_vector_free(knots12);
    }
    ~PurgedNS() {
      delete s01; 
      delete s12;
    }
  };

  class PurgedPS : public Purged {
  public:
    gsl_matrix *pmatrix01, *pmatrix12;
    double sp01, sp12;
    PurgedPS(SEXP sexp) : Purged(sexp) {
      List args = as<List>(sexp);
      pmatrix01 = as_matrix(args("pmatrix01"));
      pmatrix12 = as_matrix(args("pmatrix12"));
      sp01 = args("sp01");
      sp12 = args("sp12");
      s01 = new ps(as<double>(args("lower01")), // lower boundary
		   as<double>(args("upper01")), // upper boundary
		   false, // intercept
		   true, // centred
		   20.0, // centre
		   as<int>(args("nterm01")) // = nbeta01-2
		   );
      // knots for cessation
      s12 = new ps(as<double>(args("lower12")),
		   as<double>(args("upper12")),
		   false, // intercept
		   true, // centred
		   40.0, // centre
		   as<int>(args("nterm12")) // = nbeta12-2
		   );
    }
    double pnegll() {
      double pnegll = Purged::negll();
      gsl_vector *v01, *v12;
      double penalty01=0.0, penalty12=0.0;
      v01 = gsl_vector_alloc(beta01->size);
      gsl_blas_dgemv(CblasNoTrans, 1.0, pmatrix01, beta01, 0.0, v01);
      gsl_blas_ddot(beta01,v01,&penalty01);
      pnegll += penalty01 * sp01/2.0; // add positive penalty
      v12 = gsl_vector_alloc(beta12->size);
      gsl_blas_dgemv(CblasNoTrans, 1.0, pmatrix12, beta12, 0.0, v12);
      gsl_blas_ddot(beta12,v12,&penalty12);
      pnegll += penalty12 * sp12/2.0; // add positive penalty
      // tidy up
      gsl_vector_free(v01);
      gsl_vector_free(v12);
      return pnegll;
    }
    gsl_vector * pnegll_gradient() {
      gsl_vector * pgrad = negll_gradient();
      gsl_vector *v01, *v12;
      v01 = gsl_vector_alloc(beta01->size);
      gsl_blas_dgemv(CblasNoTrans, 1.0, pmatrix01, beta01, 0.0, v01);
      int j=1; // ignore intercept term
      for (size_t i=0; i<nbeta01; ++i, ++j)
	gsl_vector_set(pgrad, j, gsl_vector_get(pgrad,j) + sp01*gsl_vector_get(v01,i)); // add positive penalty
      j++; // ignore next intercept term
      v12 = gsl_vector_alloc(beta12->size);
      gsl_blas_dgemv(CblasNoTrans, 1.0, pmatrix12, beta12, 0.0, v12);
      for (size_t i=0; i<nbeta12; ++i, ++j)
	gsl_vector_set(pgrad, j, 
		       gsl_vector_get(pgrad,j) + sp12*gsl_vector_get(v12,i)); // add positive penalty
      // tidy up
      gsl_vector_free(v01);
      gsl_vector_free(v12);
      return pgrad;
    }
    virtual ~PurgedPS() {
      delete s01; delete s12;
      gsl_matrix_free(pmatrix01); gsl_matrix_free(pmatrix12);
    }
  };

  // define the ODE  (Model B, using observed mortality)
  int 
  funcReclassifiedB(double t, const double y[], double f[],
		    void *model)
  {
    double alpha01, alpha12, alpha20, mu0, mu1, mu2;
    Purged* P = static_cast<Purged*>(model);
    alpha01 = exp(P->int01 + P->s01->calc(t, P->beta01));
    alpha12 = exp(P->int12 + P->s12->calc(t, P->beta12));
    alpha20 = exp(P->beta20);
    mu0 = gsl_spline_eval(P->spline0, bounds(t,0.5,P->maxage), P->acc0);
    mu1 = gsl_spline_eval(P->spline1, bounds(t,0.5,P->maxage), P->acc1);
    mu2 = gsl_spline_eval(P->spline2, bounds(t,0.5,P->maxage), P->acc2);
    f[Never] = -(alpha01 + mu0)*y[Never];
    f[Current] = alpha01*y[Never]-(alpha12+mu1)*y[Current];
    f[Former] = alpha12*y[Current]-(alpha20+mu2)*y[Former];
    f[Reclassified] = alpha20*y[Former]-mu0*y[Reclassified];
    return GSL_SUCCESS;
  }
  // define the ODE (excess mortality)
  int
  funcReclassified (double t, const double y[], double f[],
		    void *model)
  {
    double alpha01, alpha12, alpha20, mu0, mu1, mu2;
    Purged* P = static_cast<Purged*>(model);
    double eps = 1e-16;
    alpha01 = exp(P->int01 + P->s01->calc(t, P->beta01));
    alpha12 = exp(P->int12 + P->s12->calc(t, P->beta12));
    alpha20 = exp(P->beta20);
    mu0 = gsl_spline_eval(P->spline0, bounds(t,0.5,P->maxage), P->acc0);
    mu1 = gsl_spline_eval(P->spline1, bounds(t,0.5,P->maxage), P->acc1);
    mu2 = gsl_spline_eval(P->spline2, bounds(t,0.5,P->maxage), P->acc2);
    f[Never] = -(alpha01)*y[Never];
    f[Current] = alpha01*y[Never]-(alpha12+std::max(eps,mu1-mu0))*y[Current];
    f[Former] = alpha12*y[Current]-(alpha20+std::max(eps,mu2-mu0))*y[Former];
    f[Reclassified] = alpha20*y[Former];
    return GSL_SUCCESS;
  }

  /** @brief
      Implementation that also calculates partials for P_ij(s,t)
   **/
  int
  funcReclassifiedExtended (double t, 
			    const double y[],   // dimension: nstate*(ncoef+1)
			    double f[],         // dimension: nstate*(ncoef+1)
			    void *model)       // type: Param*
  {
    double alpha01, alpha12, alpha20, mu0, mu1, mu2;
    double dalpha01, dalpha12, dalpha20;
    Purged* P = static_cast<Purged*>(model);
    gsl_vector *x01; x01 = gsl_vector_calloc(P->ncoef);
    gsl_vector *x12; x12 = gsl_vector_calloc(P->ncoef);
    gsl_vector *x20; x20 = gsl_vector_calloc(P->ncoef);
    gsl_vector *basis01; basis01 = gsl_vector_alloc(P->nbeta01);
    gsl_vector *basis12; basis12 = gsl_vector_alloc(P->nbeta12);
    int j = 0;
    double eps = 1e-16;
    // calculate design vectors for the splines
    P->s01->calcBasis(t, basis01);
    P->s12->calcBasis(t, basis12);
    // set x 
    // intercept for 0->1
    gsl_vector_set(x01, j++, 1);
    // spline parameters for 0->1
    for (size_t i = 0; i < P->nbeta01; ++i) {
      gsl_vector_set(x01, j++,  gsl_vector_get(basis01,i)); // unclear why the gsl_vector_get() call is needed
    }
    // intercept for 1->2
    gsl_vector_set(x12, j++, 1);
    // spline parameters for 1->2
    for (size_t i = 0; i < P->nbeta12; ++i) {
      gsl_vector_set(x12, j++, gsl_vector_get(basis12,i));
    }
    // intercept for 2->0
    gsl_vector_set(x20, j++, 1);
    // calculate rates
    alpha01 = P->mu_ij(Never,Current,t);
    alpha12 = P->mu_ij(Current,Former,t);
    alpha20 = P->mu_ij(Former,Reclassified,t);
    mu0 = P->mu_ij(Never,Death,t);
    mu1 = P->mu_ij(Current,Death,t);
    mu2 = P->mu_ij(Former,Death,t);
    // base equations
    f[Never] = -alpha01*y[Never];
    f[Current] = alpha01*y[Never]-(alpha12+std::max(eps,mu1-mu0))*y[Current];
    f[Former] = alpha12*y[Current]-(alpha20+std::max(eps,mu2-mu0))*y[Former];
    f[Reclassified] = alpha20*y[Former];
    // d/dt (partial P(s,t)/partial theta) = (partial P(s,t)/partial theta)*alpha + P(s,t)(partial alpha / partial theta)
    for (int offset=P->nstate, i=0; i < (int) P->ncoef; ++i, offset += P->nstate) {
      dalpha01 = alpha01 * gsl_vector_get(x01,i);
      dalpha12 = alpha12 * gsl_vector_get(x12,i);
      dalpha20 = alpha20 * gsl_vector_get(x20,i);
      f[Never+offset] = - (alpha01*y[Never+offset] + dalpha01*y[Never]);
      f[Current+offset] = alpha01*y[Never+offset] + dalpha01*y[Never] -
	((alpha12+std::max(eps,mu1-mu0))*y[Current+offset] + dalpha12*y[Current]);
      f[Former+offset] = alpha12*y[Current+offset]+dalpha12*y[Current] -
	((alpha20+std::max(eps,mu2-mu0))*y[Former+offset] + dalpha20*y[Former]);
      f[Reclassified+offset] = alpha20*y[Former+offset] + dalpha20*y[Former];
    }
    // tidy up
    gsl_vector_free(x01); gsl_vector_free(x12); gsl_vector_free(x20);
    gsl_vector_free(basis01); gsl_vector_free(basis12);
    return GSL_SUCCESS;
  }

  RcppExport
  SEXP call_purged_ns(SEXP sexp) {
    List args = as<List>(sexp);
    // bool debug = as<bool>(args("debug"));
    PurgedNS model(sexp);
    std::string output_type = args("output_type");
    if (output_type == "negll")
      return wrap(model.negll());
    if (output_type == "negll_gradient") {
      gsl_vector * grad = model.negll_gradient();
      SEXP value = wrap_vector(grad);
      gsl_vector_free(grad);
      return value;
    }
    if (output_type == "P_ij")
      return wrap(model.P_ij(Never,Reclassified,0.0, 50.0));
    if (output_type == "dP_ij") {
      gsl_vector * dP = model.dP_ij(Never,Reclassified,0.0, 50.0);
      SEXP value = wrap_vector(dP);
      gsl_vector_free(dP);
      return value;
    }
    Rprintf("call_purged_ns: output_type not matched (\"%s\").\n",output_type.c_str());
    return wrap(true); // default
  }

  RcppExport 
  SEXP test_ns() {
    gsl_vector * knots;
    knots = gsl_vector_alloc(2);
    gsl_vector_set(knots,0,1.0);
    gsl_vector_set(knots,1,2.0);
    ns ns1(knots);
    gsl_vector_free(knots);
    return wrap(true);
  }
  RcppExport 
  SEXP test_ps() {
    ps ps1(0.0,1.0);
    return wrap(true);
  }

  RcppExport
  SEXP call_purged_ps(SEXP sexp) {
    List args = as<List>(sexp);
    // bool debug = as<bool>(args("debug"));
    PurgedPS model(sexp);
    std::string output_type = args("output_type");
    if (output_type == "negll")
      return wrap(model.negll());
    if (output_type == "negll_gradient") {
      gsl_vector * grad = model.negll_gradient();
      SEXP value = wrap_vector(grad);
      gsl_vector_free(grad);
      return value;
    }
    if (output_type == "pnegll_gradient") {
      gsl_vector * grad = model.pnegll_gradient();
      SEXP value = wrap_vector(grad);
      gsl_vector_free(grad);
      return value;
    }
    if (output_type == "pnegll")
      return wrap(model.pnegll());
    if (output_type == "P_ij")
      return wrap(model.P_ij(Never,Reclassified,0.0, 50.0));
    if (output_type == "dP_ij") {
      gsl_vector * dP = model.dP_ij(Never,Reclassified,0.0, 50.0);
      SEXP value = wrap_vector(dP);
      gsl_vector_free(dP);
      return value;
    }
    Rprintf("call_purged_ps: output_type not matched (\"%s\").\n",output_type.c_str());
    return wrap(true); // default
  }

}
