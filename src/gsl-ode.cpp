#include <Rcpp.h>
#include <RcppGSL.h>

#include <cstdlib>
#include <map>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <boost/algorithm/cxx11/iota.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

namespace {

using namespace Rcpp;


//forward declaration

class ns;


// declare types

  // enum's cause problems with RcppGSL::vector, so instead we use int's 
  int Never=0, Current=1, Former=2, Reclassified=3, Death=3;

struct Param {
  double int01, int12, beta20, RRs, RRx;
  gsl_spline *spline0;
  gsl_interp_accel *acc0;
  ns *ns01;
  ns *ns12;
  gsl_vector *beta01;
  gsl_vector *beta12;
};

class ns {
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
      calcBasis(centre, nsC, true);
    }
    gsl_bspline_deriv_free(dw);
    gsl_matrix_free(dB);
    gsl_matrix_free(const_matrix);
    gsl_matrix_free(R);
    gsl_matrix_free(Q_full);
    gsl_vector_free(tau);
  }
  ~ns() {
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
  void calcBasis(double x, gsl_vector * basis, bool nocentering = false) {
    size_t i, j;
    size_t n = 1;
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
    if (centred && !nocentering)
      gsl_blas_daxpy (-1.0, nsC, basis);
  }
  double calc(double x, gsl_vector * beta, bool nocentering = false) {
    double value;
    calcBasis(x, _basis, nocentering);
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

class ps {
public:
  ps(double boundaryL, double boundaryU, bool intercept = false, bool centred = false, double centre = 0.0, size_t nterm = 10, size_t k = 4) : 
    intercept(intercept), boundaryL(boundaryL), boundaryU(boundaryU), nterm(nterm), centred(centred) {
    size_t i, j, degree;
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
    for (i = 0; i < ncoeffs_base; ++i) 
      {
	gsl_matrix_set(ttl, 0, i, gsl_matrix_get(dB, i, 0));
	gsl_matrix_set(ttl, 1, i, gsl_matrix_get(dB, i, 1));
      }
    gsl_bspline_deriv_eval (boundaryU, 1, dB, bw, dw);
    for (i = 0; i < ncoeffs_base; ++i) 
      {
	gsl_matrix_set(ttu, 0, i, gsl_matrix_get(dB, i, 0));
	gsl_matrix_set(ttu, 1, i, gsl_matrix_get(dB, i, 1));
      }
    if (centred) {
      calcBasis(centre, nsC, true);
    }
    gsl_bspline_deriv_free(dw);
    gsl_matrix_free(dB);
  }
  ~ps() {
    gsl_bspline_free(bw);
    gsl_matrix_free(ttl);
    gsl_matrix_free(ttu);
    gsl_vector_free(B);
    gsl_vector_free(xv);
    gsl_vector_free(knots);
    gsl_vector_free(nsC);
    gsl_vector_free(_basis);
  }
  void calcBasis(double x, gsl_vector * basis, bool nocentering = false) {
    size_t i, j;
    size_t n = 1;
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
    for (j = 0; j < ncoeffs; ++j)
      {
	gsl_vector_set(basis, j, gsl_vector_get(B, j+3+offset));
      }
    if (centred && !nocentering)
      gsl_blas_daxpy (-1.0, nsC, basis);
  }
  double calc(double x, gsl_vector * beta, bool nocentering = false) {
    double value;
    calcBasis(x, _basis, nocentering);
    int status = gsl_blas_ddot(beta, _basis, &value);
    if (status != GSL_SUCCESS)
      REprintf ("error, return value=%d\n", status);
    return value;
  }
  bool intercept, centred;
  size_t nbreak, ncoeffs_base, ncoeffs_int, ncoeffs, nterm;
  int offset;
  double boundaryL, boundaryU;
  gsl_bspline_workspace *bw;
  gsl_matrix *ttl, *ttu;
  gsl_vector *knots, *B, *xv, *nsC, *_basis;
};


// utility functions

void gsl_matrix_printf(gsl_matrix * mat) {
  int i, j;
  for (i=0; i < mat->size1; ++i) {
    for (j=0; j < mat->size2; ++j) 
      Rprintf("%f ", gsl_matrix_get(mat, i, j));
    Rprintf("\n");
  }
}

void gsl_vector_printf(gsl_vector * vec) {
  int i;
  for (i=0; i < vec->size; ++i) {
      Rprintf("%f ", gsl_vector_get(vec, i));
  }
}

// template<class T> T bounds(T x, T lo, T hi) { return x<lo ? lo : (x>hi ? hi : x); } 
double bounds(double x, double lo, double hi) { return x<lo ? lo : (x>hi ? hi : x); } 


// define the ODE
int
funcCycle (double t, const double y[], double f[],
      void *params)
{
  double alpha01, alpha12, alpha20, mu;
  Param* P = static_cast<Param*>(params);
  alpha01 = exp(P->int01 + P->ns01->calc(t, P->beta01));
  alpha12 = exp(P->int12 + P->ns12->calc(t, P->beta12));
  alpha20 = exp(P->beta20);
  mu = gsl_spline_eval(P->spline0, bounds(t,0.5,105.5), P->acc0);
  f[Never] = -(alpha01 + mu)*y[Never]+alpha20*y[Former]; // uptake for reclassified
  f[Current] = alpha01*y[Never]-(alpha12+P->RRs*mu)*y[Current];
  f[Former] = alpha12*y[Current]-(alpha20+P->RRx*mu)*y[Former];
  return GSL_SUCCESS;
}

double
mu_ijCycle (int i, int j, double t, gsl_odeiv2_driver * d)
{
  Param* P = static_cast<Param*>(d->sys->params);
  if (i==Never & j==Current)  return exp(P->int01 + P->ns01->calc(t, P->beta01));
  if (i==Current & j==Former) return exp(P->int12 + P->ns12->calc(t, P->beta12));
  if (i==Former & j==Never)  return exp(P->beta20);
  if (j==Death) {
    double mu = gsl_spline_eval(P->spline0, bounds(t,0.5,105.5), P->acc0);
    if (i==Never) return mu;
    if (i==Current) return P->RRs*mu;
    if (i==Former) return P->RRx*mu;
    }    
  return 0.0;
}

int
funcReclassified (double t, const double y[], double f[],
      void *params)
{
  double alpha01, alpha12, alpha20, mu;
  Param* P = static_cast<Param*>(params);
  alpha01 = exp(P->int01 + P->ns01->calc(t, P->beta01));
  alpha12 = exp(P->int12 + P->ns12->calc(t, P->beta12));
  alpha20 = exp(P->beta20);
  mu = gsl_spline_eval(P->spline0, bounds(t,0.5,105.5), P->acc0);
  f[Never] = -(alpha01 + mu)*y[Never];
  f[Current] = alpha01*y[Never]-(alpha12+P->RRs*mu)*y[Current];
  f[Former] = alpha12*y[Current]-(alpha20+P->RRx*mu)*y[Former];
  f[Reclassified] = -mu*y[Reclassified]+alpha20*y[Former]; // no uptake for reclassified
  return GSL_SUCCESS;
}

double
mu_ijReclassified (int i, int j, double t, gsl_odeiv2_driver * d)
{
  Param* P = static_cast<Param*>(d->sys->params);
  if (i==Never & j==Current)  return exp(P->int01 + P->ns01->calc(t, P->beta01));
  if (i==Current & j==Former) return exp(P->int12 + P->ns12->calc(t, P->beta12));
  if (i==Former & j==Reclassified)  return exp(P->beta20);
  if (j==Death) {
    double mu = gsl_spline_eval(P->spline0, bounds(t,0.5,105.5), P->acc0);
    if (i==Never) return mu;
    if (i==Current) return P->RRs*mu;
    if (i==Former) return P->RRx*mu;
    if (i==Reclassified) return mu;
  }
  REprintf("Unexpected combinations of states: %i to %j.\n",i,j);
  return 0.0;
}

// for efficiency, cache the values for P_ij
typedef boost::tuple<int,int,double,double> P_ij_key;
std::map<P_ij_key,double> P_ijs;

double
P_ij(int i, int j, double s, double t, gsl_odeiv2_driver * d) {
  // do we have the values cached?
  P_ij_key key = P_ij_key(i,j,s,t);
  if (P_ijs.count(key)>0) return P_ijs[key];
  // otherwise...
  size_t dimension = d->sys->dimension;
  double y[dimension];
  for (int k=0; k<dimension; ++k)
      y[k] = 0.0;
  y[i] = 1.0;
  double tm = s;
  int status = gsl_odeiv2_driver_apply (d, &tm, t, y);
  if (status != GSL_SUCCESS)
      REprintf ("error, return value=%d\n", status);
  // cache the results
  for (int k=0; k<dimension; ++k)
    P_ijs[P_ij_key(i,k,s,t)] = y[k];
  return y[j];
}

double
P_iK(int i, double s, double t, gsl_odeiv2_driver * d) {
  size_t dimension = d->sys->dimension;
  double total = 0.0;
  for (int k=0; k<dimension; ++k)
    total += P_ij(i,k,s,t,d);
  return total;
}

  // Note: does this likelihood ignore cycles aside from Never -> Never ?
double negllCycle(int state, double s, double t, double u, gsl_odeiv2_driver * d) {
  double ll;
  if (state == Never)
    ll = log(P_ij(Never,Never,0.0,u,d))-log(P_iK(Never,0.0,u,d)); // ignore s and t
  if (state == Current)
    ll = log(P_ij(Never,Never,0.0,s,d))+log(mu_ij(Never,Current,s,d))+log(P_ij(Current,Current,s,u,d))-
      log(P_iK(Never,0.0,u,d)); // ignore t
  if (state == Former)
    ll = log(P_ij(Never,Never,0.0,s,d))+log(mu_ij(Never,Current,s,d))+log(P_ij(Current,Current,s,t,d))+
      log(mu_ij(Current,Former,t,d))+log(P_ij(Former,Former,t,u,d))-log(P_iK(Never,0.0,u,d));
  return -ll;
}

  // the required mu_ij are not affected by the additional reclassification state
double negllReclassified(int state, double s, double t, double u, gsl_odeiv2_driver * d) {
  double ll;
  if (state == Never) // Never->Never _and_ Never->Reclassified
    ll = log(P_ij(Never,Never,0.0,u,d)+P_ij(Never,Reclassified,0.0,u,d))-log(P_iK(Never,0.0,u,d)); // ignore s and t
  if (state == Current)
    ll = log(P_ij(Never,Never,0.0,s,d))+log(mu_ij(Never,Current,s,d))+log(P_ij(Current,Current,s,u,d))-
      log(P_iK(Never,0.0,u,d)); // ignore t
  if (state == Former)
    ll = log(P_ij(Never,Never,0.0,s,d))+log(mu_ij(Never,Current,s,d))+log(P_ij(Current,Current,s,t,d))+
      log(mu_ij(Current,Former,t,d))+log(P_ij(Former,Former,t,u,d))-log(P_iK(Never,0.0,u,d));
  return -ll;
}

  RcppExport SEXP gsl_main2Cycle(SEXP _finalState,
			    SEXP _time1, SEXP _time2, SEXP _time3,
			    SEXP _int01, SEXP _int12,
			    SEXP _knots01, SEXP _knots12, 
			    SEXP _beta01, SEXP _beta12,
			    SEXP _beta20,
			    SEXP _ages0,
			    SEXP _mu0)
  {
    
    RcppGSL::vector<double> 
      finalState = _finalState,
      time1=_time1,
      time2=_time2,
      time3=_time3,
      knots01=_knots01,
      knots12=_knots12,
      beta01 = _beta01,
      beta12 = _beta12,
      ages0 = _ages0,
      mu0 = _mu0;

    double int01=as<double>(_int01), 
      int12=as<double>(_int12),
      beta20=as<double>(_beta20);
  
  gsl_spline *spline0;
  gsl_interp_accel *acc0;
  acc0 = gsl_interp_accel_alloc ();
  spline0 = gsl_spline_alloc (gsl_interp_cspline, 106);
  gsl_spline_init (spline0, ages0->data, mu0->data, 106);

  // initial parameters for initiation
  ns * ns01 = new ns(knots01, 
		    false, // intercept
		    true, // centred
		    20.0 // centre
		    );

  // knots for cessation
  ns * ns12 = new ns(knots12, 
		    false, // intercept
		    true, // centred
		    40.0 // centre
		    );

  Param beta = {int01, int12, beta20, 1.4, 1.2, spline0, acc0, ns01, ns12, beta01, beta12};

  gsl_odeiv2_system sys = {funcCycle, NULL, 3, &beta};
  gsl_odeiv2_driver * d = 
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
				  1e-8, 1e-8, 0.0);

  // clear the values whenever the parameters change
  P_ijs.clear();
  // without storing intermediate results, 10000 ODEs took ~20 sec

  double sum_negll = 0.0;
  for (int i = 0; i< finalState.size(); ++i) 
    sum_negll += negllCycle(finalState[i], time1[i], time2[i], time3[i], d);
  //sum_negll += negllCycle(Former, time1[0], time2[0], time3[0],d);

  // Rprintf("Size of Pij map = %i\n",P_ijs.size());
  // Rprintf("ll(Never,0.0,0.0,70.0)=%f\n",-negllCycle(Never,0.0,0.0,70.0,d));
  // Rprintf("ll(Never,0.0,0.0,70.0)=%f\n",-negllCycle(Never,0.0,0.0,70.0,d));
  // Rprintf("ll(Current,20.0,0.0,70.0)=%f\n",-negllCycle(Current,20.0,0.0,70.0,d));
  // Rprintf("ll(Former,20.0,50.0,70.0)=%f\n",-negllCycle(Former,20.0,50.0,70.0,d));

  delete ns01;
  delete ns12;

  finalState.free();
  time1.free();
  time2.free();
  time3.free();
  knots01.free();
  knots12.free();
  beta01.free();
  beta12.free();
  ages0.free();
  mu0.free();
  
  gsl_spline_free (spline0);
  gsl_interp_accel_free (acc0);
  gsl_odeiv2_driver_free (d);

  return wrap(sum_negll);
  }



  RcppExport SEXP gsl_main2Reclassified(SEXP _finalState,
			    SEXP _time1, SEXP _time2, SEXP _time3,
			    SEXP _int01, SEXP _int12,
			    SEXP _knots01, SEXP _knots12, 
			    SEXP _beta01, SEXP _beta12,
			    SEXP _beta20,
			    SEXP _ages0,
			    SEXP _mu0)
  {
    
    RcppGSL::vector<double> 
      finalState = _finalState,
      time1=_time1,
      time2=_time2,
      time3=_time3,
      knots01=_knots01,
      knots12=_knots12,
      beta01 = _beta01,
      beta12 = _beta12,
      ages0 = _ages0,
      mu0 = _mu0;

    double int01=as<double>(_int01), 
      int12=as<double>(_int12),
      beta20=as<double>(_beta20);
  
  gsl_spline *spline0;
  gsl_interp_accel *acc0;
  acc0 = gsl_interp_accel_alloc ();
  spline0 = gsl_spline_alloc (gsl_interp_cspline, 106);
  gsl_spline_init (spline0, ages0->data, mu0->data, 106);

  // initial parameters for initiation
  ns * ns01 = new ns(knots01, 
		    false, // intercept
		    true, // centred
		    20.0 // centre
		    );

  // knots for cessation
  ns * ns12 = new ns(knots12, 
		    false, // intercept
		    true, // centred
		    40.0 // centre
		    );

  Param beta = {int01, int12, beta20, 1.4, 1.2, spline0, acc0, ns01, ns12, beta01, beta12};

  gsl_odeiv2_system sys = {funcReclassified, NULL, 4, &beta}; // cf (func, NULL, 3, &beta)
  gsl_odeiv2_driver * d = 
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
				  1e-8, 1e-8, 0.0);

  // clear the values whenever the parameters change
  P_ijs.clear();
  // without storing intermediate results, 10000 ODEs took ~20 sec

  double sum_negll = 0.0;
  for (int i = 0; i< finalState.size(); ++i) 
     sum_negll += negllReclassified(finalState[i], time1[i], time2[i], time3[i], d);

  // Rprintf("Size of Pij map = %i\n",P_ijs.size());
  // Rprintf("P(Never,Former,20.0,50.0,70.0,d)=%f\n",P_ij(Never,Former,20.0,50.0,70.0,d));

  // Rprintf("ll(Never,0.0,0.0,70.0)=%f\n",-negllReclassified(Never,0.0,0.0,70.0,d));
  // Rprintf("ll(Never,0.0,0.0,70.0)=%f\n",-negllReclassified(Never,0.0,0.0,70.0,d));
  // Rprintf("ll(Current,20.0,0.0,70.0)=%f\n",-negllReclassified(Current,20.0,0.0,70.0,d));
  // Rprintf("ll(Former,20.0,50.0,70.0)=%f\n",-negllReclassified(Former,20.0,50.0,70.0,d));

  // Rprintf("P(Never,70)=%f\n",negllReclassified(Never, 0.0, 0.0, 70.0,d));
  // Rprintf("P_ij(Never,Never,0,50)=%f\n",P_ij(Current, Current, 0.0, 20.0, d));
  // Rprintf("mu_ij(Never,Current,20)=%f\n",mu_ij(Never, Current, 20.0, d));
  // Rprintf("P_ij(Current,Current,20,70)=%f\n",P_ij(Current, Current, 20.0, 70.0, d));
  // Rprintf("P_iK(Never,0.0,70.0)=%f\n",P_iK(Never,0.0,70.0,d));
  
  // Rprintf("negll(Current,20,70)=%f\n",negllReclassified(Current, 20.0, 0.0, 70.0,d));
  // Rprintf("negll(Former,20,50,70)=%f\n",negllReclassified(Former, 20.0, 50.0, 70.0,d));

  delete ns01;
  delete ns12;

  finalState.free();
  time1.free();
  time2.free();
  time3.free();
  knots01.free();
  knots12.free();
  beta01.free();
  beta12.free();
  ages0.free();
  mu0.free();
  
  gsl_spline_free (spline0);
  gsl_interp_accel_free (acc0);
  gsl_odeiv2_driver_free (d);

  return wrap(sum_negll);
  }
    
//   // old development code - deprecated
//   RcppExport SEXP gsl_main() {

//   // input data
//   double mu0[] = {0.00219, 0.000304, 5.2e-05, 0.000139, 0.000141, 3.6e-05, 7.3e-05, 
// 		  0.000129, 3.8e-05, 0.000137, 6e-05, 8.1e-05, 6.1e-05, 0.00012, 
// 		0.000117, 0.000183, 0.000185, 0.000397, 0.000394, 0.000585, 0.000448, 
// 		0.000696, 0.000611, 0.000708, 0.000659, 0.000643, 0.000654, 0.000651, 
// 		0.000687, 0.000637, 0.00063, 0.000892, 0.000543, 0.00058, 0.00077, 
// 		0.000702, 0.000768, 0.000664, 0.000787, 0.00081, 0.000991, 9e-04, 
// 		0.000933, 0.001229, 0.001633, 0.001396, 0.001673, 0.001926, 0.002217, 
// 		0.002562, 0.002648, 0.002949, 0.002729, 0.003415, 0.003694, 0.004491, 
// 		0.00506, 0.004568, 0.006163, 0.006988, 0.006744, 0.00765, 0.007914, 
// 		0.009153, 0.010231, 0.011971, 0.013092, 0.013839, 0.015995, 0.017693, 
// 		0.018548, 0.020708, 0.022404, 0.02572, 0.028039, 0.031564, 0.038182, 
// 		0.042057, 0.047361, 0.05315, 0.058238, 0.062619, 0.074934, 0.089776, 
// 		0.099887, 0.112347, 0.125351, 0.143077, 0.153189, 0.179702, 0.198436, 
// 		0.240339, 0.256215, 0.275103, 0.314157, 0.345252, 0.359275, 0.41768, 
// 		0.430279, 0.463636, 0.491275, 0.549738, 0.354545, 0.553846, 0.461538, 
// 		0.782609};
//   double ages0[106];
//   gsl_spline *spline0;
//   gsl_interp_accel *acc0;
//   boost::algorithm::iota(ages0, ages0+106, 0.5);
//   acc0 = gsl_interp_accel_alloc ();
//   spline0 = gsl_spline_alloc (gsl_interp_cspline, 106);
//   gsl_spline_init (spline0, ages0, mu0, 106);

//   gsl_vector *knots01, *beta01;
//   gsl_vector *knots12, *beta12;

//   // knots for initiation
//   knots01 = gsl_vector_alloc(3);
//   gsl_vector_set(knots01, 0, 10.0);
//   gsl_vector_set(knots01, 1, 20.0);
//   gsl_vector_set(knots01, 2, 30.0);
//   // initial parameters for initiation
//   ns * ns01 = new ns(knots01, 
// 		    false, // intercept
// 		    true, // centred
// 		    20.0 // centre
// 		    );
//   beta01 = gsl_vector_alloc(ns01->ncoeffs);
//   gsl_vector_set(beta01, 0, 1.0);
//   gsl_vector_set(beta01, 1, -1.0);

//   // knots for cessation
//   knots12 = gsl_vector_alloc(3);
//   gsl_vector_set(knots12, 0, 20.0);
//   gsl_vector_set(knots12, 1, 40.0);
//   gsl_vector_set(knots12, 2, 60.0);
//   ns * ns12 = new ns(knots12, 
// 		    false, // intercept
// 		    true, // centred
// 		    40.0 // centre
// 		    );
//   beta12 = gsl_vector_alloc(ns12->ncoeffs);
//   gsl_vector_set(beta12, 0, 1.0);
//   gsl_vector_set(beta12, 1, -1.0);

//   Param beta = {-3.0, -4.0, log(0.01), 1.4, 1.2, spline0, acc0, ns01, ns12, beta01, beta12};

//   gsl_odeiv2_system sys = {func, NULL, 3, &beta};
//   gsl_odeiv2_driver * d = 
//     gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
// 				  1e-8, 1e-8, 0.0);

//   // clear the values whenever the parameters change
//   P_ijs.clear();
//   // without storing intermediate results, 10000 ODEs took ~20 sec


//   Rprintf("ll(Never,0.0,0.0,70.0)=%f\n",-negll(Never,0.0,0.0,70.0,d));
//   Rprintf("ll(Never,0.0,0.0,70.0)=%f\n",-negll(Never,0.0,0.0,70.0,d));
//   Rprintf("ll(Current,20.0,0.0,70.0)=%f\n",-negll(Current,20.0,0.0,70.0,d));
//   Rprintf("ll(Former,20.0,50.0,70.0)=%f\n",-negll(Former,20.0,50.0,70.0,d));

//   // Rprintf("alpha01(10)=%f\n",mu_ij(Never,Current,10.0,d));
//   // Rprintf("alpha12(10)=%f\n",mu_ij(Current,Former,10.0,d));
//   // Rprintf("alpha20(10)=%f\n",mu_ij(Former,Never,10.0,d));
//   // Rprintf("m0(10)=%f\n",mu_ij(Never,Death,10.0,d));

//   // for (double age=10.0; age<100.0; age += 10.0) {
//   //   double never = P_ij(Never,Never,0.0,age,d);
//   //   double current = P_ij(Never,Current,0.0,age,d);
//   //   double former = P_ij(Never,Former,0.0,age,d);
//   //   double total = never + current + former;
//   //   Rprintf("Age=%f, pNever=%f, pCurrent=%f, pFormer=%f\n",age,never/total,current/total,former/total);
//   //   //Rprintf("Age=%f, never=%f, current=%f, former=%f\n",age,never,current,former);
//   // }

//   // for (int i = 0; i<10000; ++i) {
//   //   negll(Former,20.0,50.0,70.0,d);
//   // }

//   // double value = negll(Former,20.0,50.0,70.0,d);

//   delete ns01;
//   delete ns12;

//   gsl_vector_free (beta01);
//   gsl_vector_free (beta12);
//   gsl_spline_free (spline0);
//   gsl_interp_accel_free (acc0);
//   gsl_odeiv2_driver_free (d);

//   return wrap(value);
// }

}
