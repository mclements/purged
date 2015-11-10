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

#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <boost/algorithm/cxx11/iota.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

namespace { // anonymous

  using namespace Rcpp;

  //forward declaration(s)

  // Base class for splines
  class splineBasis {
  public:
    virtual double calc(double t, gsl_vector * beta, bool nocentering = false) = 0;
    virtual void calcBasis(double x, gsl_vector * basis, bool nocentering = false) = 0;
  };
  
  // declare types

  // enum's cause problems with RcppGSL::vector, so instead we use int's 
  int Never=0, Current=1, Former=2, Reclassified=3, Death=4; // Issue: Death=3 or Death=4?
  int CurrentStatus=0,Recall=1,FormerWithCessation=2;
  
  struct Param {
    double int01, int12, beta20, maxage;
    gsl_spline *spline0, *spline1, *spline2;
    gsl_interp_accel *acc0, *acc1, *acc2;
    splineBasis *s01;
    splineBasis *s12;
    gsl_vector *beta01;
    gsl_vector *beta12;
  };

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
	calcBasis(centre, nsC, true);
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
    void calcBasis(double x, gsl_vector * basis, bool nocentering = false) {
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
	calcBasis(centre, nsC, true);
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
    void calcBasis(double x, gsl_vector * basis, bool nocentering = false) {
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
    int nbreak, ncoeffs_base, ncoeffs_int, ncoeffs, nterm;
    int offset;
    double boundaryL, boundaryU;
    gsl_bspline_workspace *bw;
    gsl_matrix *ttl, *ttu;
    gsl_vector *knots, *B, *xv, *nsC, *_basis;
  };


  // utility functions

  void gsl_matrix_Rprintf(gsl_matrix * mat) {
    size_t i, j;
    for (i=0; i < mat->size1; ++i) {
      for (j=0; j < mat->size2; ++j) 
	Rprintf("%f ", gsl_matrix_get(mat, i, j));
      Rprintf("\n");
    }
  }

  void gsl_vector_Rprintf(gsl_vector * vec) {
    size_t i;
    for (i=0; i < vec->size; ++i) {
      Rprintf("%f ", gsl_vector_get(vec, i));
    }
  }

  // template<class T> T bounds(T x, T lo, T hi) { return x<lo ? lo : (x>hi ? hi : x); } 
  inline
  double bounds(double x, double lo, double hi) { return x<lo ? lo : (x>hi ? hi : x); } 


  // define the ODE  (Model B, using observed mortality)
  int
  funcReclassifiedOLD (double t, const double y[], double f[],
		    void *params)
  {
    double alpha01, alpha12, alpha20, mu0, mu1, mu2;
    Param* P = static_cast<Param*>(params);
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

  // define the ODE (Model C, using excess mortality)
  int
  funcReclassified (double t, const double y[], double f[],
		    void *params)
  {
    double alpha01, alpha12, alpha20, mu0, mu1, mu2;
    Param* P = static_cast<Param*>(params);
    alpha01 = exp(P->int01 + P->s01->calc(t, P->beta01));
    alpha12 = exp(P->int12 + P->s12->calc(t, P->beta12));
    alpha20 = exp(P->beta20);
    mu0 = gsl_spline_eval(P->spline0, bounds(t,0.5,P->maxage), P->acc0);
    mu1 = gsl_spline_eval(P->spline1, bounds(t,0.5,P->maxage), P->acc1);
    mu2 = gsl_spline_eval(P->spline2, bounds(t,0.5,P->maxage), P->acc2);
    f[Never] = -(alpha01)*y[Never];
    f[Current] = alpha01*y[Never]-(alpha12+mu1-mu0)*y[Current];
    f[Former] = alpha12*y[Current]-(alpha20+mu2-mu0)*y[Former];
    f[Reclassified] = alpha20*y[Former];
    return GSL_SUCCESS;
  }

  // Model B: observed mortality
  double
  mu_ijReclassifiedOLD (int i, int j, double t, gsl_odeiv2_driver * d)
  {
    Param* P = static_cast<Param*>(d->sys->params);
    if (i==Never && j==Current)  return exp(P->int01 + P->s01->calc(t, P->beta01));
    if (i==Current && j==Former) return exp(P->int12 + P->s12->calc(t, P->beta12));
    if (i==Former && j==Reclassified)  return exp(P->beta20);
    if (j==Death) {
      if (i==Never) return gsl_spline_eval(P->spline0, bounds(t,0.5,P->maxage), P->acc0);
      if (i==Current) return gsl_spline_eval(P->spline1, bounds(t,0.5,P->maxage), P->acc1);
      if (i==Former) return gsl_spline_eval(P->spline2, bounds(t,0.5,P->maxage), P->acc2);
      if (i==Reclassified) return gsl_spline_eval(P->spline0, bounds(t,0.5,P->maxage), P->acc0);
    }
    REprintf("Unexpected combinations of states: %i to %j.\n",i,j);
    return 0.0;
  }

  // Model C: excess mortality
  double
  mu_ijReclassified (int i, int j, double t, gsl_odeiv2_driver * d)
  {
    Param* P = static_cast<Param*>(d->sys->params);
    if (i==Never && j==Current)  return exp(P->int01 + P->s01->calc(t, P->beta01));
    if (i==Current && j==Former) return exp(P->int12 + P->s12->calc(t, P->beta12));
    if (i==Former && j==Reclassified)  return exp(P->beta20);
    if (j==Death) {
      if (i==Never) return 0.0;
      if (i==Current) return gsl_spline_eval(P->spline1, bounds(t,0.5,P->maxage), P->acc1) -
		       gsl_spline_eval(P->spline0, bounds(t,0.5,P->maxage), P->acc0);
      if (i==Former) return gsl_spline_eval(P->spline2, bounds(t,0.5,P->maxage), P->acc2) - 
		       gsl_spline_eval(P->spline0, bounds(t,0.5,P->maxage), P->acc0);
      if (i==Reclassified) return 0.0;
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
    for (size_t k=0; k<dimension; ++k)
      y[k] = 0.0;
    y[i] = 1.0;
    double tm = s;
    int status = gsl_odeiv2_driver_apply (d, &tm, t, y);
    if (status != GSL_SUCCESS)
      REprintf ("error, return value=%d\n", status);
    // cache the results
    for (size_t k=0; k<dimension; ++k)
      P_ijs[P_ij_key(i,k,s,t)] = y[k];
    return y[j];
  }

  double
  P_iK(int i, double s, double t, gsl_odeiv2_driver * d) {
    size_t dimension = d->sys->dimension;
    double total = 0.0;
    for (size_t k=0; k<dimension; ++k)
      total += P_ij(i,k,s,t,d);
    return total;
  }

  // the required mu_ij are not affected by the additional reclassification state
  double negllReclassified(int state, double s, double t, double u, gsl_odeiv2_driver * d, double freq = 1.0, int recall = 1) {
    double ll = 0.0;
    if (recall == Recall) { // recall of smoking history available
      if (state == Never) // Never->Never _and_ Never->Reclassified (equiv to current status)
	ll = log(P_ij(Never,Never,0.0,u,d)+
		 P_ij(Never,Reclassified,0.0,u,d)); // ignore s and t
      if (state == Current)
	ll = log(P_ij(Never,Never,0.0,s,d))+
	  log(mu_ijReclassified(Never,Current,s,d))+
	  log(P_ij(Current,Current,s,u,d)); // ignore t
      if (state == Former)
	ll = log(P_ij(Never,Never,0.0,s,d))+
	  log(mu_ijReclassified(Never,Current,s,d))+
	  log(P_ij(Current,Current,s,t,d))+
	  log(mu_ijReclassified(Current,Former,t,d))+
	  log(P_ij(Former,Former,t,u,d));
    }
    if (recall == CurrentStatus) {// current status only
      ll = log(P_ij(Never,state,0.0,u,d));
    }
    if (recall == FormerWithCessation && state == Former) {// recall of age quit (not initiation) for former smokers
      ll = log(P_ij(Never,Current,0.0,t,d)) + 
	log(mu_ijReclassified(Current,Former,t,d))+
	log(P_ij(Former,Former,t,u,d)); // ignores s
    }

    if (ll == 0.0) REprintf("ll==0.0? (state=%i, recall=%i)\n",state,recall);
    ll -= log(P_iK(Never,0.0,u,d));
  
    return -ll*freq;
  }


  RcppExport SEXP 
  gsl_main2Reclassified(SEXP _finalState,
			SEXP _recall,
			SEXP _time1, SEXP _time2, SEXP _time3,
			SEXP _freq,
			SEXP _int01, SEXP _int12,
			SEXP _knots01, SEXP _knots12, 
			SEXP _beta01, SEXP _beta12,
			SEXP _beta20,
			SEXP _ages0,
			SEXP _mu0,
			SEXP _mu1,
			SEXP _mu2,
			SEXP _debug)
  {
    
    bool debug = Rcpp::as<bool>(_debug);

    RcppGSL::vector<int> recall = _recall; 

    RcppGSL::vector<double> 
      finalState = _finalState,
      time1=_time1,
      time2=_time2,
      time3=_time3,
      freq = _freq,
      knots01=_knots01,
      knots12=_knots12,
      beta01 = _beta01,
      beta12 = _beta12,
      ages0 = _ages0,
      mu0 = _mu0,
      mu1 = _mu1,
      mu2 = _mu2
      ;

    double int01=as<double>(_int01), 
      int12=as<double>(_int12),
      beta20=as<double>(_beta20);
  
    gsl_spline *spline0, *spline1, *spline2;
    gsl_interp_accel *acc0, *acc1, *acc2;
    acc0 = gsl_interp_accel_alloc ();
    acc1 = gsl_interp_accel_alloc ();
    acc2 = gsl_interp_accel_alloc ();
    spline0 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    spline1 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    spline2 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    gsl_spline_init (spline0, ages0->data, mu0->data, ages0.size());
    gsl_spline_init (spline1, ages0->data, mu1->data, ages0.size());
    gsl_spline_init (spline2, ages0->data, mu2->data, ages0.size());

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

    double maxage = ages0[ages0.size()-1]+0.5;
    if (debug)
      Rprintf("maxage=%f\n",maxage);

    Param beta = {int01, int12, beta20, maxage, spline0, spline1, spline2, acc0, acc1, acc2, ns01, ns12, beta01, beta12};

    gsl_odeiv2_system sys = {funcReclassified, NULL, 4, &beta};
    gsl_odeiv2_driver * d = 
      gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
				     1e-12, 1e-12, 0.0);

    // clear the values whenever the parameters change
    P_ijs.clear();
    // without storing intermediate results, 10000 ODEs took ~20 sec

    double sum_negll = 0.0;
    for (size_t i = 0; i< finalState.size(); ++i) 
      sum_negll += negllReclassified(finalState[i], time1[i], time2[i], time3[i], d, freq[i], recall[i]);

    if (debug) {
      Rprintf("Size of Pij map = %i\n",P_ijs.size());
      Rprintf("ll(Never,0.0,0.0,70.0)=%f\n",-negllReclassified(Never,0.0,0.0,70.0,d));
      Rprintf("ll(Current,20.0,0.0,70.0)=%f\n",-negllReclassified(Current,20.0,0.0,70.0,d));
      Rprintf("ll(Former,20.0,50.0,70.0)=%f\n",-negllReclassified(Former,20.0,50.0,70.0,d));

      Rprintf("P_ij(Never,Never,0,50)=%f\n",P_ij(Never, Never, 0.0, 50.0, d));
      Rprintf("P_ij(Never,Current,0,50)=%f\n",P_ij(Never, Current, 0.0, 50.0, d));
      Rprintf("P_ij(Never,Former,0,50)=%f\n",P_ij(Never, Former, 0.0, 50.0, d));
      Rprintf("P_ij(Never,Reclassified,0,50)=%f\n",P_ij(Never, Reclassified, 0.0, 20.0, d));
      Rprintf("P_ij(Never,Death,0,50)=%f\n",1.0 - (P_ij(Never, Never, 0.0, 50.0, d)+P_ij(Never, Current, 0.0, 50.0, d)+
						   P_ij(Never, Former, 0.0, 50.0, d)+P_ij(Never, Reclassified, 0.0, 50.0, d))); // Death not included in (live) state-space
      Rprintf("mu_ij(Never,Current,50)=%f\n",mu_ijReclassified(Never, Current, 50.0, d));
      Rprintf("mu_ij(Never,Death,50)=%f\n",mu_ijReclassified(Never, Death, 50.0, d));
      Rprintf("mu_ij(Current,Former,50)=%f\n",mu_ijReclassified(Current, Former, 50.0, d));
      Rprintf("mu_ij(Current,Death,50)=%f\n",mu_ijReclassified(Current, Death, 50.0, d));
      Rprintf("mu_ij(Former,Reclassified,50)=%f\n",mu_ijReclassified(Former, Reclassified, 50.0, d));
      Rprintf("mu_ij(Former,Death,50)=%f\n",mu_ijReclassified(Former, Death, 50.0, d));
      Rprintf("mu_ij(Reclassified,Death,50)=%f\n",mu_ijReclassified(Reclassified, Death, 50.0, d));
      Rprintf("P_iK(Never,0.0,70.0)=%f\n",P_iK(Never,0.0,70.0,d));
    
      Rprintf("negll(Current,20,70)=%f\n",negllReclassified(Current, 20.0, 0.0, 70.0,d));
      Rprintf("negll(Former,20,50,70)=%f\n",negllReclassified(Former, 20.0, 50.0, 70.0,d));

    }

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
    mu1.free();
    mu2.free();
  
    gsl_spline_free (spline0);
    gsl_spline_free (spline1);
    gsl_spline_free (spline2);
    gsl_interp_accel_free (acc0);
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);
    gsl_odeiv2_driver_free (d);

    return wrap(sum_negll);
  }


  RcppExport SEXP 
  gsl_predReclassified( SEXP _maxage,
			SEXP _finalState,
			SEXP _time,
			SEXP _int01, SEXP _int12,
			SEXP _knots01, SEXP _knots12, 
			SEXP _beta01, SEXP _beta12,
			SEXP _beta20,
			SEXP _ages0,
			SEXP _mu0,
			SEXP _mu1,
			SEXP _mu2)
  {
    
    RcppGSL::vector<double> 
      finalState = _finalState,
      time=_time,
      knots01=_knots01,
      knots12=_knots12,
      beta01 = _beta01,
      beta12 = _beta12,
      ages0 = _ages0,
      mu0 = _mu0,
      mu1 = _mu1,
      mu2 = _mu2
      ;

    Rcpp::NumericVector Prob(finalState.size());

    double int01=as<double>(_int01), 
      int12=as<double>(_int12),
      beta20=as<double>(_beta20),
      maxage=as<double>(_maxage);
  
    gsl_spline *spline0, *spline1, *spline2;
    gsl_interp_accel *acc0, *acc1, *acc2;
    acc0 = gsl_interp_accel_alloc ();
    acc1 = gsl_interp_accel_alloc ();
    acc2 = gsl_interp_accel_alloc ();
    spline0 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    spline1 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    spline2 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    gsl_spline_init (spline0, ages0->data, mu0->data, ages0.size());
    gsl_spline_init (spline1, ages0->data, mu1->data, ages0.size());
    gsl_spline_init (spline2, ages0->data, mu2->data, ages0.size());

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

    Param beta = {int01, int12, beta20, maxage, spline0, spline1, spline2, acc0, acc1, acc2, ns01, ns12, beta01, beta12};

    gsl_odeiv2_system sys = {funcReclassified, NULL, 4, &beta};
    gsl_odeiv2_driver * d = 
      gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
				     1e-12, 1e-12, 0.0);

    // clear the values whenever the parameters change
    P_ijs.clear();
    // without storing intermediate results, 10000 ODEs took ~20 sec

    for (size_t i = 0; i< finalState.size(); ++i) 
      Prob[i] = P_ij(Never,finalState[i],0.0,time[i],d);

    delete ns01;
    delete ns12;

    finalState.free();
    time.free();
    knots01.free();
    knots12.free();
    beta01.free();
    beta12.free();
    ages0.free();
    mu0.free();
    mu1.free();
    mu2.free();
  
    gsl_spline_free (spline0);
    gsl_spline_free (spline1);
    gsl_spline_free (spline2);
    gsl_interp_accel_free (acc0);
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);
    gsl_odeiv2_driver_free (d);

    return wrap(Prob);
  }


  RcppExport SEXP 
  gsl_main2ReclassifiedPS(SEXP _finalState,
			  SEXP _recall,
			  SEXP _time1, SEXP _time2, SEXP _time3,
			  SEXP _freq,
			  SEXP _int01, SEXP _int12,
			  SEXP _lower01, SEXP _upper01, SEXP _nterm01, SEXP _pmatrix01, SEXP _sp01,
			  SEXP _lower12, SEXP _upper12, SEXP _nterm12, SEXP _pmatrix12, SEXP _sp12,
			  SEXP _beta01, SEXP _beta12,
			  SEXP _beta20,
			  SEXP _ages0,
			  SEXP _mu0,
			  SEXP _mu1,
			  SEXP _mu2,
			  SEXP _debug)
  {
    
    bool debug = Rcpp::as<bool>(_debug);

    RcppGSL::vector<int> recall = _recall, finalState = _finalState; 

    RcppGSL::vector<double> 
      // finalState = _finalState,
      time1=_time1,
      time2=_time2,
      time3=_time3,
      freq = _freq,
      beta01 = _beta01,
      beta12 = _beta12,
      ages0 = _ages0,
      mu0 = _mu0,
      mu1 = _mu1,
      mu2 = _mu2
      ;

    RcppGSL::matrix<double> 
      pmatrix01 = _pmatrix01, 
      pmatrix12 = _pmatrix12;

    double int01=as<double>(_int01), 
      int12=as<double>(_int12),
      beta20=as<double>(_beta20),
      lower01=as<double>(_lower01),
      upper01=as<double>(_upper01),
      sp01=as<double>(_sp01),
      lower12=as<double>(_lower12),
      upper12=as<double>(_upper12),
      sp12=as<double>(_sp12);

    int 
      nterm01=as<int>(_nterm01),
      nterm12=as<int>(_nterm12);
  
    gsl_spline *spline0, *spline1, *spline2;
    gsl_interp_accel *acc0, *acc1, *acc2;
    acc0 = gsl_interp_accel_alloc ();
    acc1 = gsl_interp_accel_alloc ();
    acc2 = gsl_interp_accel_alloc ();
    spline0 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    spline1 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    spline2 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    gsl_spline_init (spline0, ages0->data, mu0->data, ages0.size());
    gsl_spline_init (spline1, ages0->data, mu1->data, ages0.size());
    gsl_spline_init (spline2, ages0->data, mu2->data, ages0.size());

    // initial parameters for initiation
    // ps(double boundaryL, double boundaryU, bool intercept = false, bool centred = false, double centre = 0.0, size_t nterm = 10, size_t k = 4) : 

    ps * ps01 = new ps(lower01, // lower boundary
		       upper01, // upper boundary
		       false, // intercept
		       true, // centred
		       20.0, // centre
		       nterm01 // nterm
		       );

    // knots for cessation
    ps * ps12 = new ps(lower12,
		       upper12,
		       false, // intercept
		       true, // centred
		       40.0, // centre
		       nterm12
		       );

    if (debug) {
      gsl_vector_Rprintf(ps01->knots);
      Rprintf("\n");
      gsl_vector_Rprintf(ps12->knots);
      Rprintf("\n");
    }

    double maxage = ages0[ages0.size()-1]+0.5;
    if (debug)
      Rprintf("maxage=%f\n",maxage);

    Param beta = {int01, int12, beta20, maxage, spline0, spline1, spline2, acc0, acc1, acc2, ps01, ps12, beta01, beta12};

    gsl_odeiv2_system sys = {funcReclassified, NULL, 4, &beta};
    gsl_odeiv2_driver * d = 
      //gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
      gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
				     1e-12, 1e-12, 0.0);

    // clear the values whenever the parameters change
    P_ijs.clear();
    // without storing intermediate results, 10000 ODEs took ~20 sec

    double negll = 0.0;
    for (size_t i = 0; i< finalState.size(); ++i) {
      double negllComponent = negllReclassified(finalState[i], time1[i], time2[i], time3[i], d, freq[i], recall[i]); // add negative components
      if (debug) Rprintf("ll=%g; ",-negllComponent);
      negll += negllComponent;
    }

    gsl_vector *v01;
    double penalty01=0.0;
    v01 = gsl_vector_alloc(beta01->size);
    gsl_blas_dgemv(CblasNoTrans, 1.0, pmatrix01, beta01, 0.0, v01);
    gsl_blas_ddot(beta01,v01,&penalty01);
    gsl_vector_free(v01);
    double pnegll = negll + penalty01 * sp01/2.0; // add positive penalty

    gsl_vector *v12;
    double penalty12=0.0;
    v12 = gsl_vector_alloc(beta12->size);
    gsl_blas_dgemv(CblasNoTrans, 1.0, pmatrix12, beta12, 0.0, v12);
    gsl_blas_ddot(beta12,v12,&penalty12);
    gsl_vector_free(v12);
    pnegll += penalty12 * sp12/2.0; // add positive penalty

    if (debug) {
      Rprintf("Size of Pij map = %i\n",P_ijs.size());
      Rprintf("ll(Never,0.0,0.0,70.0)=%f\n",-negllReclassified(Never,0.0,0.0,70.0,d));
      Rprintf("ll(Never,0.0,0.0,70.0)=%f\n",-negllReclassified(Never,0.0,0.0,70.0,d));
      Rprintf("ll(Current,20.0,0.0,70.0)=%f\n",-negllReclassified(Current,20.0,0.0,70.0,d));
      Rprintf("ll(Former,20.0,50.0,70.0)=%f\n",-negllReclassified(Former,20.0,50.0,70.0,d));

      Rprintf("P_ij(Never,Never,0,50)=%f\n",P_ij(Current, Current, 0.0, 20.0, d));
      Rprintf("mu_ij(Never,Current,20)=%f\n",mu_ijReclassified(Never, Current, 20.0, d));
      Rprintf("P_ij(Current,Current,20,70)=%f\n",P_ij(Current, Current, 20.0, 70.0, d));
      Rprintf("P_iK(Never,0.0,70.0)=%f\n",P_iK(Never,0.0,70.0,d));
    
      Rprintf("negll(Current,20,70)=%f\n",negllReclassified(Current, 20.0, 0.0, 70.0,d));
      Rprintf("negll(Former,20,50,70)=%f\n",negllReclassified(Former, 20.0, 50.0, 70.0,d));

    }

    delete ps01;
    delete ps12;

    finalState.free();
    time1.free();
    time2.free();
    time3.free();
    beta01.free();
    beta12.free();
    ages0.free();
    mu0.free();
    mu1.free();
    mu2.free();
    pmatrix01.free();
    pmatrix12.free();
  
    gsl_spline_free (spline0);
    gsl_spline_free (spline1);
    gsl_spline_free (spline2);
    gsl_interp_accel_free (acc0);
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);
    gsl_odeiv2_driver_free (d);

    return Rcpp::List::create(Rcpp::Named("pnegll")=wrap(pnegll),
			      Rcpp::Named("negll")=wrap(negll));
  }

  
  /** @brief 
      Re-implementation using matrix inverses to calculate P(s,t). 
      This may not be numerically stable.
   **/

  std::map<P_ij_key,double> P_ijs2;

  double
  P_ij2(int i, int j, double s, double t, gsl_odeiv2_driver * d, Rcpp::NumericMatrix P, Rcpp::NumericMatrix invP) {
    // do we have the values cached?
    // P_ij_key key = P_ij_key(i,j,s,t);
    // if (P_ijs2.count(key)>0) return P_ijs2[key];
    // otherwise...
    int fulldim = d->sys->dimension+1; // assumes death state is the last state
    double pij = 0.0;
    if (s==0) { 
      pij = P(i+j*fulldim,t);
    } else {
      for (int k=0; k<fulldim-1; ++k) {
	pij += invP(k+i*fulldim,s+1) * P(j+k*fulldim,t+1); // offset=1
      }
    }
    // cache the results
    // P_ijs2[key] = pij;
    return pij;
  }

  double
  P_iK2(int i, double s, double t, gsl_odeiv2_driver * d, Rcpp::NumericMatrix P, Rcpp::NumericMatrix invP) {
    double total = 0.0;
    int fulldim = d->sys->dimension+1;
    for (int k=0; k<fulldim-1; ++k)
      total += P_ij2(i,k,s,t,d,P,invP);
    return total;
  }


  // Calculating the probability matrices cf. the log-likelihood
  RcppExport SEXP 
  gsl_main2ReclassifiedPS2(SEXP _finalState,
			  SEXP _recall,
			  SEXP _time1, SEXP _time2, SEXP _time3,
			  SEXP _freq,
			  SEXP _int01, SEXP _int12,
			  SEXP _lower01, SEXP _upper01, SEXP _nterm01, SEXP _pmatrix01, SEXP _sp01,
			  SEXP _lower12, SEXP _upper12, SEXP _nterm12, SEXP _pmatrix12, SEXP _sp12,
			  // SEXP _lower01, SEXP _upper01, SEXP _nterm01, 
			  // SEXP _lower12, SEXP _upper12, SEXP _nterm12, 
			  SEXP _beta01, SEXP _beta12,
			  SEXP _beta20,
			  SEXP _ages0,
			  SEXP _mu0,
			  SEXP _mu1,
			  SEXP _mu2,
			  SEXP _debug)
  {
    
    bool debug = Rcpp::as<bool>(_debug);

    RcppGSL::vector<int> recall = _recall, finalState = _finalState; 

    RcppGSL::vector<double> 
      // finalState = _finalState,
      time1=_time1,
      time2=_time2,
      time3=_time3,
      freq = _freq,
      beta01 = _beta01,
      beta12 = _beta12,
      ages0 = _ages0,
      mu0 = _mu0,
      mu1 = _mu1,
      mu2 = _mu2
      ;

    RcppGSL::matrix<double> 
      pmatrix01 = _pmatrix01, 
      pmatrix12 = _pmatrix12;

    double int01=as<double>(_int01), 
      int12=as<double>(_int12),
      beta20=as<double>(_beta20),
      lower01=as<double>(_lower01),
      upper01=as<double>(_upper01),
      sp01=as<double>(_sp01),
      lower12=as<double>(_lower12),
      upper12=as<double>(_upper12),
      sp12=as<double>(_sp12);

    int 
      nterm01=as<int>(_nterm01),
      nterm12=as<int>(_nterm12);
  
    gsl_spline *spline0, *spline1, *spline2;
    gsl_interp_accel *acc0, *acc1, *acc2;
    acc0 = gsl_interp_accel_alloc ();
    acc1 = gsl_interp_accel_alloc ();
    acc2 = gsl_interp_accel_alloc ();
    spline0 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    spline1 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    spline2 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    gsl_spline_init (spline0, ages0->data, mu0->data, ages0.size());
    gsl_spline_init (spline1, ages0->data, mu1->data, ages0.size());
    gsl_spline_init (spline2, ages0->data, mu2->data, ages0.size());

    // initial parameters for initiation
    // ps(double boundaryL, double boundaryU, bool intercept = false, bool centred = false, double centre = 0.0, size_t nterm = 10, size_t k = 4) : 

    ps * ps01 = new ps(lower01, // lower boundary
		       upper01, // upper boundary
		       false, // intercept
		       true, // centred
		       20.0, // centre
		       nterm01 // nterm
		       );

    // knots for cessation
    ps * ps12 = new ps(lower12,
		       upper12,
		       false, // intercept
		       true, // centred
		       40.0, // centre
		       nterm12
		       );

    size_t K = ages0.size()-1;
    double maxage = ages0[K]+0.5;

    // clear the values whenever the parameters change
    // P_ijs2.clear();

    Param beta = {int01, int12, beta20, maxage, spline0, spline1, spline2, acc0, acc1, acc2, ps01, ps12, beta01, beta12};

    gsl_odeiv2_system sys = {funcReclassified, NULL, 4, &beta};
    gsl_odeiv2_driver * d = 
      gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
				     1e-12, 1e-12, 0.0);

    size_t dimension = d->sys->dimension;
    size_t fulldim = dimension+1;
    // 3-D matrix
    Rcpp::NumericMatrix P(fulldim*fulldim,K), invP(fulldim*fulldim,K);
    double y[dimension];
    // RcppGSL::matrix<double> P(fulldim,fulldim), inv(fulldim,fulldim);
    gsl_matrix * Pmat = gsl_matrix_alloc(fulldim,fulldim);
    gsl_matrix * invPmat = gsl_matrix_alloc(fulldim,fulldim);
    gsl_permutation * perm = gsl_permutation_alloc (fulldim);
    // loop over the dimensions
    for (size_t i=0; i<dimension; ++i) {
      // initial values
      for (size_t j=0; j<dimension; ++j)
	y[j] = 0.0;
      y[i] = 1.0;
      double s, t;
      size_t k;
      for (s = 0.0, t = 0.5, k = 0; k<K; s=t, ++t, ++k) { // loop over end times
	int status = gsl_odeiv2_driver_apply (d, &s, t, y);
	if (status != GSL_SUCCESS)
	  REprintf ("error, return value=%d\n", status);
	for (size_t j=0; j<dimension; ++j) {
	  // if (y[j]<1e-7 && y[j]>-1e-7) 
	  //   P(i+fulldim*j, k) = 0; // zapsmall
	  //else 
	  P(i+fulldim*j, k) = y[j];
	}
      }
    }
    // add in column sums
    for (size_t k = 0; k<=K; ++k) {
      for (size_t i=0; i<fulldim; ++i) {
	P(i+fulldim*(fulldim-1),k) = 1;
	for (size_t j=0; j<fulldim-1; ++j) {
	  P(i+fulldim*(fulldim-1),k) -= P(i+j*fulldim,k);
	}
      }
    }
    // calculate inverses
    for (size_t k = 0; k<K; ++k) {
      for (size_t i=0; i<fulldim; ++i) {
	for (size_t j=0; j<fulldim; ++j) {
	  gsl_matrix_set(Pmat,i,j,P(i+j*fulldim,k));
	}
      }
      // if (debug) {
      // 	Rprintf("k=%i\n",k);
      // 	gsl_matrix_Rprintf(Pmat);
      // }
      int signum;
      gsl_linalg_LU_decomp (Pmat, perm, &signum);
      gsl_linalg_LU_invert (Pmat, perm, invPmat);
      for (size_t i=0; i<fulldim; ++i) {
	for (size_t j=0; j<fulldim; ++j) {
	  invP(i+j*fulldim,k) = gsl_matrix_get(invPmat,i,j);
	}
      }
    }

    double negll = 0.0;
    for (size_t i = 0; i< finalState.size(); ++i) {
      double ll = 0.0;
      int state = finalState[i];
      double s=time1[i], t=time2[i], u=time3[i];
      // double negllReclassified(int state, double s, double t, double u, gsl_odeiv2_driver * d, double freq = 1.0, int recall = 1)
      // negll += negllReclassified(finalState[i], time1[i], time2[i], time3[i], d, freq[i], recall[i]); // add negative components
      if (recall[i] == Recall) { // recall of smoking history available
	if (state == Never) // Never->Never _and_ Never->Reclassified
	  ll = log(P_ij2(Never,Never,0.0,u,d,P,invP)+
		   P_ij2(Never,Reclassified,0.0,u,d,P,invP)); // ignore s and t
	if (state == Current)
	  ll = log(P_ij2(Never,Never,0.0,s,d,P,invP))+
	    log(mu_ijReclassified(Never,Current,s,d))+
	    log(P_ij2(Current,Current,s,u,d,P,invP)); // ignore t
	if (state == Former)
	  ll = log(P_ij2(Never,Never,0.0,s,d,P,invP))+
	    log(mu_ijReclassified(Never,Current,s,d))+
	    log(P_ij2(Current,Current,s,t,d,P,invP))+
	    log(mu_ijReclassified(Current,Former,t,d))+
	    log(P_ij2(Former,Former,t,u,d,P,invP));
      }
      if (recall[i] == CurrentStatus) {// current status only
	ll = log(P_ij2(Never,state,0.0,u,d,P,invP));
      }
      if (recall[i] == FormerWithCessation && state == Former) {// recall of age quit (not initiation) for former smokers
	ll = log(P_ij2(Never,Current,0.0,t,d,P,invP)) + 
	  log(mu_ijReclassified(Current,Former,t,d))+
	  log(P_ij2(Former,Former,t,u,d,P,invP)); // ignores s
      }
      if (ll == 0.0) REprintf("ll==0.0? (state=%i, recall=%i)\n",state,recall[i]);
      ll -= log(P_iK2(Never,0.0,u,d,P,invP));
      ll *= freq[i];
      if (debug) Rprintf("ll=%g; ",ll);
      negll += -ll;
    }

    gsl_vector *v01;
    double penalty01=0.0;
    v01 = gsl_vector_alloc(beta01->size);
    gsl_blas_dgemv(CblasNoTrans, 1.0, pmatrix01, beta01, 0.0, v01);
    gsl_blas_ddot(beta01,v01,&penalty01);
    gsl_vector_free(v01);
    double pnegll = negll + penalty01 * sp01/2.0; // add positive penalty

    gsl_vector *v12;
    double penalty12=0.0;
    v12 = gsl_vector_alloc(beta12->size);
    gsl_blas_dgemv(CblasNoTrans, 1.0, pmatrix12, beta12, 0.0, v12);
    gsl_blas_ddot(beta12,v12,&penalty12);
    gsl_vector_free(v12);
    pnegll += penalty12 * sp12/2.0; // add positive penalty


    delete ps01;
    delete ps12;

    beta01.free();
    beta12.free();
    ages0.free();
    mu0.free();
    mu1.free();
    mu2.free();
  
    gsl_matrix_free (Pmat);
    gsl_matrix_free (invPmat);
    gsl_spline_free (spline0);
    gsl_spline_free (spline1);
    gsl_spline_free (spline2);
    gsl_interp_accel_free (acc0);
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);
    gsl_odeiv2_driver_free (d);
    gsl_permutation_free (perm);

    return Rcpp::List::create(Rcpp::Named("P")=wrap(P),
			      Rcpp::Named("invP")=wrap(invP),
			      Rcpp::Named("negll")=wrap(negll),
			      Rcpp::Named("pnegll")=wrap(pnegll));
  }


  /** @brief
      Implementation that also calculates partials for P_ij(s,t)
   **/

  // define the ODE (Model C, using excess mortality)
  int
  funcReclassifiedExtended (double t, 
			    const double y[],   // dimension: 4*(ntheta+1)
			    double f[],         // dimension: 4*(ntheta+1)
			    void *params)       // type: Param*
  {
    double alpha01, alpha12, alpha20, mu0, mu1, mu2;
    double dalpha01, dalpha12, dalpha20;
    Param* P = static_cast<Param*>(params);
    int ntheta = P->beta01->size + P->beta12->size + 3;
    gsl_vector *x01, *x12, *x20;
    gsl_vector *_beta01, *_beta12;
    int j = 0;

    // allocate vectors
    x01 = gsl_vector_alloc(ntheta);
    x12 = gsl_vector_alloc(ntheta);
    x20 = gsl_vector_alloc(ntheta);
    _beta01 = gsl_vector_alloc(P->beta01->size);
    _beta12 = gsl_vector_alloc(P->beta12->size);
    
    // initialise vectors to zero
    gsl_vector_set_all(x01,0.0);
    gsl_vector_set_all(x12,0.0);
    gsl_vector_set_all(x20,0.0);

    // calculate design vectors for the splines
    P->s01->calcBasis(t, _beta01, false);
    P->s12->calcBasis(t, _beta12, false);

    // set theta and x 
    // intercept
    gsl_vector_set(x01,j,1);
    // spline parameters
    for (size_t i = 0; i<P->beta01->size; ++i) {
      gsl_vector_set(x01,j,gsl_vector_get(_beta01,i));
    }
    // intercept
    gsl_vector_set(x12,j,1);
    // spline parameters
    for (size_t i = 0; i<P->beta12->size; ++i) {
      gsl_vector_set(x12, j, gsl_vector_get(_beta12, i));
    }
    gsl_vector_set(x20,j,1);

    alpha01 = exp(P->int01 + P->s01->calc(t, P->beta01));
    alpha12 = exp(P->int12 + P->s12->calc(t, P->beta12));
    alpha20 = exp(P->beta20);
    mu0 = gsl_spline_eval(P->spline0, bounds(t,0.5,P->maxage), P->acc0);
    mu1 = gsl_spline_eval(P->spline1, bounds(t,0.5,P->maxage), P->acc1);
    mu2 = gsl_spline_eval(P->spline2, bounds(t,0.5,P->maxage), P->acc2);

    // base equations
    f[Never] = -alpha01*y[Never];
    f[Current] = alpha01*y[Never]-(alpha12+mu1-mu0)*y[Current];
    f[Former] = alpha12*y[Current]-(alpha20+mu2-mu0)*y[Former];
    f[Reclassified] = alpha20*y[Former];

    // d/dt (partial P(s,t)/partial theta) = (partial P(s,t)/partial theta)*alpha + P(s,t)(partial alpha / partial theta)
    for (int offset=4, i=0; i<ntheta; ++i, offset += 4) {
      dalpha01 = alpha01 * gsl_vector_get(x01,i);
      dalpha12 = alpha12 * gsl_vector_get(x12,i);
      dalpha20 = alpha20 * gsl_vector_get(x20,i);
      f[Never+offset] = - (alpha01*y[Never+offset] + dalpha01*y[Never]);
      f[Current+offset] = alpha01*y[Never+offset] + dalpha01*y[Never] -
	((alpha12+mu1-mu0)*y[Current+offset] + dalpha12*y[Current]);
      f[Former+offset] = (alpha12+mu2-mu0)*y[Current+offset]+dalpha12*y[Current] -
	((alpha20+mu2-mu0)*y[Former+offset] + dalpha20*y[Former]);
      f[Reclassified+offset] = alpha20*y[Former+offset] + dalpha20*y[Former];
    }

    // tidy up
    gsl_vector_free(x01);
    gsl_vector_free(x12);
    gsl_vector_free(x20);
    gsl_vector_free(_beta01);
    gsl_vector_free(_beta12);

    return GSL_SUCCESS;
  }

  // for efficiency, cache the values for P_ij
  struct P_ij_report {
    std::vector<double> P; 
    std::vector<double> dPdtheta;
  };
  std::map<P_ij_key,P_ij_report> P_ij_reports;

  P_ij_report
  p_ij_report(int i, int j, double s, double t, gsl_odeiv2_driver * d) {
    // do we have the values cached?
    P_ij_key key = P_ij_key(i,j,s,t);
    if (P_ij_reports.count(key)>0) return P_ij_reports[key];
    // otherwise...
    size_t dimension = d->sys->dimension;
    size_t n_equations = 4;
    size_t M = dimension / n_equations; // number of parameters + 1 
    double y[dimension];
    for (size_t k=0; k<dimension; ++k)
      y[k] = 0.0;
    y[i] = 1.0;
    double tm = s;
    int status = gsl_odeiv2_driver_apply (d, &tm, t, y);
    if (status != GSL_SUCCESS)
      REprintf ("error, return value=%d\n", status);
    // cache the results
    for (size_t k=0; k<4; ++k) { // state
      P_ij_report out;
      for (size_t m=0; m<n_equations; ++m) 
	out.P.push_back(y[m]);
      for (size_t m=n_equations; m<M; ++m) 
	out.dPdtheta.push_back(y[m+k*M]);
      P_ij_reports[P_ij_key(i,k,s,t)] = out;
    }
    return P_ij_reports[P_ij_key(i,j,s,t)];
  }

  RcppExport SEXP 
  gsl_predReclassifiedExtended( SEXP _maxage,
			SEXP _finalState,
			SEXP _time,
			SEXP _int01, SEXP _int12,
			SEXP _knots01, SEXP _knots12, 
			SEXP _beta01, SEXP _beta12,
			SEXP _beta20,
			SEXP _ages0,
			SEXP _mu0,
			SEXP _mu1,
			SEXP _mu2)
  {
    
    RcppGSL::vector<double> 
      finalState = _finalState,
      time=_time,
      knots01=_knots01,
      knots12=_knots12,
      beta01 = _beta01,
      beta12 = _beta12,
      ages0 = _ages0,
      mu0 = _mu0,
      mu1 = _mu1,
      mu2 = _mu2
      ;

    // FIX OUTPUT
    std::vector<P_ij_report> reports;

    double int01=as<double>(_int01), 
      int12=as<double>(_int12),
      beta20=as<double>(_beta20),
      maxage=as<double>(_maxage);
  
    gsl_spline *spline0, *spline1, *spline2;
    gsl_interp_accel *acc0, *acc1, *acc2;
    acc0 = gsl_interp_accel_alloc ();
    acc1 = gsl_interp_accel_alloc ();
    acc2 = gsl_interp_accel_alloc ();
    spline0 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    spline1 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    spline2 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
    gsl_spline_init (spline0, ages0->data, mu0->data, ages0.size());
    gsl_spline_init (spline1, ages0->data, mu1->data, ages0.size());
    gsl_spline_init (spline2, ages0->data, mu2->data, ages0.size());

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

    Param beta = {int01, int12, beta20, maxage, spline0, spline1, spline2, acc0, acc1, acc2, ns01, ns12, beta01, beta12};

    size_t n_coef=3+3+3; // Is this correct?
    size_t n_equations = 4;
    size_t M = (n_coef+1) * n_equations;

    gsl_odeiv2_system sys = {funcReclassifiedExtended, NULL, M, &beta};
    gsl_odeiv2_driver * d = 
      gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
				     1e-10, 1e-10, 0.0);

    // clear the values whenever the parameters change
    P_ij_reports.clear();

    for (size_t i = 0; i< finalState.size(); ++i) {
      reports.push_back(p_ij_report(Never,finalState[i],0.0,time[i],d)); 
    }

    delete ns01;
    delete ns12;

    finalState.free();
    time.free();
    knots01.free();
    knots12.free();
    beta01.free();
    beta12.free();
    ages0.free();
    mu0.free();
    mu1.free();
    mu2.free();
  
    gsl_spline_free (spline0);
    gsl_spline_free (spline1);
    gsl_spline_free (spline2);
    gsl_interp_accel_free (acc0);
    gsl_interp_accel_free (acc1);
    gsl_interp_accel_free (acc2);
    gsl_odeiv2_driver_free (d);

    // return reports;
    return wrap(1);
  }


  // RcppExport SEXP test_gsl() {

  //   gsl_matrix *m;
  //   m = gsl_matrix_alloc(2,2);
  //   gsl_matrix_set(m,0,0,1.0);
  //   gsl_matrix_set(m,1,0,2.0);
  //   gsl_matrix_set(m,0,1,3.0);
  //   gsl_matrix_set(m,1,1,4.0);
  //   gsl_vector *b;
  //   b = gsl_vector_alloc(2);
  //   gsl_vector_set(b,0,1.0);
  //   gsl_vector_set(b,1,2.0);
  //   gsl_vector *v;
  //   double penalty=0.0;
  //   v = gsl_vector_alloc(b->size);
  //   gsl_blas_dgemv(CblasNoTrans, 1.0, m, b, 0.0, v);
  //   gsl_blas_ddot(b,v,&penalty);
  //   gsl_vector_free(v);
  //   gsl_vector_free(b);
  //   gsl_matrix_free(m);
  //   return(wrap(penalty));
  // }
  
}
