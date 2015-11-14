#include <Rcpp.h>
#include <RcppGSL.h>
//
#include <cstdlib>
#include <map>
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
  typedef RcppGSL::vector<double> Vector;
  typedef std::map<P_ij_key,Vector> dP_ij_t;
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
  void Rprint(RcppGSL::matrix<double> mat) {
    size_t i, j;
    for (i=0; i < mat->size1; ++i) {
      for (j=0; j < mat->size2; ++j) 
	Rprintf("%f ", gsl_matrix_get(mat, i, j));
      Rprintf("\n");
    }
  }
  void Rprint(Vector vec) {
    size_t i;
    for (i=0; i < vec->size; ++i) {
      Rprintf("%f ", gsl_vector_get(vec, i));
    }
  }
  // template<class T> T bounds(T x, T lo, T hi) { return x<lo ? lo : (x>hi ? hi : x); } 
  inline
  double bounds(double x, double lo, double hi) { return x<lo ? lo : (x>hi ? hi : x); } 
  
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
    size_t nstate, ncoef, nbeta01, nbeta12;
    gsl_spline *spline0, *spline1, *spline2;
    gsl_interp_accel *acc0, *acc1, *acc2;
    splineBasis *s01;
    splineBasis *s12;
    std::map<P_ij_key,double> P_ijs;
    std::map<P_ij_key,Vector> dP_ijs;
    gsl_odeiv2_system sys;
    gsl_odeiv2_driver *d;
    // data fields
    RcppGSL::vector<int> recall, finalState;
    Vector 
    time1, time2, time3, freq, beta01, beta12, ages0, mu0, mu1, mu2;
    Purged(SEXP sexp, int N, int nbeta01, int nbeta12, int nages0) : 
      nstate(4), ncoef(3+nbeta01+nbeta12), nbeta01(nbeta01), nbeta12(nbeta12),
      recall(N), finalState(N), time1(N), time2(N), time3(N), freq(N), beta01(nbeta01), 
      beta12(nbeta12), ages0(nages0), mu0(nages0), mu1(nages0), mu2(nages0) {
      List args = as<List>(sexp);
      debug = args("debug");
      finalState = args("finalState"); 
      recall = args("recall");
      time1 = args("time1"); 
      time2 = args("time2"); 
      time3 = args("time3"); 
      freq = args("freq"); 
      beta01 = args("beta01"); 
      beta12 = args("beta12"); 
      ages0 = args("ages0"); 
      mu0 = args("mu0"); 
      mu1 = args("mu1"); 
      mu2 = args("mu2"); 
      int01 = args("int01"); 
      int12 = args("int12"); 
      beta20 = args("beta20"); 
      acc0 = gsl_interp_accel_alloc ();
      acc1 = gsl_interp_accel_alloc ();
      acc2 = gsl_interp_accel_alloc ();
      spline0 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
      spline1 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
      spline2 = gsl_spline_alloc (gsl_interp_cspline, ages0.size());
      gsl_spline_init (spline0, ages0->data, mu0->data, ages0.size());
      gsl_spline_init (spline1, ages0->data, mu1->data, ages0.size());
      gsl_spline_init (spline2, ages0->data, mu2->data, ages0.size());
      maxage = ages0[ages0.size()-1]+0.5;
      if (debug)
	Rprintf("maxage=%f\n",maxage);
      sys.jacobian = 0; 
      sys.params = (void *) this;
      default_sys();
    }
    virtual ~Purged() {
      finalState.free();
      recall.free();
      time1.free();
      time2.free();
      time3.free();
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
    }
    void default_sys() {
      sys.function = funcReclassified;
      sys.dimension = nstate;
      d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
					 1e-10, 1e-10, 0.0);
    }
    void extended_sys() {
      sys.function = funcReclassifiedExtended;
      sys.dimension = nstate*(ncoef+1);
      d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
					 1e-10, 1e-10, 0.0);
    }
    double P_ij(int i, int j, double s, double t) {
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
	P_ijs[P_ij_key(i,k,s,t)] = y[k];
      return y[j];
    }
    double P_iK(int i, double s, double t) {
      double total = 0.0;
      for (size_t k=0; k<nstate; ++k)
	total += P_ij(i,k,s,t);
      return total;
    }
    Vector
    dP_ij(int i, int j, double s, double t) {
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
	Vector out(ncoef);
	P_ij_key key_ik = P_ij_key(i,k,s,t);
	P_ijs[key_ik] = y[k];
	for (size_t m=0; m<ncoef; ++m) 
	  out[m] = y[nstate*(m+1)+k]; // (y, P1', P2', ..., Pncoef')
	dP_ijs.insert(dP_ij_t::value_type(key_ik, out));
      }
      return dP_ijs.find(key)->second;
    }
    Vector dP_iK(int i, double s, double t) {
      Vector total(ncoef);
      for (size_t k=0; k<nstate; ++k) {
	gsl_vector_add(total,dP_ij(i,k,s,t));
      }
      return total;
    }
    // Excess mortality model
    double
    mu_ij (int i, int j, double t)
    {
      if (i==Never && j==Current)  return exp(int01 + s01->calc(t, beta01));
      if (i==Current && j==Former) return exp(int12 + s12->calc(t, beta12));
      if (i==Former && j==Reclassified)  return exp(beta20);
      if (j==Death) {
	if (i==Never) return 0.0;
	if (i==Current) return gsl_spline_eval(spline1, bounds(t,0.5,maxage), acc1) -
			  gsl_spline_eval(spline0, bounds(t,0.5,maxage), acc0);
	if (i==Former) return gsl_spline_eval(spline2, bounds(t,0.5,maxage), acc2) - 
			 gsl_spline_eval(spline0, bounds(t,0.5,maxage), acc0);
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
    Vector rmu_ij(int i, int j, double t)
    {
      Vector out(ncoef);
      gsl_vector_set_all(out, 0.0);
      if (i+1 != j) return out;
      Vector x01(ncoef), x12(ncoef), x20(ncoef);
      Vector _beta01(nbeta01), _beta12(nbeta12);
      int k = 0;
      // initialise vectors to zero
      gsl_vector_set_all(x01,0.0);
      gsl_vector_set_all(x12,0.0);
      gsl_vector_set_all(x20,0.0);
      // calculate design vectors for the splines
      s01->calcBasis(t, _beta01);
      s12->calcBasis(t, _beta12);
      // set theta and x 
      // intercept
      x01[k++] = 1;
      // spline parameters
      for (size_t n = 0; n<nbeta01; ++n) {
	x01[k++]=gsl_vector_get(_beta01,n);
      }
      // intercept
      x12[k++] = 1;
      // spline parameters
      for (size_t n = 0; n<nbeta12; ++n) {
	x12[k++] = gsl_vector_get(_beta12,n);
      }
      x20[k++] = 1;
      // d/dt (partial P(s,t)/partial theta) = (partial P(s,t)/partial theta)*alpha + P(s,t)(partial alpha / partial theta)
      if (i==Never && j==Current) 
	return x01;
      if (i==Current && j==Former)
	return x12;
      if (i==Former && j==Reclassified)
	return x20;
      REprintf("Unexpected combinations of states in rmu_ij: %i to %j.\n",i,j);
      return out;
    }
    double negll_component(int state, double s, double t, double u, double freq = 1.0, int recall = 1) {
      double ll = 0.0;
      if (recall == Recall) { // recall of smoking history available
	if (state == Never) // Never->Never _and_ Never->Reclassified (equiv to current status)
	  ll = log(P_ij(Never,Never,0.0,u)+
		   P_ij(Never,Reclassified,0.0,u)); // ignore s and t
	if (state == Current)
	  ll = log(P_ij(Never,Never,0.0,s))+
	    log(mu_ij(Never,Current,s))+
	    log(P_ij(Current,Current,s,u)); // ignore t
	if (state == Former)
	  ll = log(P_ij(Never,Never,0.0,s))+
	    log(mu_ij(Never,Current,s))+
	    log(P_ij(Current,Current,s,t))+
	    log(mu_ij(Current,Former,t))+
	    log(P_ij(Former,Former,t,u));
      }
      if (recall == CurrentStatus) {// current status only
	ll = log(P_ij(Never,state,0.0,u));
      }
      if (recall == FormerWithCessation && state == Former) {// recall of age quit (not initiation) for former smokers
	ll = log(P_ij(Never,Current,0.0,t)) + 
	  log(mu_ij(Current,Former,t))+
	  log(P_ij(Former,Former,t,u)); // ignores s
      }
      if (ll == 0.0) REprintf("ll==0.0? (state=%i, recall=%i)\n",state,recall);
      ll -= log(P_iK(Never,0.0,u));
      return -ll*freq;
    }
    Vector negll_gradient_component(int state, double s, double t, double u, double freq = 1.0, int recall = 1) {
      Vector gradient(ncoef);
      gsl_vector_set_all(gradient,0.0);
      if (recall == Recall) { // recall of smoking history available
	if (state == Never) // Never->Never _and_ Never->Reclassified (equiv to current status)
	  for (size_t k=0; k<ncoef; ++k)
	    gradient[k] = (dP_ij(Never,Never,0.0,u)[k]+
			   dP_ij(Never,Reclassified,0.0,u)[k]) / 
	      (P_ij(Never,Never,0.0,u)+P_ij(Never,Reclassified,0.0,u)); // ignore s and t
	if (state == Current)
	  for (size_t k=0; k<ncoef; ++k)
	    gradient[k] = dP_ij(Never,Never,0.0,s)[k] / P_ij(Never,Never,0.0,s)+
	      rmu_ij(Never,Current,s)[k]+
	      dP_ij(Current,Current,s,u)[k]/P_ij(Current,Current,s,u); // ignore t
	if (state == Former)
	  for (size_t k=0; k<ncoef; ++k)
	    gradient[k] = dP_ij(Never,Never,0.0,s)[k]/P_ij(Never,Never,0.0,s)+
	      rmu_ij(Never,Current,s)[k]+
	      dP_ij(Current,Current,s,t)[k]/P_ij(Current,Current,s,t)+
	      rmu_ij(Current,Former,t)[k]+
	      dP_ij(Former,Former,t,u)[k]/P_ij(Former,Former,t,u);
      }
      if (recall == CurrentStatus) {// current status only
	for (size_t k=0; k<ncoef; ++k)
	  gradient[k] = dP_ij(Never,state,0.0,u)[k]/P_ij(Never,state,0.0,u);
      }
      if (recall == FormerWithCessation && state == Former) {// recall of age quit (not initiation) for former smokers
	for (size_t k=0; k<ncoef; ++k)
	  gradient[k] = dP_ij(Never,Current,0.0,t)[k]/P_ij(Never,Current,0.0,t) + 
	    rmu_ij(Current,Former,t)[k]+
	    dP_ij(Former,Former,t,u)[k]/P_ij(Former,Former,t,u); // ignores s
      }
      // if (ll == 0.0) REprintf("ll==0.0? (state=%i, recall=%i)\n",state,recall);
      for (size_t k=0; k<ncoef; ++k)
	gradient[k] = gradient[k] - dP_iK(Never,0.0,u)[k] / P_iK(Never,0.0,u);
      for (size_t k=0; k<ncoef; ++k)
	gradient[k] = -gradient[k]*freq;
      return gradient;
    }
    Vector
    negll_gradient() {
      // clear the stored values
      dP_ijs.clear();
      P_ijs.clear();
      Vector total(ncoef);
      gsl_vector_set_all(total, 0.0);
      for (size_t i = 0; i< finalState.size(); ++i) 
	gsl_vector_add(total,negll_gradient_component(finalState[i], time1[i], time2[i], time3[i], freq[i], recall[i]));
      return total;
    }
    double negll() {
      // clear the stored values
      P_ijs.clear();
      double total = 0.0;
      for (size_t i = 0; i< finalState.size(); ++i) 
	total += negll_component(finalState[i], time1[i], time2[i], time3[i], freq[i], recall[i]);
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
      }
      return total;
    }
  };

  class PurgedNS : public Purged {
  public:
    PurgedNS(SEXP sexp, int N, int nbeta01, int nbeta12, int nages0) : 
      Purged(sexp, N, nbeta01, nbeta12, nages0) {
      List args = as<List>(sexp);
      // initial parameters for initiation
      s01 = new ns(as<Vector>(args("knots01")), 
		   false, // intercept
		   true, // centred
		   20.0 // centre
		   );
      // knots for cessation
      s12 = new ns(as<Vector>(args("knots12")), 
		   false, // intercept
		   true, // centred
		   40.0 // centre
		   );
    }
    virtual ~PurgedNS() {
      delete s01; 
      delete s12;
    }
  };

  class PurgedPS : public Purged {
  public:
    RcppGSL::matrix<double> pmatrix01, pmatrix12;
    double sp01, sp12;
    PurgedPS(SEXP sexp, int N, int nbeta01, int nbeta12, int nages0) : 
      Purged(sexp, N, nbeta01, nbeta12, nages0), 
      pmatrix01(nbeta01,nbeta01), pmatrix12(nbeta12, nbeta12) {
      List args = as<List>(sexp);
      pmatrix01 = args("pmatrix01");
      pmatrix12 = args("pmatrix12");
      sp01 = args("sp01");
      sp12 = args("sp12");
      s01 = new ps(as<double>(args("lower01")), // lower boundary
		   as<double>(args("upper01")), // upper boundary
		   false, // intercept
		   true, // centred
		   20.0, // centre
		   as<int>(args("nterm01")) // nterm
		   );
      // knots for cessation
      s12 = new ps(as<double>(args("lower12")),
		   as<double>(args("upper12")),
		   false, // intercept
		   true, // centred
		   40.0, // centre
		   as<int>(args("nterm12"))
		   );
    }
    double negll() {
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
    virtual ~PurgedPS() {
      delete s01; delete s12;
      pmatrix01.free(); pmatrix12.free();
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
    Vector x01(P->ncoef), x12(P->ncoef), x20(P->ncoef);
    Vector basis01(P->nbeta01), basis12(P->nbeta12);
    int j = 0;
    // initialise vectors to zero
    gsl_vector_set_all(x01,0.0);
    gsl_vector_set_all(x12,0.0);
    gsl_vector_set_all(x20,0.0);
    // calculate design vectors for the splines
    P->s01->calcBasis(t, basis01);
    P->s12->calcBasis(t, basis12);
    // set x 
    // intercept for 0->1
    x01[j++] = 1;
    // spline parameters for 0->1
    for (size_t i = 0; i < P->nbeta01; ++i) {
      x01[j++] = gsl_vector_get(basis01,i); // unclear why the gsl_vector_get() call is needed
    }
    // intercept for 1->2
    x12[j++] = 1;
    // spline parameters for 1->2
    for (size_t i = 0; i < P->nbeta12; ++i) {
      x12[j++] = gsl_vector_get(basis12,i);
    }
    // intercept for 2->0
    x20[j++] = 1;
    // calculate rates
    alpha01 = P->mu_ij(Never,Current,t);
    alpha12 = P->mu_ij(Current,Former,t);
    alpha20 = P->mu_ij(Former,Reclassified,t);
    mu0 = P->mu_ij(Never,Death,t);
    mu1 = P->mu_ij(Current,Death,t);
    mu2 = P->mu_ij(Former,Death,t);
    // base equations
    f[Never] = -alpha01*y[Never];
    f[Current] = alpha01*y[Never]-(alpha12+mu1-mu0)*y[Current];
    f[Former] = alpha12*y[Current]-(alpha20+mu2-mu0)*y[Former];
    f[Reclassified] = alpha20*y[Former];
    // d/dt (partial P(s,t)/partial theta) = (partial P(s,t)/partial theta)*alpha + P(s,t)(partial alpha / partial theta)
    for (int offset=P->nstate, i=0; i < (int) P->ncoef; ++i, offset += P->nstate) {
      dalpha01 = alpha01 * x01[i];
      dalpha12 = alpha12 * x12[i];
      dalpha20 = alpha20 * x20[i];
      f[Never+offset] = - (alpha01*y[Never+offset] + dalpha01*y[Never]);
      f[Current+offset] = alpha01*y[Never+offset] + dalpha01*y[Never] -
	((alpha12+mu1-mu0)*y[Current+offset] + dalpha12*y[Current]);
      f[Former+offset] = alpha12*y[Current+offset]+dalpha12*y[Current] -
	((alpha20+mu2-mu0)*y[Former+offset] + dalpha20*y[Former]);
      f[Reclassified+offset] = alpha20*y[Former+offset] + dalpha20*y[Former];
    }
    // tidy up - is this strictly needed?
    x01.free(); x12.free(); x20.free();
    basis01.free(); basis12.free();
    return GSL_SUCCESS;
  }

  RcppExport
  SEXP call_purged_ns(SEXP sexp) {
    List args = as<List>(sexp);
    // bool debug = as<bool>(args("debug"));
    PurgedNS model(sexp, 
		   as<int>(args("N")), 
		   as<int>(args("nbeta01")),
		   as<int>(args("nbeta12")),
		   as<int>(args("nages0")));
    std::string output_type = args("output_type");
    if (output_type == "negll")
      return wrap(model.negll());
    if (output_type == "negll_gradient")
      return wrap(model.negll_gradient());
    if (output_type == "P_ij")
      return wrap(model.P_ij(Never,Reclassified,0.0, 50.0));
    if (output_type == "dP_ij")
      return wrap(model.dP_ij(Never,Reclassified,0.0, 50.0));
    Rprintf("call_purged_ns: output_type not matched (\"%s\").\n",output_type.c_str());
    return wrap(model.negll()); // default
  }
  RcppExport
  SEXP call_purged_ps(SEXP sexp) {
    List args = as<List>(sexp);
    // bool debug = as<bool>(args("debug"));
    PurgedPS model(sexp,
		   as<int>(args("N")), 
		   as<int>(args("nbeta01")),
		   as<int>(args("nbeta12")),
		   as<int>(args("nages0")));
    std::string output_type = args("output_type");
    if (output_type == "negll")
      return wrap(model.negll());
    if (output_type == "negll_gradient")
      return wrap(model.negll_gradient());
    Rprintf("call_purged_ps: output_type not matched (%s)\n",output_type.c_str());
    return wrap(model.negll());
  }

}
