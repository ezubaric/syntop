/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef __MATH_VECTOROPS_INCLUDED
#define __MATH_VECTOROPS_INCLUDED

#include <cmath>
#include <limits>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_psi.h>
// #include "specialfunc.h"

#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

#ifndef isnan
# define isnan(x) ((x) != (x))
#endif

/*
 * take the exponent of a vector
 *
 * If the exponent is infinite, then we replace the value with a
 * suitably large max_val
 */
void vexp(const gsl_vector* v,
          gsl_vector* exp_v,
          double max_val = std::numeric_limits<double>::infinity());

/* take the exponent of a matrix */
void mexp(const gsl_matrix* m, gsl_matrix* exp_m);

/* like vexp except that it also computes sum x log x */
double vexp_entropy(const gsl_vector* v, gsl_vector* exp_v);

double ventropy(const gsl_vector* v);

double lgamma(double x);

void mlog(const gsl_matrix* m, gsl_matrix* log_m);

void vlog(const gsl_vector* v, gsl_vector* log_v);

void vlogit(const gsl_vector* v, gsl_vector* log_v);

void vsigmoid(const gsl_vector* v, gsl_vector* sig_v);

double vlog_entropy(const gsl_vector* v, gsl_vector* log_v);

double entropy(const gsl_vector* v);

void vdigamma(const gsl_vector* v, gsl_vector* digamma_v);

void vlgamma(const gsl_vector* v, gsl_vector* lgamma_v);

double gsl_blas_dsum(const gsl_vector* v);

double gsl_blas_dsum(const gsl_matrix* v);

double gsl_matrix_rowsum(const gsl_matrix* m, const int row);

double dot_product(const gsl_vector* a, const gsl_vector* b);

void uniform(gsl_vector* v);

/*
  This function takes as input a multinomial parameter vector and
  computes the "total" variance, i.e., the sum of the diagonal of the
  covariance matrix.

  If the multinomial parameter is unnormalized, then the variance of
  the normalized multinomial vector will be computed and then
  multiplied by the scale of the vector.
 */
double MultinomialTotalVariance(const gsl_vector* v);
/*
  Computes covariance using the renormalization above and adds it to
  an existing matrix.
*/
void MultinomialCovariance(double alpha,
                           const gsl_vector* v,
                           gsl_matrix* m);

double MatrixProductSum(const gsl_matrix* m1,
                        const gsl_matrix* m2);

double MatrixProductProductSum(const gsl_matrix* m1,
                               const gsl_matrix* m2,
                               const gsl_matrix* m3);

double SumLGamma(const gsl_vector* v);

void mtx_fprintf(const char* filename, const gsl_matrix * m);

void mtx_fscanf(const char* filename, gsl_matrix * m);

double vct_accum(const int i, const double contribution, gsl_vector* v);

double mtx_accum(const int i,
                 const int j,
                 const double contribution,
                 gsl_matrix* m);

void vct_fscanf(const char* filename, gsl_vector* v);

void vct_fprintf(const char* filename, gsl_vector* v);

#endif
