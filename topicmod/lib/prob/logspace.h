/*
 * logspace.h
 *
 * Copyright 2007, Jordan Boyd-Graber
 *
 * Functions for manipulating probabilities in log space
 */

#ifndef TOPICMOD_LIB_PROB_LOGSPACE_H__
#define TOPICMOD_LIB_PROB_LOGSPACE_H__
#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>

#include <cmath>
#include <limits>
#include <vector>
#include <iostream>

#ifndef isnan
# define isnan(x) \
  (sizeof(x) == sizeof(long double) ? isnan_ld(x) : sizeof(x) == sizeof(double) ? isnan_d(x) : isnan_f(x)) // NOLINT
static inline int isnan_f  (float       x) { return x != x; }
static inline int isnan_d  (double      x) { return x != x; }
static inline int isnan_ld (long double x) { return x != x; }
#endif

#ifndef isinf
# define isinf(x) \
  (sizeof(x) == sizeof(long double) ? isinf_ld(x) : sizeof(x) == sizeof(double) ? isinf_d(x) : isinf_f(x)) // NOLINT
static inline int isinf_f  (float       x) { return isnan (x - x); }
static inline int isinf_d  (double      x) { return isnan (x - x); }
static inline int isinf_ld (long double x) { return isnan (x - x); }
#endif

namespace lib_prob {

using std::cout;
using std::endl;

inline double next_uniform() {
  return (double)rand() / ((double)RAND_MAX + 1.0);
}

double gsl_blas_dsum(const gsl_vector* v);

double normalize(gsl_vector* v);

double safe_log(double x);

// Given log(a) and log(b), return log(a + b).
double log_sum(double log_a, double log_b);

// Given log(a) and log(b), return log(a - b).
double log_diff(double log_a, double log_b);

/*
 * returns the element randomly sampled from the log
 * probabilities in array (number is the number of elements)
 */
int log_vector_sample(std::vector<double> vals, int length);
int log_sample(double* vals, int length);

void gsl_matrix_print(const gsl_matrix* x);
void gsl_vector_print(const gsl_vector* x);

// Given a log vector, log a_i, compute log sum a_i.  Returns the sum.
double log_normalize(gsl_vector* x);

// Compute the log sum over all elements in the vector
double log_sum(const gsl_vector* x);

double log_dirichlet_likelihood(const double sum,
                                const double prior_sum,
                                const std::vector<int>& counts,
                                bool debug = false);

double log_dirichlet_likelihood(const double sum,
                                const double prior_scale,
                                const gsl_vector* prior,
                                const std::vector<int>& counts);

}  // namespace lib_prob

#endif  // TOPICMOD_LIB_PROB_LOGSPACE_H__
