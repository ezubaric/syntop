/*
 * Copyright 2009 Jordan Boyd-Graber
 * Date: March 2008
 *
 * This file was spun off from logspace.h in order to create a
 * logspace file that wouldn't depend on gsl.
 */
#ifndef __MATH_LOGSPACE_BASE_H__
#define __MATH_LOGSPACE_BASE_H__

#include <cstdlib>
#include <iostream>
#include <limits>
#include <assert.h>
#include <cmath>
#include <vector>

#ifndef isnan
# define isnan(x) \
  (sizeof(x) == sizeof(long double) ? isnan_ld(x) \ // NOLINT
: sizeof(x) == sizeof(double) ? isnan_d(x) \ // NOLINT
  : isnan_f(x))
static inline int isnan_f(float       x) { return x != x; }
static inline int isnan_d(double      x) { return x != x; }
static inline int isnan_ld(long double x) { return x != x; }
#endif

#ifndef isinf
# define isinf(x) \
  (sizeof(x) == sizeof(long double) ? isinf_ld(x) \ // NOLINT
    : sizeof(x) == sizeof(double) ? isinf_d(x) \ // NOLINT
  : isinf_f(x))
static inline int isinf_f(float       x) { return isnan (x - x); }
static inline int isinf_d(double      x) { return isnan (x - x); }
static inline int isinf_ld(long double x) { return isnan (x - x); }
#endif

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

/*
 * Stupid "sampling" function for deterministic testing (i.e. in unit tests)
 */
int sample_first_nonzero(double* vals, int length);
int sample_max(double* vals);

bool is_nan(double val);

#endif  // __MATH_LOGSPACE_BASE_H__
