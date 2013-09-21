/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef __MATH_LOGSPACE_H__
#define __MATH_LOGSPACE_H__
#include <cmath>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_gamma.h>

// Given two log vectors, log a_i and log b_i, compute
// log sum (a_i * b_i).
double log_dot_product(const gsl_vector* log_a, const gsl_vector* log_b);

// Given log(a) and log(b), return log(a - b).
double log_diff(double log_a, double log_b);

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

#endif  // __MATH_LOGSPACE_H__
