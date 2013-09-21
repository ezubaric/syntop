/*
 * logspace.cpp
 *
 * Copyright 2007, Jordan Boyd-Graber
 *
 * For computing operations on logspace double (usually probabilitys)
 *
 */

#include <limits>
#include <vector>

#include "topicmod/lib/prob/logspace.h"

namespace lib_prob {

  /*
   * This does not work in logspace; perhaps it should be moved
   */
double gsl_blas_dsum(const gsl_vector* v) {
  double sum = 0;
  for (unsigned int ii = 0; ii < v->size; ++ii) {
    sum += gsl_vector_get(v, ii);
  }
  return sum;
}

double normalize(gsl_vector* v) {
  double sum = gsl_blas_dsum(v);
  gsl_blas_dscal(1 / sum, v);
  return sum;
}

// Given log(a) and log(b), return log(a - b).
double log_diff(double log_a, double log_b) {
  double val;
  double dangerous_part = exp(log_b - log_a);
  assert(dangerous_part < 1.0);
  val = log_a + log(1.0 - dangerous_part);
  return val;
}

double safe_log(double x) {
  if (x <= 0) {
    return(-1e16);
  } else {
    return(log(x));
  }
}

// Given log(a) and log(b), return log(a + b).
double log_sum(double log_a, double log_b) {
  double v;

  if (log_a == -std::numeric_limits<double>::infinity() &&
      log_b == log_a) {
    return -std::numeric_limits<double>::infinity();
  } else if (log_a < log_b) {
    v = log_b + log(1 + exp(log_a - log_b));
  } else {
      v = log_a + log(1 + exp(log_b - log_a));
  }
  return(v);
}

int log_vector_sample(std::vector<double> vals, int length) {
  double normalizer = safe_log(0.0);
  int ii = 0;
  assert(length > 0 && length <= (int)vals.size());
  for (ii = 0; ii < length; ++ii) {
    normalizer = log_sum(normalizer, vals[ii]);
  }

  double val = 0, sum = 0, cutoff = (double)rand() / ((double)RAND_MAX + 1.0);
  for (ii = 0; ii < length; ++ii) {
    val = exp(vals[ii] - normalizer);
    sum += val;
    if (sum >= cutoff)
      break;
  }
  assert(ii < length);
  return ii;
}

/*
 * returns the element randomly sampled from the log
 * probabilities in array (number is the number of elements)
 */
int log_sample(double* vals, int length) {
  double normalizer = safe_log(0.0);
  int ii;
  for (ii = 0; ii < length; ++ii) {
    normalizer = log_sum(normalizer, vals[ii]);
  }

  double val = 0, sum = 0, cutoff = (double)rand() / ((double)RAND_MAX + 1.0);
  for (ii = 0; ii < length; ++ii) {
    val = exp(vals[ii] - normalizer);
    sum += val;
    if (sum >= cutoff)
      break;
  }
  assert(ii < length);
  return ii;
}

double log_dirichlet_likelihood(const double sum,
                                const double prior_sum,
                                const std::vector<int>& counts,
                                bool debug) {
  double val = 0.0;
  int length = counts.size();

  double prior_value = prior_sum / (double)length;
  val += gsl_sf_lngamma(prior_sum);
  val -= (double)length * gsl_sf_lngamma(prior_value);

  if (debug) cout << "Likelihood (" << sum << "," << prior_sum << "," <<
    prior_value << "," << length << ") = " << val << endl;

  for (int ii = 0; ii < length; ++ii) {
    if (debug) cout << "\tGAMMA(" << prior_value << " + " <<
      (double)counts[ii] << " = " << prior_value +
      (double)counts[ii] <<  ") -> " << val << endl;
    val += gsl_sf_lngamma(prior_value + (double)counts[ii]);
  }
  val -= gsl_sf_lngamma(prior_sum + sum);

  if (debug) cout << endl;

  return val;
}

double log_dirichlet_likelihood(const double sum,
                                const double prior_scale,
                                const gsl_vector* prior,
                                const std::vector<int>& counts) {
  double val = 0.0;
  int length = counts.size();

  val += gsl_sf_lngamma(prior_scale);
  for (int ii = 0; ii < length; ++ii) {
    double prior_value = gsl_vector_get(prior, ii);
    val -= gsl_sf_lngamma(prior_value);
    val += gsl_sf_lngamma(prior_value + (double)counts[ii]);
  }
  val -= gsl_sf_lngamma(prior_scale + sum);

  return val;
}

void gsl_matrix_print(const gsl_matrix* x) {
  for (int ii = 0; ii < (int)x->size1; ++ii) {
    for (int jj = 0; jj < (int)x->size2; ++jj) {
      cout << " " << gsl_matrix_get(x, ii, jj);
    }
    cout << endl;
  }
}

void gsl_vector_print(const gsl_vector* x) {
  for (int ii = 0; ii < (int)x->size; ++ii) {
    cout << gsl_vector_get(x, ii) << " ";
  }
  cout << endl;
}

double log_sum(const gsl_vector* x) {
  double sum = gsl_vector_get(x, 0);

  for (unsigned int ii = 1; ii < x->size; ii++) {
    sum = log_sum(sum, gsl_vector_get(x, ii));
  }
  return sum;
}

// Given a log vector, log a_i, compute log sum a_i.  Returns the sum.
double log_normalize(gsl_vector* x) {
  double sum = gsl_vector_get(x, 0);
  unsigned int i;

  for (i = 1; i < x->size; i++) {
    sum = log_sum(sum, gsl_vector_get(x, i));
  }

  for (i = 0; i < x->size; i++) {
    double val = gsl_vector_get(x, i);
    gsl_vector_set(x, i, val - sum);
  }
  return sum;
}
}  // namespace lib_prob_util
