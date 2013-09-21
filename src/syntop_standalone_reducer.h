/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef SYNTOP_STANDALONE_REDUCER_INCLUDED
#define SYNTOP_STANDALONE_REDUCER_INLCUDED

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_lambert.h>

#include "variational_parameters.h"
#include "syntop_parameters.pb.h"
#include "gradient.h"
#include "util.h"
#include "vectorops.h"

#include  <stdint.h>  // <--- to prevent uint64_t errors!

#include <algorithm>
#include <limits>
#include <fstream>
#include <string>

#define MAX_DOC_ITER 2
#define ALWAYS_PLOT_DOC 0

using boost::scoped_ptr;

enum LikelihoodSideEffects { VALUE_ONLY,
                             COMPUTE_RHO_DENOMS,
                             COMPUTE_TAU_DENOMS,
                             COMPUTE_ALL_PARAMS };

const string kLhoodLabel[] = {"gamma prior", "gamma var",  // NOLINT
                               "phi|gamma", "phi|nu", "gamma/nu norm",
                              "doc voc", "phi var"};

/*
 * There are too many public functions; the data is concealed, but not
 * all of these should be exposed
 *
 */
class SyntopStandaloneReducer {
 public:
  VariationalParameters* vars_;
  // const SyntopParameters* params_;
  shared_ptr<const SyntopParameters> params_;

  // shared_ptr<gsl_matrix> nu_origin;

  SyntopStandaloneReducer() {}
  explicit SyntopStandaloneReducer(shared_ptr<const SyntopParameters> params);

  ~SyntopStandaloneReducer() {}

  StringMap reduce(string type, int index, StringMap* buffer);

  // void close();

  static double GlobalWeightTerm(const gsl_matrix* nu,
                                 const gsl_vector* beta,
                                 double alpha_trans,
                                 double alpha_top,
                                 bool finite);

 protected:
  int tau_coordinate_;
  int nu_coordinate_;

  int index;

  bool lhood_;
  bool check_order_;

  StringMap output;

  std::vector<std::string> temp_string_components_;

  void Optimize();

  /*
   * Based on the optimizations performed, write the relevant statistics to be
   * read in by the driver to create the new VariationalParameters object.
   */
  void Emit(StringMap* emission_counts);

  void ProcessKey(string key, double val);

  double GlobalWeightTermPrior(const gsl_vector* beta,
                               double alpha_top,
                               bool finite);

  double GlobalWeightTermSingle(int i,
                                const gsl_matrix* nu,
                                const gsl_vector* beta,
                                double alpha_trans);

  // double GlobalWeightTerm(const gsl_matrix* nu,
  //                        const gsl_vector* beta,
  //                        double alpha_trans,
  //                        double alpha_top,
  //                          bool finite);

  double GlobalVocabularyTerm(const gsl_matrix* tau,
                              const double vocab_sigma,
                              const int number_terms,
                              const int number_topics);

  double GlobalLikelihoodTerm();

  void UpdateNu();
  void LazyUpdateNu(int topic);
};

#endif
