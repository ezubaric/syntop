/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef LENSMODEL_INCLUDED
#define LENSMODEL_INCLUDED

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/assert.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_lambert.h>

#include "util.h"
#include "specialfunc.h"
#include "vectorops.h"

#define MAX_DOC_ITER 2
#define ALWAYS_PLOT_DOC 0

enum LikelihoodSideEffects { VALUE_ONLY,
                             COMPUTE_RHO_DENOMS,
                             COMPUTE_TAU_DENOMS,
                             COMPUTE_ALL_PARAMS };

const char kLhoodLabel[] = {"gamma prior", "gamma var", "phi|gamma",
                            "phi|nu", "gamma/nu norm", "doc voc",
                            "phi var"};

/*
 * There are too many public functions; the data is concealed, but not
 * all of these should be exposed
 *
 */
class GlobalOptimizer {
 public:
  GlobalOptimizer(VariationalParameters* vars,
                  const SyntopParameters* params);

  /*
   * Merge takes in output from (key, value) pairs
   *
   * As it sees keys, it updates the boolean toggles that tell it which
   * variables to optimize.
   */
  Merge(const string& key, double val);

  /*
   * Based on observed keys, update the relevant variables
   */
  Update();

  /*
   * Based on the optimizations performed, write the relevant statistics to be
   * read in by the driver to create the new VariationalParameters object.
   */
  Emit(StringMap* emission_counts);

 private:
  int tau_coordinate_;
  int nu_coordinate_;
  bool lhood_;
  bool check_order_;

  VariationalParameters* vars_;
  const SyntopParameters* params_;

  std::vector<std::string> temp_string_components_;

  GlobalOptimizer();  // No evil constructors
};

#endif
