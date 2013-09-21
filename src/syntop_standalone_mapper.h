/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef SYNTOPMAPPER_INCLUDED
#define SYNTOPMAPPER_INCLUDED

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/assert.hpp>
#include <boost/lexical_cast.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_lambert.h>

#include "util.h"
#include "vectorops.h"
#include "variational_parameters.h"
#include "gradient.h"
#include "document_mapper.h"
#include "lens_doc.h"

#define MAX_SENT 100
#define MAX_DOC_ITER 2
#define ALWAYS_PLOT_DOC 0

using lib_prob::safe_log;

/*
 * There are too many public functions; the data is concealed, but not
 * all of these should be exposed
 *
 */

#include  <stdint.h>  // <--- to prevent uint64_t errors!

#include <algorithm>
#include <limits>
#include <fstream>
#include <string>

using std::string;
using boost::scoped_ptr;

class SyntopStandaloneMapper {
 public:
  static const float difference_tolerance = 0.001;

  scoped_ptr<VariationalParameters> vars_;
  shared_ptr<const SyntopParameters> params_;
  DocumentMapper* doc_mapper;

  explicit SyntopStandaloneMapper(shared_ptr<const SyntopParameters> params);
  SyntopStandaloneMapper() {}

  ~SyntopStandaloneMapper() {}

  void map(string line, StringMap* buffer);

 protected:
  Document* doc;
  string line;

 private:
  SyntopStandaloneMapper(const SyntopStandaloneMapper&) {}
};

#endif
