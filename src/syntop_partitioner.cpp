/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef SYNTOPPARTITIONER_INCLUDED
#define SYNTOPPARTITIONER_INCLUDED

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/assert.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_lambert.h>

using lib_prob::safe_log;

/*
 * There are too many public functions; the data is concealed, but not
 * all of these should be exposed
 *
 */

#include  <stdint.h>  // <--- to prevent uint64_t errors!

#include <string>

#include "hadoop/Pipes.hh"
#include "hadoop/TemplateFactory.hh"
#include "hadoop/StringUtils.hh"

using std::string;

class SyntopPartitioner : public HadoopPipes::Partitioner {
 public:
  SyntopPartitioner(HadoopPipes::TaskContext& context) {} // NOLINT

  virtual int partition(const std::string& key, int numberOfReducers) {
    boost::split(temp_string_components_, key, boost::is_any_of("_"));
    // return ((temp_string_components_[0].length() +
    // temp_string_components_[1].length()) % numberOfReducers);

    return (temp_string_components_[0].length()
            + boost::lexical_cast<int>(temp_string_components_[1]))
            % numberOfReducers;
  }

 private:
  std::vector<std::string> temp_string_components_;
};

#endif
