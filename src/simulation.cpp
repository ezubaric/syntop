/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#ifndef SIMULATION_INCLUDED
#define SIMULATION_INCLUDED

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/assert.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_lambert.h>

#include "util.h"
#include "vectorops.h"
#include "variational_parameters.h"
#include "gradient.h"
#include "document_mapper.h"
#include "lens_doc.h"
#include "syntop_standalone_mapper.h"
#include "syntop_standalone_reducer.h"

#define MAX_SENT 100
#define MAX_DOC_ITER 2
#define ALWAYS_PLOT_DOC 0

using lib_prob::safe_log;

#include <stdint.h>  // <--- to prevent uint64_t errors!

#include <algorithm>
#include <limits>
#include <fstream>
#include <string>

using std::string;
using std::ifstream;

DEFINE_int(truncation, 10, "Truncation level");
DEFINE_int(vocab_size, 49, "How many types in corpus");
DEFINE_int(doc_limit, -1, "How many documents we read");
DEFINE_string(directory, "", "Where we read the data");
DEFINE_string(dataset, "", "What dataset we are going to use");
DEFINE_string(suffix, "", "What suffix we are goint to append");
DEFINE_bool(ignore_doc, false, "Ignore documents");
DEFINE_bool(finite, false, "Finite model");
DEFINE_bool(ignore_trans, false, "Ignore syntax");
DEFINE_bool(shortcut_gsl, false, "Shortcut GSL with exact updates");

int main(int argc, char **argv) {
  InitFlags(argc, argv);

  SyntopParameters* params = new SyntopParameters();

  params->set_vocab_size(FLAGS_vocab_size);
  params->set_alpha_trans(1.0);
  params->set_alpha_doc(1.0);
  params->set_alpha_top(100.0);
  params->set_num_docs(FLAGS_doc_limit);
  params->set_num_topics(FLAGS_truncation);
  params->set_model_name(FLAGS_directory + "/" + FLAGS_dataset
                         + "_params_" + FLAGS_suffix + "/current");
  // params->set_model_name(FLAGS_directory + "_params/current");
  params->set_finite(FLAGS_finite);
  params->set_shortcut_gsl(FLAGS_shortcut_gsl);
  params->set_ignore_docs(FLAGS_ignore_doc);
  params->set_ignore_trans(FLAGS_ignore_trans);

  shared_ptr<const SyntopParameters> param_container(params);

  scoped_ptr<SyntopStandaloneMapper> mapper
    (new SyntopStandaloneMapper(param_container));

  string line;
  ifstream datafile((FLAGS_dataset + "/" + FLAGS_dataset + ".0.prs").c_str());
  
  StringMap buffer;

  if(true){
    cout << "update beta phase started..." << endl;

    shared_ptr<const SyntopParameters> params_ = param_container;
    scoped_ptr<VariationalParameters> vars_;
    vars_.reset(new VariationalParameters(*params_));

    if (!params_->finite()) {
      optimize_beta(vars_->beta_.get(), vars_->gamma_.get(),
                    vars_->nu_.get(), params_->ignore_docs(),
                    params_->ignore_trans(), params_->alpha_top(),
                    params_->alpha_doc(), params_->alpha_trans(),
                    params_->model_name());

      ofstream yourfile;
      yourfile.open((FLAGS_directory + "/" + FLAGS_dataset + "_params_"
                     + FLAGS_suffix + "/current.beta").c_str());
      // yourfile.open((FLAGS_directory + "_params/tmp.beta").c_str());
      for (int i = 0; i < params_->num_topics(); i++) {
	yourfile << gsl_vector_get(vars_->beta_.get(), i) << "\t";
      }
      yourfile.close();
    }

    cout << "update beta finished..." << endl;
  }

  cout << "map phase started..." << endl;

  if (datafile.is_open()) {
    while (datafile.good()) {
      getline(datafile, line);
      if (line.length() > 1) {
        mapper->map(line, &buffer);
      }
    }

    datafile.close();
  } else {
    cout << "Unable to open file" << endl;
  }

  cout << "map phase finished..." << endl;

  ofstream tokenfile;
  tokenfile.open((FLAGS_directory + "/" + FLAGS_dataset + "_params_"
                  + FLAGS_suffix + "/tmp-tokens.txt").c_str());
  StringMap::const_iterator last = buffer.end();
  for (StringMap::const_iterator itr = buffer.begin(); itr != last; itr++) {
    tokenfile << itr->first << "\t" << itr->second << endl;
  }
  tokenfile.close();

  scoped_ptr<SyntopStandaloneReducer> reducer;

  cout << "reduce tau phase started..." << endl;
  string type = "tau";
  for (int i = 0; i < FLAGS_truncation; i++) {
    reducer.reset(new SyntopStandaloneReducer(param_container));

    StringMap output = reducer->reduce(type, i, &buffer);
    ofstream paramfile;
    paramfile.open((FLAGS_directory + "/" + FLAGS_dataset + "_params_"
                    + FLAGS_suffix + "/tmp-" + type + "-"
                   + SimpleItoa(i) + ".txt").c_str());
    StringMap::const_iterator last = (output).end();
    for (StringMap::const_iterator itr = (output).begin();
         itr != last; itr++) {
      paramfile << itr->first << "\t" << itr->second << endl;
    }
    paramfile.close();
    cout << "reduce " + type + " " + SimpleItoa(i) + " finished..." << endl;
  }

  cout << "reduce nu phase started..." << endl;
  double likelihood = 0.0;
  type = "nu";
  for (int i = 0; i < FLAGS_truncation + 1; i++) {
    reducer.reset(new SyntopStandaloneReducer(param_container));
    StringMap output = reducer->reduce(type, i, &buffer);

    ofstream yourfile;
    yourfile.open((FLAGS_directory + "/" + FLAGS_dataset + "_params_"
                   + FLAGS_suffix + "/tmp-" + type + "-"
                   + SimpleItoa(i) + ".txt").c_str());
    StringMap::const_iterator last = (output).end();
    for (StringMap::const_iterator itr = (output).begin();
         itr != last; itr++) {
      if (boost::starts_with(itr->first, "lhood")) {
        cout << "lhood_1_1 detected..." << itr->first << "\t" <<
          itr->second << endl;
        likelihood += itr->second;
      } else {
        yourfile << itr->first << "\t" << itr->second << endl;
      }
    }
    yourfile.close();
    cout << "reduce " + type + " " + SimpleItoa(i) + " finished..." << endl;
  }

  {
    type = "lhood";
    reducer.reset(new SyntopStandaloneReducer(param_container));
    StringMap output = reducer->reduce(type, 1, &buffer);

    ofstream yourfile;
    yourfile.open((FLAGS_directory + "/" + FLAGS_dataset + "_params_"
                   + FLAGS_suffix + "/tmp-" + type + ".txt").c_str());
    StringMap::const_iterator last = (output).end();
    for (StringMap::const_iterator itr = (output).begin(); itr != last; itr++) {
      if (boost::starts_with(itr->first, type)) {
        yourfile << itr->first << "\t" << (itr->second + likelihood) << endl;
      }
    }
    yourfile.close();
    cout << "reduce " + type + " finished..." << endl;
  }

  {
    cout << "reduce gamma phase started..." << endl;
    type = "gamma";
    ofstream yourfile;
    yourfile.open((FLAGS_directory + "/" + FLAGS_dataset + "_params_"
                   + FLAGS_suffix + "/tmp-" + type + ".txt").c_str());
    StringMap::const_iterator last_gamma = (buffer).end();
    for (StringMap::const_iterator itr_gamma=(buffer).begin();
         itr_gamma != last_gamma; itr_gamma++) {
      // cout << itr_gamma->first << "\ttest\t" << itr_gamma->second << endl;
      if (boost::starts_with(itr_gamma->first, type)) {
	// cout << itr_gamma->first << "\ttest\t" << itr_gamma->second << endl;
        yourfile << itr_gamma->first << "\t" << itr_gamma->second << endl;
      }
    }
    yourfile.close();
    cout << "reduce " + type + " finished..." << endl;
  }

  return 0;
}

#endif
