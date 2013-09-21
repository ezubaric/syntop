/*
 * Copyright 2009 Jordan Boyd-Graber
 */

#include "syntop_standalone_mapper.h"
// This is included here to make life easier
#include "gradient.h"

#include <iostream>
#include <fstream>

/*
 * Begin implemented functions for local mode.
 */
SyntopStandaloneMapper::SyntopStandaloneMapper
(shared_ptr<const SyntopParameters> params) {
  params_ = params;
}

void SyntopStandaloneMapper::map(string line, StringMap* buffer) {
  StringMap count;

  vars_.reset(new VariationalParameters(*params_));
  doc_mapper = new DocumentMapper(vars_.get(), params_.get());

  doc = ParseDocument(line, -1);

  // shall we emit the doc_id together with the key?
  int doc_id = doc->get_id();

  doc->ResetPhi(params_->num_topics());
  doc_mapper->EStep(difference_tolerance, doc, doc_id);
  doc_mapper->Emit(doc, &count);

  for (int ii = 0; ii < params_->num_topics(); ++ii) {
    // Tau
    string key = "tau_" + SimpleItoa(ii) + "_-1";
    (*buffer)[key] += count[key];

    for (int vv = 0; vv < params_->vocab_size(); ++vv) {
      key = "tau_" + SimpleItoa(ii) + "_" + SimpleItoa(vv);
      (*buffer)[key] += count[key];
    }

    // Gamma
    // directly output to file
    key = "gamma_" + SimpleItoa(doc_id) + "_" + SimpleItoa(ii);
    (*buffer)[key] += count[key];

    // Nu ss
    for (int jj = 0; jj < params_->num_topics() + 1; ++jj) {
      // It's jj before ii because I wanted to put it inside the outer loop, but
      // the first index of nu ss needs to go from 0 to K instead of 0 to K-1
      string key_doc = "nu_" + SimpleItoa(jj) + "_" + SimpleItoa(ii) + "_doc";
      (*buffer)[key_doc] += count[key_doc];

      string key_trans = "nu_" + SimpleItoa(jj) + "_" +
        SimpleItoa(ii) + "_trans";
      (*buffer)[key_trans] += count[key_trans];
    }
  }

  // likelihood
  string key = "lhood_1_1";
  (*buffer)[key] += count[key];
  cout << "likelihood: " << key << "\t" << count[key] << endl;
}
